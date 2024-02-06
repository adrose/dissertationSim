## Load library(s)
library(markovchain)
library(msm)
library(brms)
library(foreach)
library(doParallel)
## Create functions
# simulate discrete Markov chains according to transition matrix P
## Here P will be a transition matrix -- where 
run.mc.sim <- function( P, num.iters = 50 ) {
  
  # number of possible states
  num.states <- nrow(P)
  
  # stores the states X_t through time
  states     <- numeric(num.iters)
  
  # initialize variable for first state 
  states[1]    <- sample(1:num.states, size = 1)
  
  for(t in 2:num.iters) {
    
    # probability vector to simulate next state X_{t+1}
    p  <- P[states[t-1], ]
    
    ## draw from multinomial and determine state
    states[t] <-  which(rmultinom(1, 1, p) == 1)
  }
  return(states)
}

## Now we need to get the coefficients from three levels of transition matrices
## Create the stationary matrix her 
P.stat <- matrix(c( 0.9, 0.05, 0.05, 
                    0.05, 0.9,   0.05, 
                    .05,   .05,   .9  ), nrow=3, ncol=3, byrow = TRUE)
## Create the pseudo here
P.tran <- matrix(c( 0.6, 0.2, 0.2, 
                    0.2, 0.6,   0.2, 
                    .2,   .2,   .6  ), nrow=3, ncol=3, byrow = TRUE)
## Now create the random walk here
P.rand <- matrix(c( 0.34, 0.33, 0.33, 
                    0.33, 0.34,   0.33, 
                    .33,   .33,   .34  ), nrow=3, ncol=3, byrow = TRUE)
## Now loop through each of these and grab our model coefficient estimates
trans.mats <- list(P.stat, P.tran, P.rand)
out.models <- list()
for(i in 1:length(trans.mats)){
  ## First declare our sample characteristics
  n <- 30
  p <- 100
  ## Now create a data frame to store all of the data
  pop1.vals <- data.frame(ID = paste(rep(1:n, each=p), "A", sep="_"), time=rep(1:p), 
                          state=NA, laggedState = NA)
  ## Now simulate the data from the known markov trnaisiton values
  for(N in 1:n){
    ## Identify sample index
    index.vals <- which(pop1.vals$ID==paste(N, "A", sep="_"))
    ## Now figure out a way to add variability to P1 somehow...
    state.vals <- run.mc.sim(trans.mats[[i]], p+1)
    lagged.vals <- lag(state.vals)[1:p]
    pop1.vals[index.vals,"state"] <- state.vals[-1]
    pop1.vals[index.vals,"laggedState"] <- lagged.vals
  }
  ## Now train the model
  pop1.vals$laggedState <- factor(pop1.vals$laggedState)
  pop1.vals$state <- factor(pop1.vals$state)
  cmodFE <- brm(state ~ 0 + Intercept + laggedState, data = pop1.vals, family=categorical(), 
                cores = 6, chains = 6, iter = 4000)
  ## Now store these values
  out.models[[i]] <- cmodFE
  
}
## Now check our predicted values
predict(out.models[[1]], newdata = data.frame(laggedState = c(1,1,1,2,2,2,3,3,3), pop="B", ID = "2_C"), allow_new_levels=FALSE, summary=TRUE)
predict(out.models[[2]], newdata = data.frame(laggedState = c(1,1,1,2,2,2,3,3,3), pop="B", ID = "2_C"), allow_new_levels=FALSE, summary=TRUE)
predict(out.models[[3]], newdata = data.frame(laggedState = c(1,1,1,2,2,2,3,3,3), pop="B", ID = "2_C"), allow_new_levels=FALSE, summary=TRUE)

## Now create these models WITH partiicpant level variability in the intercept
priorsMLString <- c(
  ## All mu2 values here --> intercept for 2 --> 3 values
  prior("constant(A)", coef="Intercept", class="b", dpar="mu2"),
  prior("constant(B)", coef="laggedState2", class="b", dpar="mu2"),
  prior("constant(C)", coef="laggedState3", class="b", dpar="mu2"),
  prior("constant(G)", coef="", class="sd", dpar="mu2"), ## random effect variance here
  
  ## All mu3 values here --> intercept for 2 --> 3 values
  prior("constant(D)", coef="Intercept", class="b", dpar="mu3"),
  prior("constant(E)", coef="laggedState2", class="b", dpar="mu3"),
  prior("constant(F)", coef="laggedState3", class="b", dpar="mu3"),
  prior("constant(H)", coef="", class="sd", dpar="mu3") ## random effect variance here
)
## now isolate all of the fixed effects from our pop models
fixef(out.models[[1]])
ranef.var <- c(.1, 1, 5)
out.priors <- list()
out.priors.index <- 1
for(r in 1:length(ranef.var)){
  for(i in 1:length(out.models)){
    ## First change the fixef variables
    newprior <- priorsMLString
    fixefVals <- fixef(out.models[[i]])
    for(j in 1:8){
      if(j < 7){
        newprior$prior <- gsub(pattern = LETTERS[j], replacement = fixefVals[j,"Estimate"], x = newprior$prior)
      }
      if(j > 6){
        newprior$prior <- gsub(pattern = LETTERS[j], replacement = ranef.var[r], x = newprior$prior)
      }
    }
    out.priors[[out.priors.index]] <- newprior
    out.priors.index <- out.priors.index + 1
  }
}

## Now estimate our models
mod.list <- lapply(out.priors, function(x) sampGen <- brm(state ~ 0 + Intercept + laggedState + (1|ID), data = pop1.vals, family=categorical(), 
                                                          cores = 4,sample_prior = "only", prior = x))
## Now create our samples
samp.gen <- lapply(mod.list, posterior_predict)

## Now fit our models to these data
## Lets start with 10 models at a time
out.mods2 <- list()
list.vals <- 1
mod.count <- 1
## Create a do parallel loop here
cl <- makeCluster(6) #not to overload your computer
registerDoParallel(cl)
for(list.vals in 1:length(samp.gen)){
  tmp.out <- foreach(i = 1:6, .packages = c("brms")) %dopar%{
    ## First grab the values from the data
    tmp.data <- samp.gen[[list.vals]][i,]
    ## Now attach these to our data
    pop1.vals$newData <- tmp.data
    pop1.vals$newLag <- NA
    ## Now create our lagged variables within each participant
    ## Now go through every unique ID and make the lagged variable within them
    for(l in unique(pop1.vals$ID)){
      index <- which(pop1.vals$ID == l)
      new.vals <- c(NA, lag(pop1.vals$newData[index])[1:length(index)-1])
      pop1.vals$newLag[index] <- new.vals
      pop1.vals$newLag <- factor(pop1.vals$newLag, levels = c(1,2,3))
      pop1.vals$newData <- factor(pop1.vals$newData, levels = c(1,2,3))
    }
    ## Now estimate our model
    mod.tmp <- brm(newData ~ 0 + Intercept + newLag + (1|ID), data = pop1.vals, 
                   family=categorical(),cores = 1, control = list(adapt_delta = 0.9, max_treedepth = 12), chains = 3)
    out <- summary(mod.tmp)
    rm(mod.tmp)
    out
  }
  out.mods2[[mod.count]] <- tmp.out
  mod.count <- mod.count + 1
}
doParallel::stopImplicitCluster()
out.name <- paste("./data/sample_", length(unique(pop1.vals$ID)), "_obs_", length(unique(pop1.vals$time)), ".RDS", sep='')
saveRDS(out.mods2, file = out.name)
