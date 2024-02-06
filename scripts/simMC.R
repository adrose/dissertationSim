## This firt script is going to be used to examine differeinces when estimating a markov model
## ignorninig an individual's propensity to transition
## so I will simulate at least from two populations:
## one where the population

## Load library(s)
library(markovchain)
library(msm)
library(brms)
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

## Now create a function which will sample n values and the sum of these values will be 1
n <- 4
sum(diff(c(0, sort(runif(n)), 1)))
retRow <- function(n){
  out.row <- diff(c(0, sort(runif(n)), 1))
  return(out.row)
}
retRow(n)

## Now create a function which will return a population transition matrix, and also a series of subject specific matrices
## The inputs to the function will include: a row which will 
retMat <- function(fixedTrans = retRow(2), randTrans = c("r_1", "r_1", "r_1"))



# setup transition matrix 
P1 <- matrix(c( 0.7, 0.15, 0.15, 
                 0.15, 0.7,   0.15, 
                 .15,   .15,   .7  ), nrow=3, ncol=3, byrow = TRUE)
P2 <- t(matrix(c( 0.2, 0.2, 0.6, 
                  0.1, 0.3,   0.6, 
                  0.6,   0.2,   0.1  ), nrow=3, ncol=3))
n.state <- 2
#P1 <- matrix(c(retRow(n.state), retRow(n.state)), nrow=n.state+1, ncol=n.state+1, byrow = TRUE)
#P2 <- matrix(c(retRow(n.state), retRow(n.state), retRow(n.state)), nrow=n.state+1, ncol=n.state+1, byrow = TRUE)
P2 <- P1

## Now create 50 participants from each matrix
run.mc.sim(P1, 100)
sampSize = 25
iterLength = 50
pop1.vals <- data.frame(ID = paste(rep(1:sampSize, each=iterLength), "A", sep="_"), time=rep(1:iterLength), 
                        state=NA, laggedState = NA,pop ="A")
pop2.vals <- data.frame(ID = paste(rep(1:sampSize, each=iterLength), "B", sep="_"), time=rep(1:iterLength), 
                        state=NA, laggedState = NA, pop ="B")
set.seed(16)
for(i in 1:sampSize){
  ## Identify sample index
  index.vals <- which(pop1.vals$ID==paste(i, "A", sep="_"))
  ## Now figure out a way to add variability to P1 somehow...
  state.vals <- run.mc.sim(P1, iterLength+1)
  lagged.vals <- lag(state.vals)[1:iterLength]
  pop1.vals[index.vals,"state"] <- state.vals[-1]
  pop1.vals[index.vals,"laggedState"] <- lagged.vals
  ## Now figure out a way to add variability to P2 somehow...
  state.vals <- run.mc.sim(P2, iterLength+1)
  lagged.vals <- lag(state.vals)[1:iterLength]
  pop2.vals[index.vals,"state"] <- state.vals[-1]
  pop2.vals[index.vals,"laggedState"] <- lagged.vals
}

## Now combine these
pop.vals <- dplyr::bind_rows(pop1.vals, pop2.vals)

## Now train msm model?
## Now make a q matrix
pq <- matrix(c(retRow(2),retRow(2),retRow(2)), nrow=3, ncol=3, byrow = TRUE)
mod <- msm(formula = state ~ time, data = pop.vals, qmatrix = pq, subject = ID, gen.inits = FALSE, 
           control=list(fnscale=4000, maxit=200000), covariates = ~ pop)
pmatrix.msm(mod, covariates = list(pop="A"), t=1)
pmatrix.msm(mod, covariates = list(pop="B"), t=1)

## Now try to run a BRM model using these values
pop.vals$state <- factor(pop.vals$state)
pop.vals$laggedState <- factor(pop.vals$laggedState)
#rand.markov.mod.est <- brm(state ~ laggedState + (1|ID), data = pop.vals, family=categorical(), cores = 4)
correct.mod <- brm(state ~ -1 + laggedState + pop + (1|ID), data = pop.vals, family=categorical(), cores = 4)

out.predB <- predict(correct.mod, newdata = data.frame(laggedState = c(1,2,3), pop="B", ID = "1_B"), allow_new_levels=TRUE)
out.predBA <- pmatrix.msm(mod, covariates = list(pop="B"), t = 1)
out.predA <- predict(correct.mod, newdata = data.frame(laggedState = c(1,2,3), pop="A", ID = "2_C"), allow_new_levels=TRUE)
out.predAA <- pmatrix.msm(mod, covariates = list(pop="B"), t = 1)

## Now try to generate data using brms
r1 <- get_prior(state ~ -1 + laggedState + pop + (1|ID), data = pop.vals, family=categorical(),sample_prior = "only")
priors1 <- c(
  ## All mu2 values here --> intercept for 2 --> 3 values
  #prior("constant(0)", class="b", dpar="mu2"),
  prior("constant(0)", coef="laggedState1", class="b", dpar="mu2"),
  prior("constant(10.06)", coef="laggedState2", class="b", dpar="mu2"),
  prior("constant(0)", coef="laggedState3", class="b", dpar="mu2"),
  prior("constant(0.04)", coef="popB", class="b", dpar="mu2"),
  prior("constant(0)", coef="", class="sd", dpar="mu2"), ## random effect variance here
  #prior("constant(.1)", coef="", class="sd", dpar="mu2", group="ID"), ## random effect variance here
  #prior("constant(.1)", coef="Intercept", class="sd", dpar="mu2", group="ID"), ## random effect variance here

  ## All mu3 values here --> intercept for 2 --> 3 values
  #prior("constant(0)", class="b", dpar="mu3"),
  prior("constant(0)", coef="laggedState1", class="b", dpar="mu3"),
  prior("constant(10.52)", coef="laggedState2", class="b", dpar="mu3"),
  prior("constant(3.06)", coef="laggedState3", class="b", dpar="mu3"),
  prior("constant(0.04)", coef="popB", class="b", dpar="mu3"),
  prior("constant(0)", coef="", class="sd", dpar="mu3") ## random effect variance here
  #prior("constant(.1)", coef="", class="sd", dpar="mu3", group="ID"), ## random effect variance here
  #prior("constant(.1)", coef="Intercept", class="sd", dpar="mu3", group="ID") ## random effect variance here
)
sampGen1 <- brm(state ~ -1 + laggedState + pop + (1|ID), data = pop.vals, family=categorical(), 
                sample_prior = "only", prior = priors1)
predVals <- predict(sampGen)
newY <- predict(sampGen1, newdata = data.frame(laggedState = c(1,2,3), pop="B", ID = c("1_B", "2_B", "3_B")), allow_new_levels=TRUE)
## Now use posterioer predict to get stat evalues
postPredVals <- posterior_predict(sampGen1)
test.dat <- cbind(predVals, pop.vals)

## Now try to change these constants and see where that gets us
priors2 <- c(
  ## All mu2 values here --> intercept for 2 --> 3 values
  #prior("constant(0)", class="b", dpar="mu2"),
  prior("constant(-1.51)", coef="Intercept", class="b", dpar="mu2"),
  prior("constant(3.02)", coef="laggedState2", class="b", dpar="mu2"),
  prior("constant(1.47)", coef="laggedState3", class="b", dpar="mu2"),
  prior("constant(0)", coef="popB", class="b", dpar="mu2"),
  prior("constant(1)", coef="", class="sd", dpar="mu2"), ## random effect variance here
  #prior("constant(.1)", coef="", class="sd", dpar="mu2", group="ID"), ## random effect variance here
  #prior("constant(.1)", coef="Intercept", class="sd", dpar="mu2", group="ID"), ## random effect variance here
  
  ## All mu3 values here --> intercept for 2 --> 3 values
  #prior("constant(0)", class="b", dpar="mu3"),
  prior("constant(-1.54)", coef="Intercept", class="b", dpar="mu3"),
  prior("constant(1.52)", coef="laggedState2", class="b", dpar="mu3"),
  prior("constant(3.1)", coef="laggedState3", class="b", dpar="mu3"),
  prior("constant(0)", coef="popB", class="b", dpar="mu3"),
  prior("constant(1)", coef="", class="sd", dpar="mu3") ## random effect variance here
  #prior("constant(.1)", coef="", class="sd", dpar="mu3", group="ID"), ## random effect variance here
  #prior("constant(.1)", coef="Intercept", class="sd", dpar="mu3", group="ID") ## random effect variance here
)

sampGen2 <- brm(state ~ 0 + Intercept  + laggedState + pop + (1|ID), data = pop.vals, family=categorical(), 
                sample_prior = "only", prior = priors2)
predict(sampGen2, newdata = data.frame(laggedState = c(1,2,3), pop="B", ID = "2_B"), allow_new_levels=TRUE)
## Now take these values and refit the model
post.ests <- posterior_predict(sampGen2)
## Now attach one of these back to the model
pop.vals$newEsts <- post.ests[1,]
## Fit model to these new values
fit2 <- brm(newEsts ~ 0 + Intercept  + laggedState + pop + (1|ID), data = pop.vals, family=categorical(), cores = 4)


## Now see if we can translate these models into a transition matrix
## Get the population fixed effects
out.predB <- predict(correct.mod, newdata = data.frame(laggedState = c(1,2,3), pop="B", ID = "2_C"), allow_new_levels=TRUE)
out.predBA <- pmatrix.msm(mod, covariates = list(pop="B"), t = 1)
out.predA <- predict(correct.mod, newdata = data.frame(laggedState = c(1,2,3), pop="A", ID = "2_C"), allow_new_levels=TRUE)
out.predAA <- pmatrix.msm(mod, covariates = list(pop="B"), t = 1)
## Now return all of the values of interest

## Now try to simulate data in a brms model
priors <- prior_summary(correct.mod)
sim.brm <- brm(state ~ laggedState + pop + (1|ID), data = pop.vals, family=categorical(), 
               sample_prior="only", prior = priors)
newY = predict()
bprior <- c(
  prior(normal(-0.12,0.04),class = "Intercept"),
  prior(normal(-0.14,0.01),class="b",coef="Coef1"),
  prior(normal(0.66,0.03),class="sd",coef="Intercept",group="ParticipantID"),
  prior(normal(0.07,0.01),class="sd",coef="Coef1",group="ParticipantID"),
  prior(normal(0.58,0.01),class="sigma"))


## ok I am going to try to create the coefficients from a known transiiton matrix -- without any mixed effect values here
P1 <- matrix(c( 0.9, 0.05, 0.05, 
                0.05, 0.9,   0.05, 
                .05,   .05,   .9  ), nrow=3, ncol=3, byrow = TRUE)
P2 <- P1

## Now create 50 participants from each matrix
sampSize = 2
iterLength = 5000
pop1.vals <- data.frame(ID = paste(rep(1:sampSize, each=iterLength), "A", sep="_"), time=rep(1:iterLength), 
                        state=NA, laggedState = NA,pop ="A")
pop2.vals <- data.frame(ID = paste(rep(1:sampSize, each=iterLength), "B", sep="_"), time=rep(1:iterLength), 
                        state=NA, laggedState = NA, pop ="B")
set.seed(16)
for(i in 1:sampSize){
  ## Identify sample index
  index.vals <- which(pop1.vals$ID==paste(i, "A", sep="_"))
  ## Now figure out a way to add variability to P1 somehow...
  state.vals <- run.mc.sim(P1, iterLength+1)
  lagged.vals <- lag(state.vals)[1:iterLength]
  pop1.vals[index.vals,"state"] <- state.vals[-1]
  pop1.vals[index.vals,"laggedState"] <- lagged.vals
  ## Now figure out a way to add variability to P2 somehow...
  state.vals <- run.mc.sim(P2, iterLength+1)
  lagged.vals <- lag(state.vals)[1:iterLength]
  pop2.vals[index.vals,"state"] <- state.vals[-1]
  pop2.vals[index.vals,"laggedState"] <- lagged.vals
}

## Now combine these
pop.vals <- dplyr::bind_rows(pop1.vals, pop2.vals)
msm::statetable.msm(data = pop.vals, state = state, subject = ID)
table(pop.vals$state, pop.vals$laggedState)

### Now model these
pop.vals$laggedState <- factor(pop.vals$laggedState)
pop.vals$state <- factor(pop.vals$state)
correct.modFE <- brm(state ~ 0 + Intercept + laggedState, data = pop.vals, family=categorical(), cores = 4)
## Now grab the coefficient values
predict(correct.mod, newdata = data.frame(laggedState = c(1,1,1,2,2,2,3,3,3), pop="B", ID = "2_C"), allow_new_levels=FALSE, summary=TRUE)
## Now build these out for specific matrices

## Now see how posterior predict changes these things
correct.mod <- get_prior(state ~ 0 + Intercept + laggedState, data = pop.vals, family=categorical())
priorsFE <- c(
  ## All mu2 values here --> intercept for 2 --> 3 values
  #prior("constant(0)", class="b", dpar="mu2"),
  prior("constant(-2.8)", coef="Intercept", class="b", dpar="mu2"),
  prior("constant(5.7)", coef="laggedState2", class="b", dpar="mu2"),
  prior("constant(2.8)", coef="laggedState3", class="b", dpar="mu2"),

  ## All mu3 values here --> intercept for 2 --> 3 values
  #prior("constant(0)", class="b", dpar="mu3"),
  prior("constant(-3)", coef="Intercept", class="b", dpar="mu3"),
  prior("constant(2.96)", coef="laggedState2", class="b", dpar="mu3"),
  prior("constant(5.9)", coef="laggedState3", class="b", dpar="mu3")
)
sampGen <- brm(state ~ 0 + Intercept + laggedState, data = pop.vals, family=categorical(), cores = 4,
               sample_prior = "only", prior = priorsFE, algorithm = "fixed_param")


## Now create 50 participants from each matrix
P1 <- matrix(c( 0.7, 0.15, 0.15, 
                0.15, 0.7,   0.15, 
                .15,   .15,   .7  ), nrow=3, ncol=3, byrow = TRUE)
P2 <- P1
sampSize = 2
iterLength = 5000
pop1.vals <- data.frame(ID = paste(rep(1:sampSize, each=iterLength), "A", sep="_"), time=rep(1:iterLength), 
                        state=NA, laggedState = NA,pop ="A")
pop2.vals <- data.frame(ID = paste(rep(1:sampSize, each=iterLength), "B", sep="_"), time=rep(1:iterLength), 
                        state=NA, laggedState = NA, pop ="B")
set.seed(16)
for(i in 1:sampSize){
  ## Identify sample index
  index.vals <- which(pop1.vals$ID==paste(i, "A", sep="_"))
  ## Now figure out a way to add variability to P1 somehow...
  state.vals <- run.mc.sim(P1, iterLength+1)
  lagged.vals <- lag(state.vals)[1:iterLength]
  pop1.vals[index.vals,"state"] <- state.vals[-1]
  pop1.vals[index.vals,"laggedState"] <- lagged.vals
  ## Now figure out a way to add variability to P2 somehow...
  state.vals <- run.mc.sim(P2, iterLength+1)
  lagged.vals <- lag(state.vals)[1:iterLength]
  pop2.vals[index.vals,"state"] <- state.vals[-1]
  pop2.vals[index.vals,"laggedState"] <- lagged.vals
}
pop.vals <- dplyr::bind_rows(pop1.vals, pop2.vals)
pop.vals$laggedState <- factor(pop.vals$laggedState)
pop.vals$state <- factor(pop.vals$state)
correct.modFE <- brm(state ~ 0 + Intercept + laggedState, data = pop.vals, family=categorical(), cores = 4)
priorsFE <- c(
  ## All mu2 values here --> intercept for 2 --> 3 values
  #prior("constant(0)", class="b", dpar="mu2"),
  prior("constant(-1.51)", coef="Intercept", class="b", dpar="mu2"),
  prior("constant(3.02)", coef="laggedState2", class="b", dpar="mu2"),
  prior("constant(1.47)", coef="laggedState3", class="b", dpar="mu2"),
  
  ## All mu3 values here --> intercept for 2 --> 3 values
  #prior("constant(0)", class="b", dpar="mu3"),
  prior("constant(-1.54)", coef="Intercept", class="b", dpar="mu3"),
  prior("constant(1.52)", coef="laggedState2", class="b", dpar="mu3"),
  prior("constant(3.1)", coef="laggedState3", class="b", dpar="mu3")
)
sampGen <- brm(state ~ 0 + Intercept + laggedState, data = pop.vals, family=categorical(), cores = 4,
               sample_prior = "only", prior = priorsFE, algorithm = "fixed_param")
## Now create new data sets with idenitical values


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
  n <- 100
  p <- 250
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
mod.list <- lapply(out.priors, function(x) sampGen <- brm(state ~ 0 + Intercept + laggedState + (1|ID), data = pop.vals, family=categorical(), 
                                                          cores = 4,sample_prior = "only", prior = x))
## Now create our samples
samp.gen <- lapply(mod.list, posterior_predict)


## Now fit our models to these data
## Lets start with 10 models at a time
out.mods2 <- list()
list.vals <- 1
mod.count <- 1
for(list.vals in 1:9){
for(i in 1:5){
  ## First grab the values from the data
  tmp.data <- samp.gen[[list.vals]][i,]
  ## Now attach these to our data
  pop.vals$newData <- tmp.data
  pop.vals$newLag <- NA
  ## Now create our lagged variables within each participant
  ## Now go through every unique ID and make the lagged variable within them
  for(l in unique(pop.vals$ID)){
    index <- which(pop.vals$ID == l)
    new.vals <- c(NA, lag(pop.vals$newData[index])[1:length(index)-1])
    pop.vals$newLag <- new.vals
    pop.vals$newLag <- factor(pop.vals$newLag)
    pop.vals$newData <- factor(pop.vals$newData)
  }
  ## Now estimate our model
  mod.tmp <- brm(newData ~ 0 + Intercept + newLag + (1|ID), data = pop.vals, 
                 family=categorical(),cores = 4, control = list(adapt_delta = 0.9, max_treedepth = 12))
  out.mods2[[mod.count]] <- summary(mod.tmp)
  mod.count <- mod.count + 1
}
}
