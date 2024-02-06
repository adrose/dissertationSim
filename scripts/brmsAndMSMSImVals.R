## Oringal code
library(msm)
library(brms)
library(foreach)
library(doParallel)
library(dplyr)

## Set seed... if that does anything?
seedVal <- 16
set.seed(16)

sim.df <- data.frame(subject = rep(1:100, rep(13,100)), time = rep(seq(0, 24, 2), 100))
qmatrix <- rbind(c(-0.11,   0.1,  0.01 ),
                 c(0.05,   -0.15,  0.1 ),
                 c(0.02,   0.07, -0.09))
simmulti.msm(sim.df, qmatrix)


## Adon explore
sim.df <- data.frame(subject = 1, time = 1:100 + runif(100, min = 0, max = .9))
qmatrix <- rbind(c(-0.2,   0.1,  0.1 ),
                 c(0.1,   -0.2,  0.1 ),
                 c(0.1,   0.1, -0.2))
sim.vals <- simmulti.msm(sim.df, qmatrix)
## now estimate the model
tmp.mod <- msm(state ~ time, data = sim.vals, gen.inits = TRUE, qmatrix = qmatrix)


## Now combine the msm package and the data simulation via brms code into one script so 
## we know the population trnasiiton parameters??

## First simulate a time sreires that is long enough to know the population fixed effects
## with a known intensitie transition matrix -- I will make three of these
## one with equal prob of trnasitioning -- one with moderate stationarity -- and one with very stable transition patterns
qmatrix.Rand <- rbind(c(-0.6,   0.3,  0.3 ),
                      c(0.3,   -0.6,  0.3 ),
                      c(0.3,   0.3, -0.6))
## Now simulate the data for a single subject
sim.dfR <- data.frame(subject = 1, timeR = runif(10001, min = 0, max = 10))
sim.dfR$time <- cumsum(sim.dfR$timeR)
sim.vals <- simmulti.msm(sim.dfR, qmatrix.Rand)
## Now create our transition type varaible
sim.vals$lagState <- dplyr::lag(sim.vals$state)
sim.vals$transType <- paste(sim.vals$lagState, sim.vals$state, sep="")
sim.vals$timeIn <- dplyr::lag(sim.vals$time)
## Now remove intrastate transitions
sim.vals <- sim.vals[-which(sim.vals$transType %in% c("11", "22", "33")),]
sim.vals$timeIn <- c(NA, diff(sim.vals$time))
sim.valsR <- na.omit(sim.vals)

## Now calculate time spent in each observation
mod.estR <- brm(timeIn ~ -1 + transType, data = sim.valsR, family = weibull(), cores=4)

qmatrix.Mod <- rbind(c(-0.3,   0.15,  0.15 ),
                      c(0.15,   -0.3,  0.15 ),
                      c(0.15,   0.15, -0.3))

## Now simulate the data for a single subject
sim.dfR <- data.frame(subject = 1, timeR = runif(10001, min = 0, max = 10))
sim.dfR$time <- cumsum(sim.dfR$timeR)
sim.vals <- simmulti.msm(sim.dfR, qmatrix.Mod)
## Now create our transition type varaible
sim.vals$lagState <- dplyr::lag(sim.vals$state)
sim.vals$transType <- paste(sim.vals$lagState, sim.vals$state, sep="")
sim.vals$timeIn <- dplyr::lag(sim.vals$time)
sim.vals <- sim.vals[-which(sim.vals$transType %in% c("11", "22", "33")),]
sim.vals$timeIn <- c(NA, diff(sim.vals$time))
sim.valsM <- na.omit(sim.vals)

## Now calculate time spent in each observation
mod.estM <- brm(timeIn ~ -1 + transType, data = sim.valsM, family = weibull(), cores=4)


qmatrix.Stab <- rbind(c(-0.02,   0.01,  0.01 ),
                      c(0.01,   -0.02,  0.01 ),
                      c(0.01,   0.01, -0.02))

## Now simulate the data for a single subject
sim.dfR <- data.frame(subject = 1, timeR = runif(10001, min = 0, max = 10))
sim.dfR$time <- cumsum(sim.dfR$timeR)
sim.vals <- simmulti.msm(sim.dfR, qmatrix.Stab)
## Now create our transition type varaible
sim.vals$lagState <- dplyr::lag(sim.vals$state)
sim.vals$transType <- paste(sim.vals$lagState, sim.vals$state, sep="")
sim.vals <- sim.vals[-which(sim.vals$transType %in% c("11", "22", "33")),]
sim.vals$timeIn <- c(NA, diff(sim.vals$time))
sim.valsS <- na.omit(sim.vals)

## Now estimate our model
mod.estS <- brm(timeIn ~ -1 + transType, data = sim.valsS, family = weibull(), cores=4)

# ## Now calculate time spent in each observation
# prior.vals1 <- get_prior(timeIn ~ -1 + transType, data = sim.valsS, family = weibull())
# prior.vals2 <- get_prior(timeIn ~ -1 + transType + (1 | subject), data = sim.valsS, family = weibull())
# prior.vals3 <- get_prior(timeIn ~ -1 + transType + (transType | subject), data = sim.valsS, family = weibull())

## Now simulate these data 
## First do a random transition data
n <- 30
minObs <- 40
maxObs <- 100
p <- sample(minObs:maxObs, size = n, replace = TRUE)
all.data <- NULL
for(i in 1:n){
  sim.tmp <- data.frame(subject=i, time = cumsum(runif(p[i]+1, min = .5, max = 10)))
  sim.vals <- simmulti.msm(sim.tmp, qmatrix.Rand)
  sim.vals$lagState <- dplyr::lag(sim.vals$state)
  sim.vals$transType <- paste(sim.vals$lagState, sim.vals$state, sep="")
  sim.vals <- sim.vals[-which(sim.vals$transType %in% c("11", "22", "33")),]
  sim.vals$timeIn <- c(NA, diff(sim.vals$time))
  sim.vals <- na.omit(sim.vals)
  all.data <- dplyr::bind_rows(all.data, sim.vals)
}

## Now model these data
tmp.mod <- brm(timeIn ~ -1 + transType + (1|subject), data = all.data, family = weibull(), cores=4)
tmp.priors <- get_prior(timeIn ~ -1 + transType + (1|subject), data = all.data, family = weibull(), cores=4)
## Now assign the priors
priorsMLString <- c(
  #prior("constant(0)", coef=""),
  #prior("constant(.2)", coef="transType11"),
  prior("constant(A)", coef="transType12"),
  prior("constant(B)", coef="transType13"),
  prior("constant(C)", coef="transType21"),
  #prior("constant(.2)", coef="transType22"),
  prior("constant(D)", coef="transType23"),
  prior("constant(E)", coef="transType31"),
  prior("constant(F)", coef="transType32"),
  #prior("constant(.2)", coef="transType33"),
  #prior("constant(1)", coef="", class="sd", group=="", coef=""), 
  #prior("constant(1)", coef="", class="sd", group=="subject", coef=""),
  prior("constant(G)", class="sd", group="subject", coef="Intercept"), 
  prior("constant(1.51)", class="shape", lb="0")## Weibull shape parameter here
)

#samp.gen <- brm(timeIn ~ -1 + transType + (1|subject), data = all.data, family = weibull(), cores=4, prior = priorsMLString,sample_prior = "only")
#sampVals <- posterior_predict(samp.gen)
## Now estimate these data
#all.data$gen1 <- sampVals[1,]
#samp.test <- brm(gen1 ~ -1 + transType + (1|subject), data = all.data, family = weibull(), cores=4,
#                 control = list(adapt_delta = 0.9, max_treedepth = 12), thin = 5, iter = 12000, warmup = 4000)

## Now loop through each of these models, and change the participant intercept varaince, and potentially scale up to the
## slopes for each transition type as well
out.models <- list(mod.estR, mod.estM, mod.estS)
ranef.var <- c(.1, 1, 5)
out.priors <- list()
out.priors.index <- 1
for(r in 1:length(ranef.var)){
  for(i in 1:length(out.models)){
    ## First change the fixef variables
    newprior <- priorsMLString
    fixefVals <- fixef(out.models[[i]])
    for(j in 1:7){
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
mod.list <- lapply(out.priors, function(x) sampGen <- brm(timeIn ~ -1 + transType + (1|subject), data = all.data, family = weibull(), cores=4, prior = x,sample_prior = "only"))
## Now create our samples
samp.gen <- lapply(mod.list, posterior_predict)

## Now estimate the models
out.mods2 <- list()
list.vals <- 1
mod.count <- 1
## Create a do parallel loop here
cl <- makeCluster(3) #not to overload your computer
registerDoParallel(cl)
for(list.vals in 1:list.vals){
  tmp.out <- foreach(i = 1:12, .packages = c("brms", "dplyr")) %dopar%{
    ## First grab the values from the data
    tmp.data <- samp.gen[[list.vals]][i,]
    ## Now attach these to our data
    all.data$newData <- tmp.data
    ## Now estimate our model
    mod.tmp <- brm(newData ~ -1 + transType + (1|subject), data = all.data, 
                   family=weibull(),cores = 1, control = list(adapt_delta = 0.9, max_treedepth = 12), chains = 3)
    out <- summary(mod.tmp)
    outFixed <- out$fixed
    outRand <- bind_rows(out$random)
    trueVals <- out.priors[[list.vals]]
    trueVals$Val <- gsub(x = trueVals$prior, pattern = "[^0-9.-]", replacement = "")
    ## Now prep the output
    out2 <- data.frame(coef=c(rownames(outFixed)), paramEst = NA, lowerEst = NA, upperEst = NA, paramTrue = NA)
    out2$paramEst <- outFixed[out2$coef,"Estimate"]
    out2$lowerEst <- outFixed[out2$coef,"l-95% CI"]
    out2$upperEst <- outFixed[out2$coef,"u-95% CI"]
    out2$paramTrue <- trueVals[which(trueVals$coef %in% out2$coef),"Val"]
    randRow <- data.frame(coef=rownames(outRand), paramEst = outRand$Estimate, lowerEst = outRand$`l-95% CI`, 
                          upperEst = outRand$`u-95% CI`, paramTrue = trueVals$Val[which(trueVals$group=="subject")])
    out2 <- bind_rows(out2, randRow)
    out2$Row <- i
    ## Now combine the true values and the estimated values 
    rm(mod.tmp, out, outFixed, outRand, trueVals)
    ## Now organize the data
    out2
  }
  out.mods2[[mod.count]] <- bind_rows(tmp.out)
  mod.count <- mod.count + 1
}
stopCluster(cl)

out.name <- paste("./data/sample_", length(unique(all.data$subject)), "_sampMin_", length(unique(all.data$time)), ".RDS", sep='')
saveRDS(out.mods2, file = out.name)



## Now go through this process while also modifying the slope parameters
tmp.priors <- get_prior(timeIn ~ -1 + transType + (transType|subject), data = all.data, family = weibull(), cores=4)
## Now assign the priors
priorsMLString <- c(
  #prior("constant(0)", coef=""),
  #prior("constant(.2)", coef="transType11"),
  prior("constant(A)", coef="transType12"),
  prior("constant(B)", coef="transType13"),
  prior("constant(C)", coef="transType21"),
  #prior("constant(.2)", coef="transType22"),
  prior("constant(D)", coef="transType23"),
  prior("constant(E)", coef="transType31"),
  prior("constant(F)", coef="transType32"),
  #prior("constant(.2)", coef="transType33"),
  #prior("constant(1)", coef="", class="sd", group=="", coef=""), 
  #prior("constant(1)", coef="", class="sd", group=="subject", coef=""),
  prior("lkj(1)", class="cor", group="subject"),
  prior("constant(G)", class="sd", group="subject", coef="Intercept"),
  prior("constant(H)", class="sd", group="subject", coef="transType13"),
  prior("constant(H)", class="sd", group="subject", coef="transType21"),
  prior("constant(H)", class="sd", group="subject", coef="transType23"),
  prior("constant(H)", class="sd", group="subject", coef="transType31"),
  prior("constant(H)", class="sd", group="subject", coef="transType32"),
  prior("constant(1.51)", class="shape", lb="0")## Weibull shape parameter here
)

#samp.gen <- brm(timeIn ~ -1 + transType + (1|subject), data = all.data, family = weibull(), cores=4, prior = priorsMLString,sample_prior = "only")
#sampVals <- posterior_predict(samp.gen)
## Now estimate these data
#all.data$gen1 <- sampVals[1,]
#samp.test <- brm(gen1 ~ -1 + transType + (1|subject), data = all.data, family = weibull(), cores=4,
#                 control = list(adapt_delta = 0.9, max_treedepth = 12), thin = 5, iter = 12000, warmup = 4000)

## Now loop through each of these models, and change the participant intercept varaince, and potentially scale up to the
## slopes for each transition type as well
out.models <- list(mod.estR, mod.estM, mod.estS)
ranef.var <- c(0, .1, .5, 1)
ranef.var2 <- c(0, .01, .05, .1)
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
      if(j == 7){
        newprior$prior <- gsub(pattern = LETTERS[j], replacement = ranef.var[r], x = newprior$prior)
      }
      if(j == 8 ){
        newprior$prior <- gsub(pattern = LETTERS[j], replacement = ranef.var2[r], x = newprior$prior)
      }
    }
    out.priors[[out.priors.index]] <- newprior
    out.priors.index <- out.priors.index + 1
  }
}
## Now estimate our models
mod.list <- lapply(out.priors, function(x) sampGen <- brm(timeIn ~ -1 + transType + (transType|subject), data = all.data, family = weibull(), cores=4, prior = x,sample_prior = "only", seed = seedVal))
## Now create our samples
samp.gen <- lapply(mod.list, posterior_predict)

## Now estimate these models ignoring the random effect coefficients
cl <- makeCluster(5) #not to overload your computer
registerDoParallel(cl)
out.mods2 <- list()
list.vals <- 1
mod.count <- 1
for(list.vals in 1:length(samp.gen)){
  tmp.out <- foreach(i = 1:100, .packages = c("brms", "dplyr")) %dopar%{
    ## First grab the values from the data
    tmp.data <- samp.gen[[list.vals]][i,]
    ## Now attach these to our data
    all.data$newData <- tmp.data
    ## Now estimate our model
    mod.tmp <- brm(newData ~ -1 + transType, data = all.data, 
                   family=weibull(),cores = 1, control = list(adapt_delta = 0.9, max_treedepth = 12), chains = 3)
    out <- summary(mod.tmp)
    outFixed <- out$fixed
    #outRand <- bind_rows(out$random)
    trueVals <- out.priors[[list.vals]]
    trueVals$Val <- gsub(x = trueVals$prior, pattern = "[^0-9.-]", replacement = "")
    ## Now prep the output
    out2 <- data.frame(coef=c(rownames(outFixed)), paramEst = NA, lowerEst = NA, upperEst = NA, paramTrue = NA)
    out2$paramEst <- outFixed[out2$coef,"Estimate"]
    out2$lowerEst <- outFixed[out2$coef,"l-95% CI"]
    out2$upperEst <- outFixed[out2$coef,"u-95% CI"]
    out2$paramTrue <- trueVals[which(trueVals$coef[which(trueVals$group=="")] %in% out2$coef),"Val"]
    randRow <- data.frame(paramTrue = trueVals$Val[which(trueVals$group=="subject" & trueVals$coef!="")],
                          coef = trueVals$coef[which(trueVals$group=="subject" & trueVals$coef!="")])
    out2 <- merge(out2, randRow, by="coef", suffixes = c("", "_Rand"), all=TRUE)
    #out2 <- bind_rows(out2, randRow)
    out2$Row <- i
    ## Now combine the true values and the estimated values 
    rm(mod.tmp, out, outFixed, trueVals)
    ## Now organize the data
    out2
  }
  out.mods2[[mod.count]] <- bind_rows(tmp.out)
  mod.count <- mod.count + 1
}
stopCluster(cl)

