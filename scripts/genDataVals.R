## This script will be used to create a series of population observations for the weibull model estimates
## These data will come from models where I will vary the:
## 1. Transition intensities
## 2. random effect variance (only intercept?)
## 3. random coefficient variance?
## 4. sample size; number of observations
## 5. the percent of individuals that come from each model varying the above ^^



## Begin by loading packages and setting seed
library(msm)
library(brms)
library(foreach)
library(doParallel)
library(dplyr)

## Set seed... if that does anything?
seedVal <- 16
set.seed(16)


## This next portion of code estimates the fixed effect parameters that we need to use to generate
## the data --- this can easily be replaced by just allowing the parameters to be fixed, but these shouldn't change very much

## First simulate a time series that is long enough to know the population fixed effects
## with a known intensity transition matrix -- I will make three of these
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

########
########
########
######## This is where the data generation process actually begins
######## The n, min, and max obs are going to need to be varried to alter the sample
######## the other characteristics (the fixed parameters) will not vary across samples... duh
######## Lets see if we can loop through these values and produce each one of these across a range of population levels
nAll <- c(10,50,100)
minObsAll <- c(10, 50, 100)
maxObsAll <- c(50, 100, 200)

for(z in 1:length(nAll)){
  seedVal = sample(1:100000, size = 1)
  n <- nAll[z]
  minObs <- minObsAll[z]
  maxObs <- maxObsAll[z]
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
  
  ## Now collapse all of our values
  data.out <- lapply(samp.gen, t)
  ## Now go through each of these and fix the column names
  for(i in 1:length(data.out)){
    colnames(data.out[[i]]) <- paste("mod_", i, "_obs_", 1:nrow(samp.gen[[i]]), sep='')
  }
  ## Now attach these to the all data
  data.out <- bind_cols(data.out)
  all.data.out <- bind_cols(all.data, data.out)
  
  ## now save this with the correct indicator variables
  out.name <- paste("./data/sampleSize_", n, "_minObs_", minObs, "_maxObs_", maxObs,"_seed_",seedVal, ".csv", sep='')
  write.csv(x = all.data.out, file = out.name, quote = FALSE, row.names = FALSE)
}
