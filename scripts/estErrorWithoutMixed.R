

## Begin by loading packages and setting seed
library(msm)
library(brms)
library(foreach)
library(doParallel)
library(dplyr)

## First simulate a time sreires that is long enough to know the population fixed effects
## with a known intensitie transition matrix -- I will make three of these
## one with equal prob of trnasitioning -- one with moderate stationarity -- and one with very stable transition patterns
qmatrix.Rand <- rbind(c(-0.6,   0.4,  0.2 ),
                      c(0.4,   -0.6,  0.2 ),
                      c(0.2,   0.4, -0.6))
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
mod.estR <- brm(timeIn ~ 1  + transType, data = sim.valsR, family = weibull(), cores=4)

qmatrix.Mod <- rbind(c(-0.3,   0.2,  0.1 ),
                     c(0.15,   -0.3,  0.15 ),
                     c(0.1,   0.2, -0.3))

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
mod.estM <- brm(timeIn ~ 1 + transType, data = sim.valsM, family = weibull(), cores=4)


qmatrix.Stab <- rbind(c(-0.02,   0.015,  0.005 ),
                      c(0.01,   -0.02,  0.01 ),
                      c(0.005,   0.015, -0.02))

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
mod.estS <- brm(timeIn ~ 1 + transType, data = sim.valsS, family = weibull(), cores=4)

# ## Now calculate time spent in each observation
# prior.vals1 <- get_prior(timeIn ~ -1 + transType, data = sim.valsS, family = weibull())
# prior.vals2 <- get_prior(timeIn ~ -1 + transType + (1 | subject), data = sim.valsS, family = weibull())
# prior.vals3 <- get_prior(timeIn ~ -1 + transType + (transType | subject), data = sim.valsS, family = weibull())

## Now simulate these data 
## First do a random transition data
nAll <- c(10,50,200)
minObsAll <- c(20, 40, 80)
maxObsAll <- c(90, 140, 180)
all.samp.params <- expand.grid(nAll, minObsAll, maxObsAll)
all.params.out <- NULL
for(z in 1:10){#nrow(all.samp.params)){
  ## Now go through this process while also modifying the slope parameters
  #tmp.priors <- get_prior(timeIn ~ 1 + transType + (1 + transType|subject), data = all.data, family = weibull(), cores=4)
  ## Now assign the priors
  priorsMLString <- c(
    #prior("constant(0)", coef=""),
    #prior("constant(.2)", coef="transType11"),
    prior("constant(A)", class="Intercept"),
    #prior("constant(A)", coef="transType12"),
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
  ranef.var <- c(0,.1, .5, 1)
  ranef.var2 <- c(0, .01,.05, .1)
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
  #mod.list <- lapply(out.priors, function(x) sampGen <- brm(timeIn ~ 1 + transType + (1 + transType|subject), data = all.data, family = weibull(), cores=4, prior = x,sample_prior = "only", seed = seedVal))
  ## Now create our samples
  #samp.gen <- lapply(mod.list, posterior_predict)
  
  ## Now estimate these models ignoring the random effect coefficients
  cl <- makeCluster(8) #not to overload your computer
  registerDoParallel(cl)
  out.mods2 <- list()
  list.vals <- 1
  mod.count <- 1
  for(list.vals in 1:length(out.priors)){
    tmp.out <- foreach(t = 1:50, .packages = c("brms", "dplyr", "msm"), .errorhandling = "remove") %dopar%{
      ## Idenitfy our unique file identifier
      out.file <- paste("./data/individualSims/",z,"_", list.vals, "_", t, ".csv", sep="")
      if(file.exists(out.file)){
        out2 <- read.csv(out.file)
      }else{
        ## First create our new dataset
        seedVal = sample(1:100000, size = 1)
        n <- all.samp.params[z,1]
        minObs <- all.samp.params[z,2]
        maxObs <- all.samp.params[z,3]
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
        ## Now sample priors from these data -- only need about 100?
        sampGen <- brm(timeIn ~ 1 + transType + (1 + transType|subject), data = all.data, family = weibull(), 
                       cores=1, prior = out.priors[[list.vals]],sample_prior = "only", seed = seedVal, iter = 100, chains = 1)
        ## Now grab our estimated values
        tmp.data <- posterior_predict(sampGen)
        ## Now attach these to our data
        all.data$newData <- tmp.data[sample(1:50, size = 1),]
        ## Now estimate our model
        mod.tmp <- brm(newData ~ 1 + transType, data = all.data, 
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
        out2$paramTrue <- trueVals[1:nrow(out2),"Val"]
        randRow <- data.frame(paramTrue = trueVals$Val[which(trueVals$group=="subject" & trueVals$coef!="")],
                              coef = trueVals$coef[which(trueVals$group=="subject" & trueVals$coef!="")])
        out2 <- merge(out2, randRow, by="coef", suffixes = c("", "_Rand"), all=TRUE)
        #out2 <- bind_rows(out2, randRow)
        out2$Row <- t
        out2$n <- n
        out2$minP <- minObs
        out2$maxP <- maxObs
        ## Now combine the true values and the estimated values 
        rm(mod.tmp, out, outFixed, trueVals)
        ## Now organize the data
        write.csv(out2, file = out.file, quote=F, row.names = F)
      }
      out2
    }
    out.mods2[[mod.count]] <- tmp.out
    mod.count <- mod.count + 1
  }
  stopCluster(cl)
  out.mods2<- bind_rows(out.mods2)
  all.params.out <- rbind(all.params.out, out.mods2)
  print(z)
}
saveRDS(all.params.out, file = "./data/estErrorWithoutMixed.RDS")


## Now perform our ANOA on these results here
all.params.out <- readRDS("./data/estErrorWithoutMixed.RDS")


## Now organize the data
all.params.out$estError <- as.numeric(all.params.out$paramEst) - as.numeric(all.params.out$paramTrue)
## Now make factors where factors need to be
all.params.out$coef <- factor(all.params.out$coef)
all.params.out$minP <- factor(all.params.out$minP)
all.params.out$maxP <- factor(all.params.out$maxP)
all.params.out$n <- factor(all.params.out$n)
all.params.out$paramTrue_Rand <- factor(all.params.out$paramTrue_Rand)
## Now add an indicator variable if these are coming from a stable MM, moderately stable or random MM

all.params.mod <- all.params.out[-which(all.params.out$coef=="Intercept"),]
## Drop intercept term
mod <- lm(estError ~ (coef + minP + paramTrue_Rand)^3, data = all.params.out)
## Now examine only intercept term
