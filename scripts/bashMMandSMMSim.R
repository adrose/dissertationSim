### Now make an R script which can be called by the shell and will run one row 
### from the simulateWeibul.R script -- I think this should help to avoid the following error:
### Error in dyn.load(libLFile) : 
### unable to load shared object '/tmp/RtmphUpD1w/filebef1be2acc8.so':
###  `maximal number of DLLs reached...

## First identify the row of the simulation we want to run
seedVal=as.numeric(commandArgs(TRUE)[1])
n <- c(20, 100)
minObsAll <- c(20, 40)
n.states <- c(3)
matrixType <- c("mod", "rand")
scaleVals <- c(".5:8","8:16")
shapeVals <- c(1, 5)
mainEffectVals <- c(0,.6)
rand.var <- c(0.1,3)
iter.vals <- 1:1000
all.parms <- expand.grid(n, minObsAll, n.states, matrixType, scaleVals, shapeVals, mainEffectVals,rand.var,iter.vals)
set.seed(seedVal)


for(i in sample(which(all.parms[,9]==seedVal), replace = FALSE)){
  ## Now declare output file
  out.file <- paste("./data/individualSimsMM_SMM-3/seedVal_", seedVal,"/rowVal_",i, "_seedVal_", seedVal, ".RDS", sep='')
  out.dir <- paste("./data/individualSimsMM_SMM-3/seedVal_", seedVal,sep='')
  if(!dir.exists(out.dir)){
    dir.create(out.dir)
    print("make dir")
  }
  if(!file.exists(out.file)){
    
    ## Load library(s)
    library(brms)
    source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
    source("./scripts/weiFuncs.R")
    
    ## Now declare all of the simulation components
    #mainEffectMag <- c(0, .2, .4)
    # Should be: rand.var.int <- c(.05, .01)
    # Should be: rand.var.slope <- c(.05, .1)
    ## Now go ahead and loop through all of these params real quick
    ## First create the data
    all.mods <- list()
    mod.count <- 1
    tmp.dat <- simFunc(n = all.parms[i,1], minObs = all.parms[i,2],maxObs = all.parms[i,2]*2,nState = all.parms[i,3], 
                       matrixType = all.parms[i,4], scaleRange = all.parms[i,5], shapeVal = all.parms[i,6], meMag = all.parms[i,7], addLag = FALSE)
    if(sum(is.na(tmp.dat$sampleData$timeIn))>0){
      sum.check <- TRUE
      while(sum(is.na(tmp.dat$sampleData$timeIn))>0){
        tmp.dat <- simFunc(n = all.parms[i,1], minObs = all.parms[i,2],maxObs = all.parms[i,2]*2,nState = all.parms[i,3], 
                           matrixType = all.parms[i,4], scaleRange = all.parms[i,5], shapeVal = all.parms[i,6], meMag = all.parms[i,7], addLag = FALSE)
        if(sum(is.na(tmp.dat$sampleData$timeIn))==0){
          sum.check <- FALSE
        }
        print("sumCheck")
      }
      
    }
    ## Modify the transTYpe values to reflect the model names
    tmp.dat$transVals$coefName <- paste("transType", gsub(x = tmp.dat$transVals$transType, pattern=" ", replacement=""), sep='')
    ## First get the priors of interest
    priorVar <- get_prior(timeIn ~ -1 + transType + effectOfInt + (1|part), data = tmp.dat$sampleData, family=weibull())
    ## Now fix the priors here
    ## Use the weibull parameters from the tmp.dat$transVals
    for(w in 1:nrow(tmp.dat$transVals)){
      ## First identify the true value from the simulated data
      modIndex <- which(tmp.dat$transVals$coefName[w] == priorVar$coef & priorVar$class=="b")
      ## Get the prior value we need
      priorValue <- log(tmp.dat$transVals$scale[w])
      newString <- tmp.dat$transVals$coefName[w]
      randVarVal <- all.parms[i,8]
      if(randVarVal == 0){
        randVarVal <- .1
      }
      if(randVarVal == 1){
        randVarVal <- 3
      }
      newPrior <- prior(normal(QRT,.3), class = "b", coef = XYZ)
      ## Now change the values to what we need
      newPrior$prior <- gsub(x = newPrior$prior, replacement = priorValue, pattern = "QRT")
      newPrior$prior <- gsub(x = newPrior$prior, replacement = randVarVal, pattern = "ABC")
      newPrior$coef <- gsub(x = newPrior$coef, replacement = newString, pattern = "XYZ")
      
      ## Now put this back into the priorVar
      priorVar[modIndex,] <- newPrior
    }
    ## Now do the effect of interest
    ## Find the row index
    mePrior <- prior(normal(QRT, .3), class = "shape", lb=.02, ub=9)
    ## Now change the values to what we need
    mePrior$prior <- gsub(x = mePrior$prior, replacement = all.parms[i,6], pattern = "QRT")
    priorVar[which(priorVar$class=="shape"),] <- mePrior
    ## Now do the main effect of interest
    mePrior <- prior(normal(QRT,.4), class = "b", coef = XYZ)
    ## Now change the values to what we need
    mePrior$prior <- gsub(x = mePrior$prior, replacement = all.parms[i,7], pattern = "QRT")
    mePrior$coef <- gsub(x = mePrior$coef, pattern = "XYZ", replacement = "effectOfInt")
    priorVar[which(priorVar$coef=="effectOfInt"),] <- mePrior
    ## Now do the random intercept
    randEfInt <- prior(constant(QRT), class="sd", coef="Intercept", group="part")
    randEfInt$prior <- gsub(x = randEfInt$prior, replacement = all.parms[i,8], pattern = "QRT")
    priorVar[which(priorVar$group=="part" & priorVar$coef=="Intercept"),] <- randEfInt 
    
    ## Now make the four priors we need -- two with frailty terms -- two without
    sampGen <- brm(timeIn ~ -1 + transType + effectOfInt + (1|part), data = tmp.dat$sampleData, family = weibull(), 
                   cores=1, prior = priorVar,sample_prior = "only", seed = 16, iter = 500, chains = 1)
    tmp.data <- predict(sampGen,newdata = tmp.dat$sampleData ,summary=FALSE, ndraws=250)
    tmp.dat$sampleData$genVals <- NA
    ## Now grab some values from each of these
    for(sampleIndiv in unique(tmp.dat$sampleData$part)){
      ## grab the index values
      part.count <- which(tmp.dat$sampleData$part==sampleIndiv)
      tmp.dat$sampleData$genVals[part.count] <- tmp.data[sample(1:250, size = 1),part.count]
    }
    ## Now modify the trans type to ignore the lagged influence
    tmp.dat$sampleData$transType1 <- paste(tmp.dat$sampleData$stateFrom, tmp.dat$sampleData$stateTo)
    ## SMM no frailty term
    mod.est1 <- brm(genVals ~ -1  + transType1 + effectOfInt, data = tmp.dat$sampleData, family = weibull(), cores=1,
                    control = list(adapt_delta = 0.9, max_treedepth = 12), iter = 5000, init = 2000, thin = 3, chains = 2)
    all.mods[[mod.count]] <- mod.est1
    mod.count <- mod.count + 1
    ## SMM with frailty term
    mod.est2 <- brm(genVals ~ -1  + transType1 + effectOfInt + (1|part), data = tmp.dat$sampleData, family = weibull(), cores=1,
                    control = list(adapt_delta = 0.9, max_treedepth = 12), iter = 5000, init = 2000, thin = 3, chains = 2)
    all.mods[[mod.count]] <- mod.est2
    mod.count <- mod.count + 1
    ## Now do the markov model terms here
    mod.est3 <- brm(genVals ~ -1  + transType1 + effectOfInt, data = tmp.dat$sampleData, family = brmsfamily("exponential"),
                    cores=1,control = list(adapt_delta = 0.9, max_treedepth = 12), iter = 5000, init = 2000, thin = 3, chains = 2)
    all.mods[[mod.count]] <- mod.est3
    mod.count <- mod.count + 1
    ## Now do MM with frail
    mod.est4 <- brm(genVals ~ -1  + transType1 + effectOfInt + (1|part) , data = tmp.dat$sampleData, family = brmsfamily("exponential"), 
                    cores=1,control = list(adapt_delta = 0.9, max_treedepth = 12), iter = 5000, init = 2000, thin = 3, chains = 2)
    all.mods[[mod.count]] <- mod.est4
    mod.count <- mod.count + 1
    
    ## Now do the priors
    all.mods[[mod.count]]  <- priorVar
    ## Now write all.out
    saveRDS(all.out, file = out.file)
    
  }else{
    print("Job done")
  }
}

