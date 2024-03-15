### Now make an R script which can be called by the shell and will run one row 
### from the simulateWeibul.R script -- I think this should help to avoid the following error:
### Error in dyn.load(libLFile) : 
### unable to load shared object '/tmp/RtmphUpD1w/filebef1be2acc8.so':
###  `maximal number of DLLs reached...

## First identify the row of the simulation we want to run
i=as.numeric(commandArgs(TRUE)[1])
subVal <- as.numeric(commandArgs(TRUE)[2])
n <- c(20, 100)
minObsAll <- c(20, 50)
maxObsAll <- c(80, 200)
n.states <- c(3,4)
matrixType <- c("mod", "rand")
scaleVals <- c("2:8")
shapeVals <- c(1, 2)
mainEffectVals <- c(0, .2, .4)
rand.var <- c(0, 1)
iter.vals <- 1:300
if(!is.null(subVal)){
  i <- i + subVal*1000
}
all.parms <- expand.grid(n, minObsAll, n.states, matrixType, scaleVals, shapeVals, mainEffectVals,rand.var,iter.vals)
seedVal <- all.parms[i,9]
set.seed(all.parms[i,9])

## Now declare output file
out.file <- paste("./data/individualSimsMM_SMM/rowVal_", i, "_seedVal_", seedVal, ".csv", sep='')
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
                     matrixType = all.parms[i,4], scaleRange = all.parms[i,5], shapeVal = all.parms[i,6], meMag = all.parms[i,7])
  ## Modify the transTYpe values to reflect the model names
  tmp.dat$transVals$coefName <- paste("transType", gsub(x = tmp.dat$transVals$transType, pattern=" ", replacement=""), sep='')
  ## First get the priors of interest
  priorVar <- get_prior(timeIn ~ -1 + transType + effectOfInt + (-1 + transType|part), data = tmp.dat$sampleData, family=weibull())
  ## Now fix the priors here
  ## Use the weibull parameters from the tmp.dat$transVals
  for(w in 1:nrow(tmp.dat$transVals)){
    ## First idenitfy the true value from the simulated data
    modIndex <- which(tmp.dat$transVals$coefName[w] == priorVar$coef & priorVar$class=="b")
    ## Get the prior value we need
    priorValue <- log(tmp.dat$transVals$scale[w])
    newString <- tmp.dat$transVals$coefName[w]
    newPrior <- prior(normal(QRT,.1), class = "b", coef = XYZ)
    ## Now change the values to what we need
    newPrior$prior <- gsub(x = newPrior$prior, replacement = priorValue, pattern = "QRT")
    newPrior$coef <- gsub(x = newPrior$coef, replacement = newString, pattern = "XYZ")
    ## Now put this back into the priorVar
    priorVar[modIndex,] <- newPrior
  }
  ## Now do the same for the random effect
  for(w in 1:nrow(tmp.dat$transVals)){
    ## First idenitfy the true value from the simulated data
    modIndex <- which(tmp.dat$transVals$coefName[w] == priorVar$coef & priorVar$class=="sd")
    ## Get the prior value we need
    priorValue <- all.parms[i,8]
    newString <- tmp.dat$transVals$coefName[w]
    newPrior <- prior(constant(QRT), class = "sd", coef = XYZ, group=part)
    ## Now change the values to what we need
    newPrior$prior <- gsub(x = newPrior$prior, replacement = priorValue, pattern = "QRT")
    newPrior$coef <- gsub(x = newPrior$coef, replacement = newString, pattern = "XYZ")
    ## Now put this back into the priorVar
    priorVar[modIndex,] <- newPrior
  }
  ## Now do the effect of interest
  ## Find the row index
  mePrior <- prior(constant(QRT), class = "shape", lb=0)
  ## Now change the values to what we need
  mePrior$prior <- gsub(x = mePrior$prior, replacement = all.parms[i,6], pattern = "QRT")
  priorVar[which(priorVar$class=="shape"),] <- mePrior
  ## Now do the main effect of interest
  mePrior <- prior(normal(QRT, .1), class = "b", coef = XYZ)
  ## Now change the values to what we need
  mePrior$prior <- gsub(x = mePrior$prior, replacement = all.parms[i,7], pattern = "QRT")
  mePrior$coef <- gsub(x = mePrior$coef, pattern = "XYZ", replacement = "effectOfInt")
  priorVar[which(priorVar$coef=="effectOfInt"),] <- mePrior
  
  ## Now make the four priors we need -- two with frailty terms -- two without
  sampGen <- brm(timeIn ~ -1 + transType + effectOfInt + (-1 + transType|part), data = tmp.dat$sampleData, family = weibull(), 
                 cores=1, prior = priorVar,sample_prior = "only", seed = 16, iter = 10, chains = 1)
  tmp.data <- predict(sampGen,newdata = tmp.dat$sampleData ,summary=FALSE, ndraws=3)
  tmp.dat$sampleData$genVals <- tmp.data[1,]
  ## SMM no frailty term
  mod.est1 <- brm(genVals ~ -1  + transType + effectOfInt, data = tmp.dat$sampleData, family = weibull(), cores=1,
                  control = list(adapt_delta = 0.9, max_treedepth = 12), iter = 5000, init = 2000, thin = 3, chains = 2)
  all.mods[[mod.count]] <- mod.est1
  mod.count <- mod.count + 1
  ## SMM with frailty term
  mod.est2 <- brm(genVals ~ -1  + transType + effectOfInt + (1|part), data = tmp.dat$sampleData, family = weibull(), cores=1,
                  control = list(adapt_delta = 0.9, max_treedepth = 12), iter = 5000, init = 2000, thin = 3, chains = 2)
  all.mods[[mod.count]] <- mod.est2
  mod.count <- mod.count + 1
  ## Now do the markov model terms here
  ## Constrain the prior to be 1 for the shape term
  priorMMShape <- get_prior(genVals ~ -1  + transType + effectOfInt, data = tmp.dat$sampleData, family=weibull())
  priorMMShape[which(priorMMShape$class=="shape"),"prior"] <- "constant(1)"
  ## Now estimate the model -- again MM without frailty
  mod.est3 <- brm(genVals ~ -1  + transType + effectOfInt, data = tmp.dat$sampleData, family = weibull(), 
                  cores=1,control = list(adapt_delta = 0.9, max_treedepth = 12), prior = priorMMShape, iter = 5000, init = 2000, thin = 3, chains = 2)
  all.mods[[mod.count]] <- mod.est3
  mod.count <- mod.count + 1
  ## Now do MM with frail
  priorMMShape <- get_prior(genVals ~ -1  + transType + effectOfInt + (1|part), data = tmp.dat$sampleData, family=weibull())
  priorMMShape[which(priorMMShape$class=="shape"),"prior"] <- "constant(1)"
  mod.est4 <- brm(genVals ~ -1  + transType + effectOfInt + (1|part) , data = tmp.dat$sampleData, family = weibull(), 
                  cores=1,control = list(adapt_delta = 0.9, max_treedepth = 12), prior = priorMMShape, iter = 5000, init = 2000, thin = 3, chains = 2)
  all.mods[[mod.count]] <- mod.est4
  mod.count <- mod.count + 1

  ## Now prep all of the output
  all.out <- NULL
  for(m in 1:length(all.mods)){
    outFixed <- data.frame(summary(all.mods[[m]])$fixed)
    outFixed$expEstimate <- NA
    outFixed$expLower <- NA
    outFixed$expUpper <- NA
    ## Put these back in the original units
    for(r in 1:nrow(outFixed)){
      outFixed$expEstimate[r] <- exp(sum(outFixed$Estimate[c(r)]))
      outFixed$expLower[r] <-   exp(sum(outFixed$l.95..CI[c(r)]))
      outFixed$expUpper[r] <-   exp(sum(outFixed$u.95..CI[c(r)]))
    }
    ## Now prepare all output
    out2 <- data.frame(transType=c(rownames(outFixed)), paramEst = NA, lowerEst = NA, upperEst = NA,
                       randVar = all.parms[i,8], rowAllParam = i, n = all.parms[i,1], minObs = all.parms[i,2],
                       nState = all.parms[i,3], matrixType = all.parms[i,4], scaleRange = all.parms[i,5], shapeVal = all.parms[i,6],
                       mainEffectMag = all.parms[i,7],seedVal = all.parms[i,9], modCount = m)
    out2$paramEst <- outFixed[out2$transType,"expEstimate"]
    out2$lowerEst <- outFixed[out2$transType,"expLower"]
    out2$upperEst <- outFixed[out2$transType,"expUpper"]
    ## Now add the shape value
    out2Shape <- out2[1,]
    out2Shape$transType <- "shapeParam"
    shapeEst <- summary(all.mods[[m]])$spec_pars 
    out2Shape$paramEst <- shapeEst[,1]
    out2Shape$lowerEst <- shapeEst$`l-95% CI`
    out2Shape$upperEst <- shapeEst$`u-95% CI`
    out2 <- rbind(out2, out2Shape)
    ## Now see if we have a random effect and add this to the parameter estimates
    all.out <- rbind(all.out, out2)
  }
  ## Now attach the real values
  real.vals <- tmp.dat$transVals
  real.vals$transType <- paste("transType", real.vals$Var1, real.vals$Var2, sep='')
  real.vals <- real.vals[,c("transType", "scale")]
  all.out <- merge(all.out, real.vals, by=c("transType"), all.x=TRUE, suffixes = c("", "_True"))
  ## Now write all.out
  write.csv(all.out, file = out.file, quote = FALSE, row.names = FALSE)
}else{
  print("Job done")
}