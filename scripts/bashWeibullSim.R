### Now make an R script which can be called by the shell and will run one row 
### from the simulateWeibul.R script -- I think this should help to avoid the following error:
### Error in dyn.load(libLFile) : 
### unable to load shared object '/tmp/RtmphUpD1w/filebef1be2acc8.so':
###  `maximal number of DLLs reached...

## First idenitfy the row of the simulation we want to run
i=as.numeric(commandArgs(TRUE)[1])
seedVal=as.numeric(commandArgs(TRUE)[2])
set.seed(seedVal)
## Now declare output file
out.file <- paste("./data/individualSims/rowVal_", i, "_seedVal_", seedVal, ".csv")
if(!file.exists(out.file)){
  
  ## Load library(s)
  library(brms)
  source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
  source("./scripts/weiFuncs.R")
  
  ## Now declare all of the simulation components
  n <- c(10, 50, 100)
  minObsAll <- c(20, 40, 80)
  maxObsAll <- c(90, 140, 180)
  n.states <- c(2,3,4)
  matrixType <- c("stable", "mod", "rand")
  scaleVals <- c("2:6", "4:8")
  shapeVals <- c(1.6, 2.2)
  all.parms <- expand.grid(n, minObsAll, maxObsAll, n.states, matrixType, scaleVals, shapeVals)
  rand.var.int <- c(.05, .1)
  rand.var.slope <- c(.01, .05)
  ## Now go ahead and loop through all of these params real quick
  mod.count <- 1
  ## First create the data
  tmp.dat <- simFunc(n = all.parms[i,1], minObs = all.parms[i,2],maxObs = all.parms[i,3],nState = all.parms[i,4], 
                     matrixType = all.parms[i,5], scaleRange = all.parms[i,6], shapeVal = all.parms[i,7])
  mod.est <- brm(timeIn ~ 1  + transType, data = tmp.dat$sampleData, family = weibull(), cores=1,
                 control = list(adapt_delta = 0.9, max_treedepth = 12))
  for.sim <- data.frame(fixef(mod.est))
  for.sim$transType <- rownames(for.sim)
  for.sim$transType[for.sim$transType=="Intercept"] <- "transType11"
  new.vec <- rep(0, nrow(for.sim))
  for.sim.tmp <- tmp.dat$transVals
  for.sim.tmp$transType <- paste("transType", for.sim.tmp$Var1, for.sim.tmp$Var2, sep='')
  ## Now merge the two
  for.sim <- merge(for.sim, for.sim.tmp)
  for.sim$Estimate <- for.sim$scale
  all.mods <- list()
  all.mods[[mod.count]] <- mod.est
  mod.count <- mod.count + 1
  for(r in 1:length(rand.var.int)){
    priorVar <- retSampPriors(fixEfModel = fixef(mod.est), modShape = tmp.dat$transVals$shape[1], groupName = "part",
                              ranVarInt = rand.var.int[r], ranVarSlope = rand.var.slope[r])
    sampGen <- brm(timeIn ~ 1 + transType + (1 + transType|part), data = tmp.dat$sampleData, family = weibull(), 
                   cores=1, prior = priorVar,sample_prior = "only", seed = 16, iter = 10, chains = 1)
    tmp.data <- posterior_predict(sampGen)
    tmp.dat$sampleData$genVals <- tmp.data[1,]
    mod.est2 <- brm(genVals ~ 1  + transType, data = tmp.dat$sampleData, family = weibull(), cores=1,
                    control = list(adapt_delta = 0.9, max_treedepth = 12))
    all.mods[[mod.count]] <- mod.est2
    mod.count <- mod.count + 1
  }
  ## Now prep all of the output
  ## Now loop through all models
  rand.var.count <- 1
  all.out <- NULL
  for(m in 1:length(all.mods)){
    outFixed <- data.frame(summary(all.mods[[m]])$fixed)
    outFixed$expEstimate <- NA
    outFixed$expLower <- NA
    outFixed$expUpper <- NA
    ## Put these back in the original units
    for(r in 1:nrow(outFixed)){
      if(r == 1){
        outFixed$expEstimate[r] <- exp(outFixed$Estimate[r])
        outFixed$expLower[r] <- exp(outFixed$l.95..CI[r])
        outFixed$expUpper[r] <- exp(outFixed$u.95..CI[r])
      }else{
        outFixed$expEstimate[r] <- exp(sum(outFixed$Estimate[c(1,r)]))
        outFixed$expLower[r] <- exp(sum(outFixed$l.95..CI[c(1,r)]))
        outFixed$expUpper[r] <- exp(sum(outFixed$u.95..CI[c(1,r)]))
      }
    }
    ## Now prepare all output
    out2 <- data.frame(transType=c(rownames(outFixed)), paramEst = NA, lowerEst = NA, upperEst = NA,
                       randVar = 0, rowAllParam = i, n = all.parms[i,1], minObs = all.parms[i,2], maxObs = all.parms[i,3],
                       nState = all.parms[i,4], matrixType = all.parms[i,5], scaleRange = all.parms[i,6], shapeVal = all.parms[i,7],
                       allParmsRow = i)
    out2$paramEst <- outFixed[out2$transType,"expEstimate"]
    out2$lowerEst <- outFixed[out2$transType,"expLower"]
    out2$upperEst <- outFixed[out2$transType,"expUpper"]
    ## Now see if we have a random effect and add this to the parameter estimates
    if(m > 1){
      out2$randVar <- rand.var.slope[rand.var.count]
      out2$randVar[out2$transType=="Intercept"] <- rand.var.int[rand.var.count]
      rand.var.count <- rand.var.count + 1
    }
    all.out <- rbind(all.out, out2)
  }
  
  ## Now write all.out
  write.csv(all.out, file = out.file, quote = FALSE, row.names = FALSE)
}else{
  print("Job done")
}