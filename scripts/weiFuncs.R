### Quick script to understand the hazard function of the weibull distribution

## I will need to define several functions -- perhaps all will just
## call the built in weibull distribution tools
## but I will need a survival function, the cumulative distribution function
## which requires the probability distribution function (pweibull)


## First grab the cumulative distribution function here
## this is the integral of the probability distribution function from 0 --> t

cdfWei <- function(time=1, shape=1, scale=1){
  cdfOut <- pweibull(time, shape = shape, scale = scale, lower.tail = TRUE)
  return(cdfOut)
}

## now do the survival function
survWei <- function(time, shape, scale){
  survOut <- 1 - cdfWei(time, shape, scale)
  return(survOut)
}

## Now do the hazard function
## This is the cumulative hazard 
hazardWei <- function(time, shape, scale){
  #hazardOut <- shape / scale * (time / scale)^shape - 1
  hazardOut <- -pweibull(q = time, shape = shape, scale = scale, log.p = TRUE, lower.tail = FALSE)
  return(hazardOut)
}

## Here is the instantaneous hazard
instHazardWei <- function(time, shape, scale){
  ## Take the ratio of hazard & surv function
  outHazardInst <- ifelse(time < 0, 0, shape * (time/scale)^(shape - 1)/scale)
  return(outHazardInst)
}

## First declare some functions
estSurvivialFunc <- function(time = 1, scale.val = 1, shape.val = 1){
  survVal <- exp( - (time / scale.val)^shape.val )
  return(survVal)
}


estIntent <- function(time =1, scale.val = 1, shape.val = 1, trans.prob = .5){
  ## The output of this function should be the transition intensity for
  ## a transfer from state i --> j
  ## I think I want to use the following formula:
  ## p_{ij} * exp(- (t/scale) ^ shape ) /
  new.scale <- trans.prob ^ (-1/shape.val) * scale.val
  out.int <- shape.val / new.scale * ((time / new.scale) ^ (scale.val -1))
  return(out.int)
}

## I want to make a function which will perform all o fthis -- I want to provide it with
## the output brms model, the time that I want to estimate the transition intensities
allParam <- function(brmsMod, timeVal=1){
  ## This function will take the estimated weibull brm model
  ## and return several outputs, these include:
  ## 1. discrete transition patterns
  ## 2. transition intensity @ time(t)
  ## 3. probability of transition at time(t)
  ## 4. Static probability transition matrix (integrate out t)
  ## 5. model parameters
  
  
  ### Begin by identifying the transition values
  ### This will return output #1
  tmp.dat <- brmsMod$data
  from <- strSplitMatrixReturn(tmp.dat[,2], " ")[,1]
  to <- strSplitMatrixReturn(tmp.dat[,2], " ")[,2]
  trans.mat1 <- table(from, to)
  trans.mat <- trans.mat1
  diag(trans.mat) <- 0
  trans.names <- reshape2::melt(trans.mat1)
  trans.names <- paste(trans.names$from, trans.names$to)
  trans.matP <- trans.mat / rowSums(trans.mat)
  trans.matPR <- trans.mat1 / rowSums(trans.mat1)
  ### Output # 4: model parameters
  tmp.dat <- data.frame(transType = factor(trans.names),
                        timeIn = 0)
  n.state <- sqrt(dim(tmp.dat)[1])
  tmp.dat <- model.matrix(timeIn ~ 1 + transType, data = tmp.dat)
  ## now multiply these with the fixef from the model
  scale.vals <- tmp.dat %*% fixef(brmsMod)[,1]
  scale.vals <- exp(scale.vals) #/ gamma(1 + 1/scale.vals)
  scale.vals <- matrix(scale.vals, nrow=n.state, ncol=n.state)
  shape.val <- summary(brmsMod)$spec_pars[,1]
  
  ## Now work on output #2
  int.mat <- matrix(0, ncol=n.state, nrow=n.state)
  for(i in 1:n.state){ ## row
    for(j in 1:n.state){ ## col
      new.scale.val <- trans.matP[i,j]^(-1/shape.val) * scale.vals[i,j]
      int.mat[i,j] <- trans.matP[i,j]*hazardWei(time = timeVal, shape = shape.val, scale = scale.vals[i,j])
    }
  }
  diag(int.mat) <- -1 * rowSums(int.mat, na.rm = TRUE)
  # diag(int.mat2) <- -1 * rowSums(int.mat2, na.rm = TRUE)
  
  ## Now work on output #3
  p.mat <- matrix(0, ncol=n.state, nrow=n.state)
  for(i in 1:n.state){
    for(j in 1:n.state){
      if(i == j){
        sumVal <- 0
        iterVals <- (1:n.state)[-i]
        for(acrossVals in iterVals){
          sumValTmp <- trans.matP[i,acrossVals] * estSurvivialFunc(time = timeVal, scale.val = scale.vals[i,acrossVals],shape = shape.val)
          sumVal <- sumVal + sumValTmp
        }
        p.mat[i,j] <- sumVal
      }
    }
  }
  ## Now do the off diagonal values
  for(i in 1:n.state){
    for(j in 1:n.state){
      if(i != j){
        p.mat[i,j] <- (1 - p.mat[i,i]) * int.mat[i,j] / -int.mat[i,i]
      }
    }
  }
  ## Now work on static transition rates
  ######## WORK ON THIS SECTION
  # It is formula 23 from the Estimation of semi-Markov multi-state models manuscript
  stat.trans <- matrix(0, nrow=n.state, ncol=n.state)
  for(i in 1:n.state){
    for(j in 1:n.state){
      if(i !=j ){
        col.index <- 1:ncol(stat.trans)
        col.index <- col.index[-i]
        denom.val <- sum(scale.vals[i,col.index]^shape.val)
        stat.trans[i,j] <- (scale.vals[i,j]^shape.val) / denom.val
      }
    }
  }
  
  ## Now return all of these values
  all.out <- list(weibullScales = scale.vals, weibullShape = shape.val,
                  transInt = int.mat, probMat = p.mat, observedTransWO = trans.matP, 
                  estimatedTrans = stat.trans)
  return(all.out)
}


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

## Now make a function which will return a markov chain matrix
## given the number of states, and if the matrix should be 
## random in nature (all states have equal probability of being chose)
## stable in nature (pii is the most likley state to be chose)
## or moderatly stable (pii is more than twice as likley to be chose)
retMMMat <- function(nstate, matPat = "stable"){
  if(!matPat %in% c("stable", "mod", "rand") & is.character(matPat)){
    errorCondition(paste("matPat must be one of: stable, mod, rand"))
  }
  ## Now determine our values
  if(matPat == "stable"){
    transIn <- .9
    transOut <- (1 - transIn) / (nstate -1)
  }
  if(matPat == "mod"){
    transIn <- .5
    transOut <- (1 - transIn) / (nstate -1)
  }
  if(matPat == "rand"){
    transIn <- transOut <- 1 / nstate
  }
  if(is.numeric(matPat)){
    transIn <- matPat
    transOut <- (1 - transIn) / (nstate -1)
  }
  ## Now assign the matrix
  outputMat <- matrix(transOut, nrow=nstate, ncol=nstate)
  diag(outputMat) <- transIn
  ## Now return matrix
  return(outputMat)
}


## Now create a function which will return the sampling priors given the fixed effects and 
## desired random effect variance from a weibull model
retSampPriors <- function(fixEfModel, modShape = 1.6, groupName = "part", ranVarInt = .3, ranVarSlope = 0, addMainEffect=NULL){
  ## First create all of the strings we need
  ## start with the fixed effects
  all.fix.effect.vals <- NULL
  for(i in 1:nrow(fixEfModel)){
    coef.name <- rownames(fixEfModel)[i]
    coef.val <- fixEfModel[i,"Estimate"]
    ## Now create the prior string
    if(i == 1){
      prior.string <- prior("constant(XYZ)", class="ABC")
      prior.string$prior <- gsub(pattern = "XYZ", replacement = coef.val, x = prior.string$prior)
      prior.string$class <- coef.name
      
    }else{
      prior.string <- prior("constant(XYZ)", coef="ABC")
      prior.string$prior <- gsub(pattern = "XYZ", replacement = coef.val, x = prior.string$prior)
      prior.string$coef <- coef.name
    }
    ## Now modify
    all.fix.effect.vals <- rbind(all.fix.effect.vals, prior.string)
  }
  ## Add model shape
  prior.string <- prior("constant(XYZ)", class="shape", lb="0")
  prior.string$prior <- gsub(pattern = "XYZ", replacement = modShape, x = prior.string$prior)
  all.fix.effect.vals <- rbind(all.fix.effect.vals, prior.string)
  ## Now do the random effects
  all.mix.effect.vals <- NULL
  for(i in 1:nrow(fixEfModel)){
    #Start with the variance of the intercept
    if(i == 1){
      ## Begin with initializing the cor across random effects
      tmp.p1 <- prior("lkj(1)", class="cor", group="subject")
      tmp.p1$group <- groupName
      
      ## Now do the intercept
      tmp.p2 <- prior("constant(XYZ)", class="sd", group="ABC", coef="Intercept")
      tmp.p2$prior <- gsub(pattern = "XYZ", replacement = ranVarInt, x = tmp.p2$prior)
      tmp.p2$group <- groupName
      all.mix.effect.vals <- rbind(all.mix.effect.vals, tmp.p1)
      all.mix.effect.vals <- rbind(all.mix.effect.vals, tmp.p2)
    }else{
      tmp.p1 <- prior("constant(XYZ)", class="sd", group="ABC", coef="DEF")
      tmp.p1$prior <- gsub(pattern = "XYZ", replacement = ranVarSlope, x = tmp.p1$prior)
      tmp.p1$group <- gsub(pattern = "ABC", replacement = groupName, x = tmp.p1$group)
      tmp.p1$coef <- gsub(pattern = "ABC", replacement = groupName, x = rownames(fixEfModel)[i])
      all.mix.effect.vals <- rbind(all.mix.effect.vals, tmp.p1)
    }
  }
  out.priors <- rbind(all.fix.effect.vals, all.mix.effect.vals)
}


## Now create a function which will run through a set of these params and estimate everything we need
simFunc <- function(n=10,minObs=20, maxObs=90, nState=3, 
                    matrixType="stable", scaleRange="2:6", shapeVal=1.6, meMag = NULL, addLag = FALSE){
  ## First load all required packages
  source("./scripts/weiFuncs.R")
  ## Now create our data
  # first create our number of sample vectors
  obsPat <- sample(minObs:maxObs, size = n,replace = TRUE)
  ## Now create a data frame to house all of our variables
  samp.dat <- data.frame(part = rep(1:n, times = obsPat), order = NA,timeIn = NA,stateFrom = NA,stateTo=NA, stateFrom2=NA, stateFrom3=NA)
  ## Check for ME value
  if(!is.null(meMag)){
    samp.dat$effectOfInt <- rep(rnorm(n), times = obsPat)
  }
  ## We have to get our discrete MM first
  discreteMMVals <- retMMMat(nstate = nState, matPat = matrixType)
  ## Now go through each participant and estimate their discrete markov matrix
  for(i in 1:n){
    samp.dat[which(samp.dat$part==i),"stateTo"] <- run.mc.sim(discreteMMVals, num.iters = obsPat[i])
    ## Now identify our lagged pattern
    samp.dat[which(samp.dat$part==i),"stateFrom"] <- dplyr::lag(samp.dat[which(samp.dat$part==i),"stateTo"])
    samp.dat[which(samp.dat$part==i),"stateFrom2"] <- dplyr::lag(samp.dat[which(samp.dat$part==i),"stateTo"], n = 2)
    samp.dat[which(samp.dat$part==i),"stateFrom3"] <- dplyr::lag(samp.dat[which(samp.dat$part==i),"stateTo"], n = 3)
    
    ## Now fix the NA value
    samp.dat[which(samp.dat$part==i)[1],"stateFrom"] <- sample(1:nState, size = 1)
    ## Now fix the order variable
    samp.dat[which(samp.dat$part==i),"order"] <- 1:obsPat[i]
  }
  ## Now paste these together
  samp.dat$transType <- paste(samp.dat$stateFrom, samp.dat$stateTo)
  
  ## Now create our transition type values
  if(addLag){
    ntrans <- expand.grid(1:nState, 1:nState, 1:nState)
    ntrans$transType <- paste(ntrans[,1], ntrans[,2], ntrans[,3])
    samp.dat$transType <- paste(samp.dat$stateFrom2, samp.dat$stateFrom, samp.dat$stateTo)
  }else{
    ntrans <- expand.grid(1:nState, 1:nState)
    ntrans$transType <- paste(ntrans[,1], ntrans[,2])
  }
  ## Now create a shape and scale variable for each of these
  shape.val <- shapeVal
  scaleMin <- as.numeric(strSplitMatrixReturn(charactersToSplit = scaleRange, ":")[,1])
  scaleMax <- as.numeric(strSplitMatrixReturn(charactersToSplit = scaleRange, ":")[,2])
  scale.val <- runif(nrow(ntrans), min = scaleMin, max = scaleMax)
  ntrans$shape <- shape.val
  ntrans$scale <- scale.val
  samp.dat <- merge(samp.dat, ntrans, by=c("transType"))
  samp.dat <- samp.dat[order(samp.dat$part, samp.dat$order),]
  ## Now go through each of these and sample our weibull times
  if(! is.null(meMag)){
    for(i in 1:nrow(samp.dat)){
      samp.dat$timeIn[i] <- rweibull(1, shape = samp.dat$shape[i], scale = samp.dat$scale[i]+samp.dat$effectOfInt*meMag)
    }
  }else{
    for(i in 1:nrow(samp.dat)){
      samp.dat$timeIn[i] <- rweibull(1, shape = samp.dat$shape[i], scale = samp.dat$scale[i])
    }
  }  
  ## Now return everything we need
  outList <- list(sampleData = samp.dat, transVals = ntrans)
  return(outList)
}
