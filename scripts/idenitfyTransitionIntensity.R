## Begin by loading packages and setting seed
library(msm)
library(brms)
library(foreach)
library(doParallel)
library(dplyr)
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
## First simulate a time sreires that is long enough to know the population fixed effects
## with a known intensitie transition matrix -- I will make three of these
## one with equal prob of trnasitioning -- one with moderate stationarity -- and one with very stable transition patterns
qmatrix.Rand <- rbind(c(-0.6,   0.3,  0.3 ),
                      c(0.3,   -0.6,  0.3 ),
                      c(0.3,   0.3, -0.6))
## Now simulate the data for a single subject
sim.dfR <- data.frame(subject = 1, timeR = runif(201, min = 0, max = 10))
sim.dfR$time <- cumsum(sim.dfR$timeR)
sim.vals <- simmulti.msm(sim.dfR, qmatrix.Rand)
## Now create our transition type varaible
sim.vals$lagState <- dplyr::lag(sim.vals$state)
sim.vals$transType <- paste(sim.vals$lagState, sim.vals$state, sep="")
sim.vals$timeIn <- dplyr::lag(sim.vals$time)
## Now remove intrastate transitions
#sim.vals <- sim.vals[-which(sim.vals$transType %in% c("11", "22", "33")),]
sim.vals$timeIn <- c(NA, diff(sim.vals$time))
sim.valsR <- na.omit(sim.vals)

## Now calculate time spent in each observation
mod.estR <- brm(timeIn ~ 1  + transType, data = sim.valsR, family = weibull(), cores=4)
msm.estR <- msm(state ~ time, data = sim.vals, qmatrix = qmatrix.Rand, subject = subject, 
                exacttimes = TRUE,control=list(fnscale=4000))

## Ok I now need to translate these paramweter estimates into varoius outcomes of a usrivial model -- such as
## density, surivial function & hazard rate -- see thiese here: https://www.degruyter.com/document/doi/10.1515/ijb-2020-0083/html
## First start with identifying the scale parameters for each transition type
## start these by declaring a model matrix with intercept for all trans types
tmp.dat <- data.frame(transType = factor(names(table(sim.valsR$transType))),
                  timeIn = 0)
tmp.dat <- model.matrix(timeIn ~ 1 + transType, data = tmp.dat)
## now multiply these with the fixef from the model
scale.vals <- tmp.dat %*% fixef(mod.estR)[,1]
scale.vals <- exp(scale.vals) / gamma(1 + 1/scale.vals)
## NOw grab the shape value
shape.val <- summary(mod.estR)$spec_pars[,1]

## Now we have the shape and scale values for each transition
## start with density function here
densVals <- dweibull(seq(0,8,.1), shape = shape.val, scale = scale.vals[1])
## now do survivial
survVals <- exp( - (seq(.1,8,.1) / scale.vals[1])^shape.val )
## Now do the hazard rate
hazardRate <- scale.vals[1] / shape.val * (seq(.1,8,.1) / shape.val[1])^(scale.vals[1]-1)
## Now try to identify the transition intensity values
## Again this is taken from: https://www.degruyter.com/document/doi/10.1515/ijb-2020-0083/html
#Ok this is going to have to be done in some kind of loop 
## FIrst estimate our intensities
## I am going to write a function to perofrm this
## First create a function which will idenitfy the size of state matrix when given the off diagonal
## scale parameters
retMatDim <- function(scaleParams){
  dimGuess <- 2
  doneVal <- TRUE
  while(doneVal){
    ## First create a matrix with the dim of the dimGuess
    tmp.mat <- matrix(NA, nrow=dimGuess, ncol = dimGuess)
    ## Now idenitfy the length of the off diagonal
    off.diag.length <- length(tmp.mat[row(tmp.mat)!=col(tmp.mat)])
    if(off.diag.length==length(scaleParams)){
      doneVal <- FALSE
    }else{
      dimGuess <- dimGuess + 1
    }
  }
  return(dimGuess)
}

retIntMat <- function(shapeParam = 1, scaleParams = runif(6), t=1, withinT=FALSE){
  ## First prepare the output matrix
  if(withinT){
    out.mat.dim <- sqrt(length(scaleParams))
  }else{
    out.mat.dim <- retMatDim(scaleParams = scaleParams)
  }
  proto.mat <- matrix(0, ncol = out.mat.dim, nrow=out.mat.dim)
  ## Now assign the scale params to their respective elements in the proto mat
  scale.mat <- proto.mat
  if(withinT){
    scale.mat <- matrix(scaleParams, nrow = out.mat.dim, ncol=out.mat.dim, byrow = FALSE)
  }else{
    scale.mat[row(scale.mat)!=col(scale.mat)] <- scaleParams
  }
  ## Now go through every off diag element and estimate probability of transition
  intentMat <- proto.mat
  for(i in 1:nrow(scale.mat)){
    for(j in 1:ncol(scale.mat)){
      ## First estimate the state specific intensity
      intentMat[i,j] <-  scale.mat[i,j] * shapeParam * (t^(scale.mat[i,j]-1))
    }
  }
  ## Now make the diagonal the sum of the off diagonal within each row
  diag(intentMat) <- -1 * rowSums(intentMat, na.rm = TRUE)
  return(intentMat)
}




## Now do the probability values
retProbMat <- function(shapeParam=1, scaleParams = runif(6), withinT=FALSE){
  ## First prepare the output matrix
  if(withinT){
    out.mat.dim <- sqrt(length(scaleParams))
  }else{
    out.mat.dim <- retMatDim(scaleParams = scaleParams)
  }
  proto.mat <- matrix(0, ncol = out.mat.dim, nrow=out.mat.dim)
  ## Now prep our scale values
  scale.mat <- proto.mat
  if(withinT){
    scale.mat <- matrix(scaleParams, nrow = out.mat.dim, ncol=out.mat.dim, byrow = FALSE)
  }else{
    scale.mat[row(scale.mat)!=col(scale.mat)] <- scaleParams
  }
  ## Now run through each value and idenitfy our prob of transition
  prob.mat <- proto.mat
  for(i in 1:nrow(prob.mat)){
    for(j in 1:ncol(prob.mat)){
      if(i !=j){
        ## First estimate the state specific intensity
        prob.mat[i,j] <-  scale.mat[i,j] ^ shapeParam
      }
    }
  }
  ## Now take the sum for all elements where i!=j
  prob.mat2 <- prob.mat
  for(i in 1:nrow(prob.mat)){
    for(j in 1:ncol(prob.mat)){
      if(i != j){
        prob.mat2[i,j] <- prob.mat[i,j] / rowSums(prob.mat)[i]
      }
    }
  }
  return(prob.mat2)
}

## Now estimate our intensity matrix
intMat <- retIntMat(shapeParam = shape.val, scaleParams = scale.vals, t=1, withinT = TRUE)
proMat <- retProbMat(shapeParam = shape.val, scaleParams = scale.vals, withinT = TRUE)


## Now try to create our weibull times from a weibull distirbution sampling procedure
nstate <- 3
ntrans <- expand.grid(1:nstate, 1:nstate)
ntrans$transType <- paste(ntrans[,1], ntrans[,2])
## Now create a shape and scale variable for each of these
shape.val <- runif(1, min=1.5, max=3)
scale.val <- runif(nrow(ntrans), min = 2, max = 8)
ntrans$shape <- shape.val
ntrans$scale <- scale.val
q.mat <- retIntMat(shapeParam = shape.val, scaleParams = scale.val, withinT = TRUE)
p.mat <- retProbMat(shapeParam = shape.val, scaleParams = scale.val, withinT = TRUE)
## Now sample from each of these params
state.trans.vals <- sample(1:nstate, replace = TRUE, size = 5001, prob = diff(c(0, sort(runif(nstate-1)), 1)))
state.trans.valsPat <- dplyr::lag(state.trans.vals)
sim.dat <- data.frame(transTo = state.trans.vals, transFrom = state.trans.valsPat,time = NA)
sim.dat$transType <- paste(sim.dat$transFrom, sim.dat$transTo)
sim.dat$order <- 1:nrow(sim.dat)
sim.dat <- merge(sim.dat, ntrans)
sim.dat <- sim.dat[order(sim.dat$order),]
## Now go through each of these and simulate a weibull time
for(i in 1:nrow(sim.dat)){
  sim.dat$time[i] <- rweibull(1, shape = sim.dat$shape[i], scale = sim.dat$scale[i])
}
sim.dat$timeR <- cumsum(sim.dat$time)
sim.dat <- na.omit(sim.dat)
## Now train model
mod.estR <- brm(time ~ 1  + transType, data = sim.dat, family = weibull(), cores=4)
msm.estR <- msm(transTo ~ timeR, data = sim.dat, qmatrix = qmatrix.Rand,
                exacttimes = TRUE,control=list(fnscale=4000))
## Now estimate our scale vals 
tmp.dat <- data.frame(transType = factor(names(table(sim.dat$transType))),
                      timeIn = 0)
tmp.dat <- model.matrix(timeIn ~ 1 + transType, data = tmp.dat)
## now multiply these with the fixef from the model
scale.vals <- tmp.dat %*% fixef(mod.estR)[,1]
scale.vals <- exp(scale.vals) / gamma(1 + 1/scale.vals)
shape.val <- summary(mod.estR)$spec_pars[,1]
## Now do our p mat
q.mat2 <- retIntMat(shapeParam = shape.val, scaleParams = scale.val, withinT = TRUE, t = 1)
p.mat2 <- retIntMat(shapeParam = shape.val, scaleParams = scale.vals, withinT = TRUE)

## Now compare these estimates
pmatrix.msm(msm.estR, t = 1)
-q.mat2[1,3]/q.mat2[1,1]
