## This script will be used to translate weibull parameters into semi markov processes 
## It closely follows this paper: https://www.degruyter.com/document/doi/10.1515/ijb-2020-0083/html
## Specifically, sections 2 & Approach 1: sojourn times
## Sepcfically, I will need to translate the survivial times into probabilities of state transitions

## clear cache
rm(list=ls())

## Load library(s)
library(brms)
library(msm)
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
source("./scripts/weiFuncs.R")

## Now lets use these to simulate some tests
nstate <- 3
ntrans <- expand.grid(1:nstate, 1:nstate)
ntrans$transType <- paste(ntrans[,1], ntrans[,2])
## Now create a shape and scale variable for each of these
shape.val <- runif(1, min=.5, max=5)
scale.val <- runif(nrow(ntrans), min = 2, max = 8)
ntrans$shape <- shape.val
ntrans$scale <- scale.val
prob.vals <- diff(c(0, sort(runif(nstate-1)), 1))
prob.vals <- c(.33, .33, .34)
state.trans.vals <- sample(1:nstate, replace = TRUE, size = 501, prob = prob.vals)
state.trans.valsPat <- dplyr::lag(state.trans.vals)
sim.dat <- data.frame(transTo = state.trans.vals, transFrom = state.trans.valsPat,time = NA)
sim.dat$transType <- paste(sim.dat$transFrom, sim.dat$transTo)
sim.dat$order <- 1:nrow(sim.dat)
sim.dat <- merge(sim.dat, ntrans)
sim.dat <- sim.dat[order(sim.dat$order),]
trueScale <- matrix(scale.val, ncol=nstate, nrow=nstate)
trueShape <- shape.val
## Now go through each of these and simulate a weibull time
for(i in 1:nrow(sim.dat)){
  sim.dat$time[i] <- rweibull(1, shape = sim.dat$shape[i], scale = sim.dat$scale[i])
}
sim.dat$timeR <- cumsum(sim.dat$time)
sim.dat <- na.omit(sim.dat)
mod.estR <- brm(time ~ 1  + transType, data = sim.dat, family = weibull(), cores=4)
allParam(mod.estR)
## Now see if we can't recreate our prob of trans values
for(i in seq(0,3,.1)){
  probMat <- allParam(mod.estR, t = i)$probMat
}

