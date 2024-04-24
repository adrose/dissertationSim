## Load library(s)
library(brms)
library(doParallel)
library(tidyverse)

## CP all data over from schooner
#system("rsync -au --info=progress2 arosen@schooner.oscer.ou.edu:/scratch/arosen/dissertationSim/data/individualSimsMM_SMM/ /home/arosen/Documents/dissertationSim/data/individualSimsMM_SMM/")
#system("rsync -au --info=progress2 arosen@schooner.oscer.ou.edu:/scratch/arosen/dissertationSim/data/individualSimsMM_SMM-3/ /home/arosen/Documents/dissertationSim/data/individualSimsMM_SMM-3/")
#system("rsync -au --info=progress2 arosen@schooner.oscer.ou.edu:/scratch/arosen/dissertationSim/data/individualSimsMM_SMM-3/ /media/arosen/adonsBU/dissertationSimData/individualSimsMM_SMM-3/")

## Now declare all of the simulation components
n <- c(10, 100)
minObsAll <- c(10, 40)
n.states <- c(3)
matrixType <- c("mod", "rand")
scaleVals <- c(".5:8","15:25")
shapeVals <- c(1, 7)
mainEffectVals <- c(0,2)
rand.var <- c(0.1,4)
iter.vals <- 1:1000
all.parms <- expand.grid(n, minObsAll, n.states, matrixType, scaleVals, shapeVals, mainEffectVals,rand.var,iter.vals)
## Now go through and load all of the data
all.files <- list.files("/scratch/arosen/dissertationSimData/individualSimsMM_SMM-3/", pattern = "*RDS",full.names = TRUE, recursive = TRUE)

## Run this in parallel
totalCore <- detectCores()
cl <- makeCluster(totalCore)
registerDoParallel(cl)
target <- foreach(i = all.files) %dopar%{
  tmp.in <- readRDS(i)
  ## First grab the shape parameters that were estimated
  tmp.in1 <- dplyr::bind_rows(lapply(tmp.in[1:2], function(x) summary(x)$spec_pars))
  tmp.in <- dplyr::bind_rows(lapply(tmp.in[1:4], function(x) summary(x)$fixed["effectOfInt",]))
  ## Now find the correct anova params
  rowVals <-  strSplitMatrixReturn(basename(i), "_")[,2]
  tmp.in <- dplyr::bind_cols(tmp.in, all.parms[rowVals,])
  tmp.in$modCount <- 1:4
  tmp.in$shapeEst <- c(tmp.in1$Estimate, 1, 1)
  tmp.in$shapeEstLC <- c(tmp.in1$`l-95% CI`, 0, 0)
  tmp.in$shapeEstUC <- c(tmp.in1$`u-95% CI`, 0, 0)
  ## Return these values
  tmp.in
}
stopCluster(cl)

saveRDS(target, file = "./data/targetLoadVals.RDS")