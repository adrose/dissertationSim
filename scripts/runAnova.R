## Load library(s)
library(brms)
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
source("./scripts/weiFuncs.R")

## Now declare all of the simulation components
n <- c(10, 100)
minObsAll <- c(20, 80)
maxObsAll <- c(90, 180)
n.states <- c(3,4)
matrixType <- c("mod", "rand")
scaleVals <- c("2:7", "5:10")
shapeVals <- c(1.6, 2.7)
all.parms <- expand.grid(n, minObsAll, maxObsAll, n.states, matrixType, scaleVals, shapeVals)

## Now go through and load all of the data
all.files <- list.files("./data/individualSims/", full.names = TRUE)

## Now go through and add all of these data to one another
target <- read.csv(all.files[1])

for(i in all.files[-c(1)]){
  tmp.in <- read.csv(i)
  target <- dplyr::bind_rows(target, tmp.in)
}

## Now estimate our param estimation error
target$error <- target$paramEst - target$scale
## Now do root mean squared error
target$rse <- sqrt((target$scale - target$paramEst)^2)

## Now model this error
target$n <- factor(target$n)
target$minObs <- factor(target$minObs)
target$maxObs <- factor(target$maxObs)
target$nState <- factor(target$nState)
target$shapeVal <- factor(target$shapeVal)
## Now chnage the name of the intercept variable
target$transType <- plyr::revalue(target$transType, c("Intercept" = "transType11"))
## Now change the randVar value
target$randVarChar <- target$randVar
target$randVarChar[target$randVar==0] <- "None"
target$randVarChar[which(target$randVar==.05 & target$transType=="transType11")] <- "Mod"
target$randVarChar[which(target$randVar==.1 & target$transType=="transType11")] <- "Large"
target$randVarChar[which(target$randVar==.01 & target$transType!="transType11")] <- "Mod"
target$randVarChar[which(target$randVar==.05 & target$transType!="transType11")] <- "Large"

# mod.one <- lm(error ~ (randVarChar+ n + minObs + maxObs + nState + matrixType + scaleRange + shapeVal)^2, data = target)
# anova.one <- car::Anova(mod.one)
# effectsize::eta_squared(anova.one, partial = FALSE)
mod.two <- lm(rse ~ -1 + (randVarChar+ n + minObs + maxObs + nState + matrixType + scaleRange + shapeVal), data = target)
anova.two <- car::Anova(mod.two)
ef.vals <- effectsize::eta_squared(anova.two, partial = TRUE)
effectsize::eta_squared(anova.two, partial = TRUE)
## Now visualize this?
library(ggplot2)
visreg::visreg(mod.two, "randVarChar", gg=TRUE) + theme_bw() + xlab("Random_Variance") + ylab("RMSE")

# Now plot these eta squared values
for.plot <- data.frame(ef.vals)
for.plot <- for.plot[order(for.plot$Eta2_partial, decreasing = TRUE),]
for.plot$rank <- 1:nrow(for.plot)
for.plot$CharVal <- plyr::revalue(for.plot$Parameter, c("randVarChar" = "Random_Variance",
                                                        "scaleRange" = "Scale_Range",
                                                        "n" = "n",
                                                        "nState" = "n_States",
                                                        "shapeVal" = "Shape_Param",
                                                        "maxObs" = "Max_Obs_Length",
                                                        "minObs" = "Min_Obs_Length",
                                                        "matrixType" = "Matrix_Type"))
library(ggplot2)
library(ggrepel)
ggplot(for.plot, aes(x=rank, y=Eta2_partial)) +
  geom_point() +
  geom_label_repel(aes(label=CharVal)) +
  theme_bw() +
  ylab("Eta2") +
  xlab("Rank")
