## Clean slate
rm(list = ls())

## Load library(s)
library("ggplot2")
library("depmixS4")
library("mHMMbayes")
library("tidyverse")
source("~/Documents/Github/adroseHelperScripts/R/afgrHelpFunc.R")
library("reshape2")
library("msm")
library("visreg")
library("brms")

## Load data
more.dat1 <- foreign::read.spss("~/Documents/oregonDPICS/data/W1_rawDataEpochNums.sav", to.data.frame = TRUE)
more.dat2 <- foreign::read.spss("~/Documents/oregonDPICS/data/W2_rawDataEpochNums.sav", to.data.frame = TRUE)
more.dat3 <- foreign::read.spss("~/Documents/oregonDPICS/data/W3_rawDataEpochNums.sav", to.data.frame = TRUE)
caps.data <- readRDS("~/Documents/eegOregon/data/capsData.RDS")

## Clean the .sav files
more.dat1$Behavior <- str_squish(more.dat1$Behavior)
more.dat2$Behavior <- str_squish(more.dat2$Behavior)
more.dat3$Behavior <- str_squish(more.dat3$Behavior)
more.dat1$Behavior <- gsub(pattern = " ", replacement = "_", x = more.dat1$Behavior)
more.dat2$Behavior <- gsub(pattern = " ", replacement = "_", x = more.dat2$Behavior)
more.dat3$Behavior <- gsub(pattern = " ", replacement = "_", x = more.dat3$Behavior)
## Now make the noncompliable commands a distinct flavor
more.dat1$Behavior[which(more.dat1$Modifier_1=="No Opportunity")] <- paste(more.dat1$Behavior[which(more.dat1$Modifier_1=="No Opportunity")], "_NO", sep='') 
more.dat2$Behavior[which(more.dat2$Modifier_1=="No Opportunity")] <- paste(more.dat2$Behavior[which(more.dat2$Modifier_1=="No Opportunity")], "_NO", sep='') 
more.dat3$Behavior[which(more.dat3$Modifier_1=="No Opportunity")] <- paste(more.dat3$Behavior[which(more.dat3$Modifier_1=="No Opportunity")], "_NO", sep='') 


## Now try to organize data, so the columns are the actions that can be taken from the parent (indirect command, neutral talk.. et cetera), the rows are still the time, but it is a binary indiciator for when an event occurs
wide.form1.T <- reshape2::dcast(more.dat1, Observation + Time_Relative_sf ~ Behavior, fun.aggregate = length)
## Some values == 2?
## Observation CAPS181-CLUP-H2 ## Time: 292.08600
more.dat1[more.dat1$Observation=="CAPS181-CLUP-H2" & more.dat1$Time_Relative_sf =="292.086",]
#       Observation Time_Relative_sf Time_Relative_hmsf       Behavior     Modifier_1  FAMILY DPICS_Task EPOCH
# 31874 CAPS181-CLUP-H2          292.086              292.1 Direct_Command No Opportunity CAPS181       CLUP    10
# 31875 CAPS181-CLUP-H2          292.086              292.1 Direct_Command No Opportunity CAPS181       CLUP    10
## Looks like I do have identical rows -- I am just going to take all unique rows
wide.form1.T <- reshape2::dcast(unique(more.dat1), Observation + Time_Relative_sf ~ Behavior, fun.aggregate = length)
wide.form2.T <- reshape2::dcast(unique(more.dat2), Observation + Time_Relative_sf ~ Behavior, fun.aggregate = length)

## Still seeing values > 1 with the unique fix for values from 2
## Observation: CAPS195-CLUP-H2 & Time_Relative_sf: 5.57212
more.dat2[more.dat2$Observation=="CAPS195-CLUP-H2" & more.dat2$Time_Relative_sf=="5.57212",]
#     Observation Time_Relative_sf Time_Relative_hmsf       Behavior     Modifier_1  FAMILY DPICS_Task EPOCH
# 8270 CAPS195-CLUP-H2          5.57212                5.6 Direct_Command Comply         CAPS195       CLUP     1
# 8426 CAPS195-CLUP-H2          5.57212                5.6 Direct_Command Comply         CAPS194       CLUP     1
## Looks like the family got messed up
df.2.index <- duplicated(more.dat2[,c("Observation", "Time_Relative_sf", "Behavior", "FAMILY")])
wide.form2.T <- reshape2::dcast(unique(more.dat2), Observation + Time_Relative_sf ~ Behavior, fun.aggregate = length)

wide.form3.T <- reshape2::dcast(unique(more.dat3), Observation + Time_Relative_sf ~ Behavior, fun.aggregate = length)

## go through and correct the 2 values
iter.vals <- colnames(wide.form2.T)[3:15]
for(i in iter.vals){
  length.fix <- sum(wide.form2.T[,i] > 1)
  if(length.fix > 0){
    wide.form2.T[which(wide.form2.T[,i]>1),i] <- 1
  }
}


## Identify the block of the study here
wide.form1.T$Block <- NA
wide.form1.T$Block[grep(pattern = "CLED", x = wide.form1.T$Observation)] <- "CLED"
wide.form1.T$Block[grep(pattern = "PLED", x = wide.form1.T$Observation)] <- "PLED"
wide.form1.T$Block[grep(pattern = "PELD", x = wide.form1.T$Observation)] <- "PLED"
wide.form1.T$Block[grep(pattern = "CLUP", x = wide.form1.T$Observation)] <- "CLUP"
table(wide.form1.T$Block, useNA = "ifany")
wide.form2.T$Block <- NA
wide.form2.T$Block[grep(pattern = "CLED", x = wide.form2.T$Observation)] <- "CLED"
wide.form2.T$Block[grep(pattern = "PLED", x = wide.form2.T$Observation)] <- "PLED"
wide.form2.T$Block[grep(pattern = "CLUP", x = wide.form2.T$Observation)] <- "CLUP"
table(wide.form2.T$Block, useNA = "ifany")
wide.form3.T$Block <- NA
wide.form3.T$Block[grep(pattern = "CLED", x = wide.form3.T$Observation)] <- "CLED"
wide.form3.T$Block[grep(pattern = "CLEAD", x = wide.form3.T$Observation)] <- "CLED"
wide.form3.T$Block[grep(pattern = "PLED", x = wide.form3.T$Observation)] <- "PLED"
wide.form3.T$Block[grep(pattern = "CLUP", x = wide.form3.T$Observation)] <- "CLUP"
table(wide.form3.T$Block, useNA = "ifany")

## Now see if I can do it by participant
wide.form1.T$participant <- str_extract(wide.form1.T$Observation, "CAPS[0-9][0-9][0-9]")
wide.form2.T$participant <- str_extract(wide.form2.T$Observation, "CAPS[0-9][0-9][0-9]")
wide.form3.T$participant <- str_extract(wide.form3.T$Observation, "CAPS[0-9][0-9][0-9]")

## Now identify session code
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
#First remove all blank characters in the strings
wide.form1.T$Observation <- str_trim(wide.form1.T$Observation)
## Now remove the Block and the participant from the Observation string
rm.string <- function(x, rmVals = c("CAPS[0-9][0-9][0-9]", "CLED", "PLED", "PELD", "CLUP", "CLEAD", "--", "-")){
  out.vals <- x
  for(i in rmVals){
    out.vals <- gsub(pattern = i, x = out.vals, replacement = "")
  }
  return(out.vals)
}
wide.form1.T$sessionCode <- rm.string(wide.form1.T$Observation)
wide.form2.T$sessionCode <- rm.string(wide.form2.T$Observation)
wide.form3.T$sessionCode <- rm.string(wide.form3.T$Observation)

## Now melt and plot
plot.dat1 <- melt(wide.form1.T, id.vars = c("Time_Relative_sf", "Block", "Observation", "participant", "sessionCode", "Var.3"))
plot.dat2 <- melt(wide.form2.T, id.vars = c("Time_Relative_sf", "Block", "Observation", "participant", "sessionCode", "Var.3"))
plot.dat3 <- melt(wide.form3.T, id.vars = c("Time_Relative_sf", "Block", "Observation", "participant", "sessionCode", "Var.3"))
plot.dat <- bind_rows(plot.dat1, plot.dat2, plot.dat3)
## Now collapse variables into the 3 class system: PRIDE, Dont, Neutral
collpaseVars <- list()
collpaseVars[[1]] <- c("Behavior_Description", "Labeled_Praise", 
                       "Reflection", "Unlabeled_Praise") ## Pride
collpaseVars[[2]] <- c("Neutral_Talk") ## Neutral
collpaseVars[[3]] <- c("Direct_Command", "Direct_Command_NO", 
                       "Indirect_Command", "Indirect_Command_NO",
                       "Negative_Talk", "Question", "Command") ## Don't
collpaseVars[[4]] <- c("Behavior_Description", "Labeled_Praise", 
                       "Reflection", "Unlabeled_Praise",
                       "Direct_Command", "Indirect_Command", "Command") # rowSums(all.counts[,c("LP", "UP","Re", "BD", "DC", "IC")])
collpaseVars[[5]] <- c("Neutral_Talk", "Question")
collpaseVars[[6]] <- c("Direct_Command_NO", "Indirect_Command_NO","Negative_Talk") # rowSums(all.counts[,c("NeT", "DCn","ICn")])

## Now do this for all participants
wide.form1.T$Wave <- 1
wide.form3.T$Wave <- 3
for.msm <- merge(wide.form1.T, wide.form3.T, all=TRUE)
for.msm$state <- NA
for.msm$diffTime <- NA
for.msm$participant <- paste(for.msm$participant, for.msm$Wave, sep="_")
## Trun all "Commands" into a 0 value
for.msm$Command[is.na(for.msm$Command)] <- 0
for(t in unique(for.msm$participant)){
  for(i in 1:6){
    if(i %in% c(1:3)){
      use.vars <- collpaseVars[[i]]
      index <- which(for.msm$Block=="CLED" & for.msm$participant==t)
      ## Now calculate the diff within participant
      for.msm$diffTime[index] <- c(NA, diff(for.msm$Time_Relative_sf[index]))
      for(z in index){
        for(u in use.vars){
          if(for.msm[z,u]==1){for.msm$state[z] <- i}
        }
      }
    }
    if(i %in% c(4:6)){
      use.vars <- collpaseVars[[i]]
      index <- which(for.msm$Block!="CLED" & for.msm$participant==t)
      for.msm$diffTime[index] <- c(NA, diff(for.msm$Time_Relative_sf[index]))
      for(z in index){
        for(u in use.vars){
          if(for.msm[z,u]==1){for.msm$state[z] <- i-3}
        }
      }
    }
  }
}

## Now create a matrix with all of the times we need and an ability to include state
out.states <- NULL
check.val <- 1000
for(i in unique(for.msm$participant)){
  for(b in unique(for.msm$Block)){
    index <- which(for.msm$participant==i & for.msm$Block==b)
    MX <- NULL
    count.val <- length(index) + 1
    for(z in 1:count.val){
      #print(z)
      if(z == 1){
        tmp.row <- c(i,b,0,4,0)
        MX <- rbind(MX, tmp.row)
      }else{
        back.row <- z -1
        tmp.time <- as.numeric(for.msm[index[back.row],c("Time_Relative_sf")])
        tmp.time.diff <- tmp.time - as.numeric(MX[back.row,3])
        ## Now identify tmp.state
        if(b=="CLED"){
          for(p in 1:3){
            sum.val <- sum(for.msm[index[back.row],collpaseVars[[p]]])
            if(sum.val>0){tmp.state <- p}
          }
        }else{
          for(p in 4:6){
            sum.val <- sum(for.msm[index[back.row],collpaseVars[[p]]])
            if(sum.val>0){tmp.state <- p-3}
          }
        }
        tmp.row <- NULL
        if(tmp.time.diff>check.val & MX[back.row,4] !=4){
          insert.val <- as.numeric(MX[back.row,3]) + 2
          tmp.out <- rbind(c(i,b,insert.val,4), c(i,b,tmp.time,tmp.state))
        }else{
          tmp.out <- c(i,b,tmp.time,tmp.state,tmp.time.diff)
        }
        MX <- rbind(MX, tmp.out)
      }
    }
    out.states <- rbind(out.states, MX)
  }
}

## Clean up outstates
out.states <- data.frame(out.states)
out.states$X3 <- as.numeric(out.states$X3)
out.states$X4 <- as.numeric(out.states$X4)
out.states$X5 <- as.numeric(out.states$X5)
out.states$Part <- paste(out.states$X1, out.states$X2)
out.states$Block <- as.numeric(factor(out.states$X2, levels=c("CLED", "PLED", "CLUP")))

## Now add the disengaged state here
id.rep <- unique(out.states$X1)
block.rep <- 1:3
unengage.time <- 7
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}
out.states.UN <- NULL
for(d in id.rep){
  #print(d)
  for(b in block.rep){
    #print(b)
    tmp.data <- out.states[which(out.states$X1 == d & out.states$Block == b),]
    diff.times <- c(0,diff(tmp.data$X3))
    runDiff <- FALSE
    if(length(which(diff.times>unengage.time))>0){
      runDiff <- TRUE
    }
    while(runDiff){
      index.one <- which(diff.times>unengage.time)[1]
      tmp.row <- index.one - 1
      new.row <- tmp.data[tmp.row,]
      new.row$X3 <- new.row$X3 + (unengage.time-.01)
      new.row$X4 <- 4
      tmp.data <- insertRow(tmp.data, new.row, index.one)
      diff.times <- c(0, diff(tmp.data$X3))
      if(length(which(diff.times>unengage.time))>0){
        runDiff <- TRUE
      }else{
        runDiff <- FALSE
      }
    }
    out.states.UN <- bind_rows(out.states.UN, tmp.data)  
  }
}

## Now organize the counts
out.states2 <- out.states.UN %>% group_by(X1, Block) %>% 
  arrange(X1, Block,X3)
## Drop subjects where count is only 1
mod.dat <- out.states2 %>% group_by(Part) %>% 
  dplyr::mutate(countVal = n()) %>% 
  filter(countVal > 12) %>% 
  filter(X3 < 300)
## Now drop observations with identical time
mod.dat <- mod.dat[-which(duplicated(mod.dat[,c("X1", "X2", "X3", "Part")])),]
## Now make it one continous time series wsithin each individual
mod.dat$OrigTime <- mod.dat$X3
mod.dat$X3[mod.dat$X2=="PLED"] <- mod.dat$X3[mod.dat$X2=="PLED"] + 300
mod.dat$X3[mod.dat$X2=="CLUP"] <- mod.dat$X3[mod.dat$X2=="CLUP"] + 600

## Plot a single partiicpant's time series
mod.dat[which(mod.dat$X1=="CAPS111_1"),] %>% ggplot(., aes(x=X3, y=X4)) +
  #geom_line() +
  geom_point() + 
  facet_grid(. ~ Block, scales = "free")

## Now estimate time spent in each state
state.time <- NULL
for(d in id.rep){
  tmp.out <- data.frame(id=d, block=1:3, time_1=NA, time_2 = NA, time_3 = NA, time_4 = NA)
  for(b in 1:3){
    tmp.dat <- mod.dat[which(mod.dat$X1==d & mod.dat$Block==b),]
    ## Identify all values
    all.trans.vals <- cbind(tmp.dat$X4, c(0,diff(tmp.dat$X3)), cumsum(c(0, diff(tmp.dat$X3))))
    ## Go through all values and add the time wihtin each state
    for(s in 1:4){
      index.vals <- which(all.trans.vals[,1]==s)
      if(length(index.vals)>0){
        tot.count <- 0
        for(i in index.vals){
          if(i+1 <= nrow(all.trans.vals)){
            tot.count <- tot.count + diff(c(all.trans.vals[i,3],all.trans.vals[i+1,3]))
          }else{
            tot.count <- tot.count + diff(c(all.trans.vals[i,3],300))
          }
        }
        tmp.col <- paste("time_", s, sep='')
        tmp.out[b,tmp.col] <- tot.count
      }
    }
  }
  state.time <- bind_rows(state.time, tmp.out)
}
## Add a wave effect for the state time
state.time$wave <- strSplitMatrixReturn(state.time$id, "_")[,2]
state.time$part <- strSplitMatrixReturn(state.time$id, "_")[,1]
findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}
plot.dat <- melt(state.time, id.vars = c("id","block","wave", "part")) %>% 
  group_by(variable) %>% 
  filter(!is.na(value)) %>% 
  mutate(value = as.numeric(value)) %>% 
  mutate(outlier = ifelse(findoutlier(value), part, NA))

## Now save these times to disk
state.time$wave <- paste("Wave_", state.time$wave, sep='')
state.time$block <- paste("Block_", state.time$block, sep='')
out.vals <- tidyr::pivot_wider(state.time[,-1], id_cols = c("part"), names_from = c("block","wave"), values_from = c("time_1", "time_2", "time_3", "time_4"))

## Now attach the child behavior to each of these and use this as a time variant covars
wide.form1.C <- reshape2::dcast(more.dat1, Observation + Time_Relative_sf ~ Modifier_1, fun.aggregate = length)
wide.form3.C <- reshape2::dcast(more.dat3, Observation + Time_Relative_sf ~ Modifier_1, fun.aggregate = length)
## Now see if we can merge these into the states dataset?
wide.form1.C <- wide.form1.C[,c(1,2,4,6,7)]
names(wide.form1.C) <- gsub(patter=" ", replacement = "", x = names(wide.form1.C))
wide.form3.C <- wide.form3.C[,c(1,2,4,6,7)]
names(wide.form3.C) <- gsub(patter=" ", replacement = "", x = names(wide.form3.C))
## Now merge these onto the state dataset?
wide.form1.C$Block <- NA
wide.form1.C$Block[grep(pattern = "CLED", x = wide.form1.C$Observation)] <- "CLED"
wide.form1.C$Block[grep(pattern = "PLED", x = wide.form1.C$Observation)] <- "PLED"
wide.form1.C$Block[grep(pattern = "PELD", x = wide.form1.C$Observation)] <- "PLED"
wide.form1.C$Block[grep(pattern = "CLUP", x = wide.form1.C$Observation)] <- "CLUP"
wide.form3.C$Block <- NA
wide.form3.C$Block[grep(pattern = "CLED", x = wide.form3.C$Observation)] <- "CLED"
wide.form3.C$Block[grep(pattern = "CLEAD", x =wide.form3.C$Observation)] <- "CLED"
wide.form3.C$Block[grep(pattern = "PLED", x = wide.form3.C$Observation)] <- "PLED"
wide.form3.C$Block[grep(pattern = "CLUP", x = wide.form3.C$Observation)] <- "CLUP"

## Now see if I can do it by participant
wide.form1.C$participant <- str_extract(wide.form1.C$Observation, "CAPS[0-9][0-9][0-9]")
wide.form3.C$participant <- str_extract(wide.form3.C$Observation, "CAPS[0-9][0-9][0-9]")

## Now merge these?
wide.form1.C$Time_Relative_sf <- round(wide.form1.C$Time_Relative_sf, 2)
wide.form1.C$Wave <- 1
wide.form3.C$Time_Relative_sf <- round(wide.form3.C$Time_Relative_sf, 2)
wide.form3.C$Wave <- 3
mod.dat <- data.frame(mod.dat)
mod.dat$timeMerge <- as.numeric(round(mod.dat$OrigTime, 2))
mod.dat$Wave <- strSplitMatrixReturn(mod.dat$X1, "_")[,2]
## Now merge these into foo
wide.form1.C2 <- merge(wide.form1.C, wide.form3.C, all=TRUE)
mod.dat$Part <- strSplitMatrixReturn(mod.dat$X1, "_")[,1]
mod.dat$MergeCol <- paste(mod.dat$Part, mod.dat$X2, mod.dat$timeMerge, mod.dat$Wave)
wide.form1.C2$MergeCol <- paste(wide.form1.C2$participant, wide.form1.C2$Block, wide.form1.C2$Time_Relative_sf, wide.form1.C2$Wave)
mod.dat2 <- merge(mod.dat, wide.form1.C2, by.x=c("Part", "X2", "timeMerge", "Wave"), by.y=c("participant", "Block", "Time_Relative_sf", "Wave"), all.x = TRUE)

## Now convert NA values --> 0
mod.dat2$Comply[is.na(mod.dat2$Comply)] <- 0
mod.dat2$Comply[mod.dat2$Comply==2] <- 1
mod.dat2$NoOpportunity[is.na(mod.dat2$NoOpportunity)] <- 0
mod.dat2$NoOpportunity[mod.dat2$NoOpportunity==2] <- 1
mod.dat2$Noncomply[is.na(mod.dat2$Noncomply)] <- 0
mod.dat2$Noncomply[mod.dat2$Noncomply==2] <- 1

mod.dat2 <- mod.dat2 %>% group_by(Part, Block) %>% 
  arrange(X1, Block,X3)

mod.dat2$Wave[mod.dat2$Wave==1] <- 1
mod.dat2$Wave[mod.dat2$Wave==3] <- 2
brms.dat <- NULL
for(d in id.rep){
  #print(d)
  for(b in block.rep){
    tmp.data <- mod.dat2[which(mod.dat2$X1 == d & mod.dat2$Block == b),]
    ## Remove any duplicated time elements
    ## Trigger duplicated warning here
    if(sum(duplicated(tmp.data$X3))>0){
      tmp.data <- tmp.data[-which(duplicated(tmp.data$X3)),]
    }
    if(dim(tmp.data)[1]>0){
      ## Now cycle through every row and identify the stay length, current state, next state, and previous state
      out.dat <- NULL
      for(z in 2:nrow(tmp.data)){
        row.prev <- z - 1
        row.next <- z + 1
        stayLength <- tmp.data$X3[z] - tmp.data$X3[row.prev]
        currentState <- tmp.data$X4[z]
        nextState <- tmp.data$X4[row.next]
        if(is.na(nextState)){nextState <- currentState}
        transType <- paste(currentState, nextState)
        statePrev <- tmp.data$X4[row.prev]
        if(z > 2){
          timePrev <- out.dat[nrow(out.dat),5]
        }else{
          timePrev <- 0.5
        }
        ## Now grab childs behavior
        complyPrev <- tmp.data$Comply[row.prev]
        complyCurrent <- tmp.data$Comply[z]
        noOpPrev <- tmp.data$NoOpportunity[row.prev]
        noOpCur <- tmp.data$NoOpportunity[z]
        nonComplyPrev <- tmp.data$Noncomply[row.prev]
        nonComplyCurrent <- tmp.data$Noncomply[z]
        out.row <- c(d,b, row.prev, row.next, stayLength, currentState, nextState, transType, statePrev, timePrev,
                     complyPrev, complyCurrent, noOpPrev, noOpCur, nonComplyPrev, nonComplyCurrent)
        out.dat <- rbind(out.dat, out.row)
      }
      brms.dat <- rbind(brms.dat, out.dat)
    }
  }
}
## Make everything the correct data type
brms.dat <- data.frame(brms.dat)
brms.dat$Block <- factor(brms.dat$V2,levels = c(1,2,3))
brms.dat$V5 <- as.numeric(brms.dat$V5)
brms.dat$stayLength <- as.numeric(brms.dat$V5)
## Drop all state 4 values
brms.dat <- brms.dat[-which(brms.dat$V6==4),]
brms.dat <- brms.dat[-which(brms.dat$V7==4),]
#brms.dat <- brms.dat[-which(brms.dat$V8==4),]
brms.dat <- brms.dat[-which(brms.dat$V9==4),]
brms.dat$currentState <- factor(brms.dat$V6, levels = c(1,2,3))
brms.dat$nextState <- factor(brms.dat$V7, levels = c(1,2,3))
brms.dat$transType <- factor(brms.dat$V8)
brms.dat$statePrev <- factor(brms.dat$V9, levels=c(1,2,3))
brms.dat$timePrev <- as.numeric(brms.dat$out.row)
brms.dat$subject <- strSplitMatrixReturn(brms.dat$V1, "_")[,1]
brms.dat$wave <- strSplitMatrixReturn(brms.dat$V1, "_")[,2]
brms.dat$childBehavior <- "NoBehavior"
for(i in 1:nrow(brms.dat)){
  if(brms.dat$V12[i]==1){brms.dat$childBehavior[i] <- "Comply"}
  #if(brms.dat$V14[i]==1){brms.dat$childBehavior[i] <- "NoOp"}
  if(brms.dat$V16[i]==1){brms.dat$childBehavior[i] <- "nonComply"}
}
## Fix the edge cases
brms.dat$childBehavior[which(brms.dat$childBehavior=="Comply" & brms.dat$transType=="2 3")] <- "NoBehavior"
brms.dat$childBehavior[which(brms.dat$childBehavior=="Comply" & brms.dat$transType=="2 4")] <- "NoBehavior"
brms.datAll <- brms.dat[which(brms.dat$V2==3),]
brms.datAll <- merge(brms.datAll, caps.data, by.x="subject", by.y = "FAMILY", all.x=TRUE)
brms.datmod <- brms.datAll[which(brms.datAll$Block==3 & brms.datAll$subject %in% unique(brms.datAll$subject)),]

priors <- get_prior(stayLength ~ -1 + (transType), data = brms.datmod, family=weibull())
priors2 <- get_prior(stayLength ~ -1 + (transType)  + (1|subject), data = brms.datmod, family=weibull())
## Now drop all state 4 values


mod.est <- brm(stayLength ~ -1 + (transType), data = brms.datmod, family=weibull(), cores = 4) ## The most basic markov model here
## now try the msm model here
priorsMLString <- c(
  #prior("constant(0)", coef=""),
  prior("constant(.2)", coef="transType11"),
  prior("constant(2)", coef="transType12"),
  prior("constant(1)", coef="transType13"),
  prior("constant(1)", coef="transType21"),
  prior("constant(.2)", coef="transType22"),
  prior("constant(2)", coef="transType23"),
  prior("constant(2)", coef="transType31"),
  prior("constant(1)", coef="transType32"),
  prior("constant(.2)", coef="transType33"),
  #prior("constant(1)", coef="", class="sd", group=="", coef=""), 
  #prior("constant(1)", coef="", class="sd", group=="subject", coef=""),
  prior("constant(3)", class="sd", group="subject", coef="Intercept"), 
  prior("constant(1.5)", class="shape", lb="0")## random effect variance here
)

## Now see if we can simulate data with various properties
#mod.est <- brm(stayLength ~ -1 + (transType), data = brms.datmod, family=weibull(), cores = 4)
sampGen <- brm(stayLength ~ -1 + transType + (1|subject), data = brms.datmod, chains = 1,family=weibull())#, prior = priorsMLString,sample_prior = "only")
sampGen <- brm(stayLength ~ -1 + transType + (1|subject), data = brms.datmod, chains = 1,family=weibull(), prior = priorsMLString,sample_prior = "only")
test <- posterior_predict(sampGen)
## Now estimate model
brms.datmod$testOne <- test[1,]
sampTest <- brm(testOne ~ -1 + transType + (1|subject), data = brms.datmod, cores = 4,family=weibull(),
                control = list(adapt_delta = 0.9, max_treedepth = 12), thin = 50, iter = 20000, warmup = 4000)

