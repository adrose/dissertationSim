## This script will be used to estimate weibull models under various population characteristics
## and then I will examine parameter estimation characteristics under each of these

nAll <- c(10,50,100)
minObsAll <- c(10, 50, 100)
maxObsAll <- c(50, 100, 200)
samplePercents <- expand.grid(c(0, .25, .5, .75, 1),c(0, .25, .5, .75, 1), c(0, .25, .5, .75, 1))
## Now go through each one of these, and grab the data files across each of the smapling permutations
for(i in 1:length(nAll)){
  in.dat <- list.files(path = "./data/", pattern = paste("sampleSize_", nAll[i], "_min", sep=""), full.names = TRUE)
  in.dat <- read.csv(in.dat)

}
