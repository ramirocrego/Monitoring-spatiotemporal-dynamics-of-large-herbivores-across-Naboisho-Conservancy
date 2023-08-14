###########################################################################
# Monitoring spatiotemporal dynamics of large herbivores across an African
# rangeland using hierarchical multi-species distance sampling

# Ramiro D. Crego, Harry B.M. Wells, Grant Connette, Jared A. Stabach,
# Naitareu Soit & Stewart Thompson

## Goodness of fit
library(dplyr)

#Load posterior output for parameters
#Abundance
load("./Data/MultiSpDSdata")
load("./ModelOutput/HMSDS_MC.Rdata")
load("./Data/SurveyMonths.Rdata")

Time = as.Date(paste(year_months,'01',sep='_'), format='%Y_%m_%d')

speciesList <- c("Hartebeest",'Eland','Giraffe',"Grant's gazelle",'Impala',"Thomson's gazelle",'Topi','Warthog','Wildebeest',"Plains zebra")

# Extract model values
ntimes = 30
ntransects = 8
nspecies = 10

sp <-  seq(0,(ntimes*ntransects*nspecies), (ntimes * ntransects))
num <- seq(0,(ntimes*ntransects),ntimes)

# Detection prob
pcapcol <- grepl(colnames(mc), pattern = "pcap")
pcap <- mc[,pcapcol]

postPcap <- array(NA, dim = c(nrow(mc), ntimes, ntransects, nspecies))

for(s in 1:nspecies){
  temp <- pcap[, (1+sp[s]):(sp[s+1])]
  for(i in 1:ntransects){
    postPcap[,,i,s] <- temp[,(1+num[i]):(num[i+1])]
  }}

# Lambda star
lambda.starcol <- grepl(colnames(mc), pattern = "lambda.star")
lambda.star <- mc[,lambda.starcol]

postlambdastar <- array(NA, dim = c(nrow(mc), ntimes, ntransects, nspecies))

for(s in 1:nspecies){
  temp <- lambda.star[, (1+sp[s]):(sp[s+1])]
  for(i in 1:ntransects){
    postlambdastar[,,i,s] <- temp[,(1+num[i]):(num[i+1])]
  }}


###################################################################
### Compare observed num groups with simulated number of groups ###
###################################################################

SimNumGroups <- data.frame()
for(iter in 1:nrow(mc)){
  temp2 <- data.frame()
  for(s in 1:nspecies){
    
    yTemp <- array(NA, dim = c(ntimes,ntransects))
    for(i in 1:ntimes){
      for(j in 1:ntransects){
        n <- rpois(1, postlambdastar[iter,i,j,s])
        yTemp[i,j] <- rbinom(1,n,postPcap[iter,i,j,s])
      }}
    Sum1<-apply(yTemp, 1, sum)
    temp<- data.frame(Ngr = Sum1, Time = 1:30, Species = speciesList[s], Sim = iter)
    temp2 <- rbind(temp2, temp)
  }
  SimNumGroups <- rbind(SimNumGroups, temp2)
}
head(SimNumGroups)

resultssim <- SimNumGroups %>% group_by(Species, Time) %>% dplyr::summarise(mean(Ngr), SD = sd(Ngr))  %>% as.data.frame()
resultssim
colnames(resultssim) <- c("Species","Month","Mean", "SD")

obserdata <- data.frame()
data <- DSdata$y
for(s in 1:nspecies){
  temp <- apply(data[,,s], 1, sum)
  temp2 <- data.frame(Species = speciesList[s], Time = 1:30, ObNgr = temp)
  obserdata <- rbind(obserdata, temp2)
}


obserdata <- arrange(obserdata, Species)

FinalTable <- cbind(resultssim,obserdata)
FinalTable$Delta <- FinalTable$ObNgr - FinalTable$Mean
hist(FinalTable$Delta)
mean(FinalTable$Delta)
range(FinalTable$Delta)
write.csv(FinalTable, "./ModelOutPut/SimulatedObservedData.csv")


#################################
### Goodness of fit Abundance ###
#################################

# N
Ncol <- grepl(colnames(mc), pattern = "N")
N <- mc[,Ncol]

postN <- array(NA, dim = c(nrow(mc), ntimes, ntransects, nspecies))

for(s in 1:nspecies){
  temp <- N[, (1+sp[s]):(sp[s+1])]
  for(i in 1:ntransects){
    postN[,,i,s] <- temp[,(1+num[i]):(num[i+1])]
  }}

sumResObs1 <- sumResSim1 <- array(NA, dim = c(nrow(mc), nspecies))
for(iter in 1:nrow(mc)){
  for(s in 1:nspecies){
    resSim <- array(NA, dim = c(ntimes,ntransects))
    resObs <- array(NA, dim = c(ntimes,ntransects))
    for(i in 1:ntimes){
      for(j in 1:ntransects){
        resObs[i,j] <-  (sqrt(postN[iter,i,j,s]) - sqrt(postlambdastar[iter,i,j,s]))^2 # Freeman-Tukey statistic for observed data
        nSim <- rpois(1, postlambdastar[iter,i,j,s])
        resSim[i,j] <-  sum((sqrt(nSim) - sqrt(postlambdastar[iter,i,j,s]))^2) # Freeman-Tukey statistic for simulated data
      }}
    sumResObs1[iter,s] <- sum(resObs)
    sumResSim1[iter,s] <- sum(resSim)
  }
}

FTobs <- apply(sumResObs1, 1, sum)
FTsim <- apply(sumResSim1, 1, sum)

MASS::eqscplot(FTobs, FTsim, xlim=range(FTobs, FTsim), ylim=range(FTobs, FTsim),
               xlab="Observed data", ylab="Simulated data")
abline(0, 1, lwd=2, col='red')
mean(FTobs > FTsim) # the P value
