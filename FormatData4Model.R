###########################################################################
# Monitoring spatiotemporal dynamics of large herbivores across an African
# rangeland using hierarchical multi-species distance sampling

# Ramiro D. Crego, Harry B.M. Wells, Grant Connette, Jared A. Stabach,
# Naitareu Soit & Stewart Thompson


#########################
## Format data for model

# Load libraries

library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)

# Load data
load("./Data/DistData")

# Summary
ddply(data, ~ TransectID + Species, dplyr::summarise, NrGropus = n(), NrInd = sum(GroupSize))


# Create a transect ID column
transectnames <- unique(arrange(data, TransectID)[,4])
ntransects <- length(transectnames)
transid <- seq(1,ntransects,1)

data$Site_ID <- data$TransectID
for (i in 1:length(transectnames)){
  data$Site_ID[data$Site_ID == transectnames[i]] <- transid[i]
}
data$Site_ID <- as.integer(data$Site_ID)

# Create a transect rep column
transreps <- unique(arrange(data, YearMonth)[,7])
nreps <- length(transreps)
sampleid <- seq(1,nreps,1)


data$Sample_ID <- data$YearMonth
for (i in 1:length(transreps)){
  data$Sample_ID[data$Sample_ID == transreps[i]] <- sampleid[i]
}
data$Sample_ID <- as.integer(data$Sample_ID)

# Summary
ddply(data, ~ Site_ID + Sample_ID + Species, dplyr::summarise, NrGropus = n())

##### 
data <- droplevels(data)

# Table 1
table1 <- ddply(data, ~ Species, dplyr::summarise, NrGropus = n(), abu = sum(GroupSize), meanGS = mean(GroupSize), MinGS = min(GroupSize), maxGS = max(GroupSize))
write.csv(table1, file = "./ModelOutput/Table1.csv")

##################################################
#Create observation array (rep x site x species)

#Number of sites
nsites <- as.integer(max(data$Site_ID))

#Number of reps
nreps <-  as.integer(max(data$Sample_ID))

#Number of Species
nspec <- as.integer(length(unique(data$Species)))

#Initialize observation array (rep x site x species)
y <- array(NA, dim = c(nreps,nsites,nspec))

#Vector of species
name <- unique(data$Species)

#Generate observation array
for(s in 1:nspec){
  J <- (filter(data, Species == name[s]))
  J <- group_by(J,  Site_ID, Sample_ID, Species)%>%dplyr::summarize(n())
  W <- data.frame(rep(1:nsites, rep(nreps, nsites)), rep(1:nreps, nsites))
  colnames(W) <- c("Site_ID", "Sample_ID")
  Y <- full_join(W, J, by = c("Site_ID", "Sample_ID"))
  Y$`n()`[is.na(Y$`n()`)] = 0
  X <- split(Y$`n()`, f = Y$Site_ID)
  X <- do.call(cbind, X)
  y[,,s] <- X
}

# Create distance classes
distdata <- data %>% drop_na(Distance)
hist(distdata$Distance)

# Width of distance classes
v <- 50 #meters

# Transect half-width
B <- 300 #meters


# Distance classes
dclass <- distdata$Distance %/% v + 1

# Check for observation bias. ie. larger groups are seen more frequently at farther distances
distdata$dcalss <- dclass
ddply(distdata, ~ Species + dclass, dplyr::summarise, MeanGS = mean(GroupSize), SD = sd(GroupSize))
ggplot(distdata, aes(x = dclass, y = GroupSize)) + facet_wrap(~Species, scales="free_y") + geom_point() + geom_smooth(method="lm")

# Check groups and distance histogram for all species
forhist <- ddply(distdata, ~ Species + dclass, dplyr::summarise, N = n())
ggplot(forhist, aes(x = dclass, y = N)) + geom_col() + facet_wrap(~Species, scales="free_y")

##########################
# Create IDs and indices

#Species ID
spec <- as.integer(as.factor(droplevels(distdata$Species)))

#Site ID
site <- distdata$Site_ID

#Replicate ID
rep <- distdata$Sample_ID

#Distance class midpoint ID
mdpt <- seq(25, B, v)

#Number of observations
nobs <- sum(y, na.rm = TRUE)

#Number of distance classes
nD <- length(mdpt)

# Import group size of observations
gs <- distdata$GroupSize

# Body mass
traits <- read.csv('./Data/Traits_GEC.csv')
bodymass <- traits$Mass[1:10]

BMsc <- scale(bodymass)

################
# Combine data

DSdata <- list(y, dclass, v, B, gs, site, rep, spec, mdpt, nsites, nreps, nspec, nobs, nD, BMsc)
heads <- c("y", "dclass", "v", "B", "gs", "site", "rep", "spec", "mdpt", "nsites","nreps", "nspec", "nobs", "nD", "BMsc")
DSdata <- setNames(DSdata, nm = heads)

# Save data
save(DSdata, file = "./Data/MultiSpDSdata")

