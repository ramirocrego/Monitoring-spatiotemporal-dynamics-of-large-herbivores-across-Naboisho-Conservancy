###########################################################################
# Monitoring spatiotemporal dynamics of large herbivores across an African
# rangeland using hierarchical multi-species distance sampling

# Ramiro D. Crego, Harry B.M. Wells, Grant Connette, Jared A. Stabach,
# Naitareu Soit & Stewart Thompson

# Load libraries

library(bayestestR)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
library(dplyr)
library(ggforce)
library(reshape)
library(mcmcOutput)
library(pals) # color palette

# Set ggplot theme

theme_set(theme_bw(base_size = 16))

# Load Files

# Load posterior output for parameters
load("./ModelOutput/HMSDS_MC.Rdata")
# Formatted Data
load("./Data/MultiSpDSdata")
attach(DSdata)
load("./Data/SurveyMonths.Rdata")
Time <- as.Date(paste(year_months,'01',sep='_'), format='%Y_%m_%d')

# Set colors for plots
speciesList <-  c("Hartebeest",'Eland','Giraffe',"Grant's gazelle",'Impala',"Thomson's gazelle",'Topi','Warthog','Wildebeest',"Plains zebra")
col <- c(polychrome(n = length(speciesList)), "#000000")
spelistCom <- c(speciesList, "Community")
names(col) <- spelistCom

##########################
# Detection probability

# Simulate distance across transect half-width
dist.sim <- seq(0, 300, 10)

# Harvest sigma values
mu_g <- grepl(colnames(mc), pattern = "mu_g")
probsmu_g <- mean(mc[,mu_g])
gamma0 <- grepl(colnames(mc), pattern = "^gamma0\\[")
probsG0 <- apply(mc[,gamma0],2,mean)
sigmaR1 <- c(exp(probsG0), exp(probsmu_g))

distfuncR1 <- matrix(NA, nrow = nspec+1, ncol = length(dist.sim))
for(i in 1:(nspec+1)){
  for(j in 1:length(dist.sim)){
    distfuncR1[i,j] <- exp(-dist.sim[j]*dist.sim[j]/(2*sigmaR1[i]*sigmaR1[i]))
  }
}
forplotDist <- melt(distfuncR1)
forplotDist$Species <- rep(spelistCom, length(dist.sim))
forplotDist$Dist <- rep(dist.sim, each = nspec+1)

# Plot detection probability vs observed distance

(FigDetProb <- ggplot(forplotDist, aes(x = Dist, y = value, group = Species, colour = Species)) + geom_line(size = 1.5) + 
    scale_color_manual(name = , values = rev(col)) + 
    labs(x = "Distance (m)", y = "Detection probability") +
    theme(
      legend.key.size = unit(0.2, "in"),
      legend.position = c(0.8,0.66),
      legend.title = element_blank(),
      legend.text = element_text(size = 16))
)

ggsave(file = "./Figures/DetectionProb.jpg", plot = FigDetProb, width = 8, height = 6)

# Body mass effect
median(mc[,"gamma1"])
ci(mc[,"gamma1"], 0.89)
mean(mc[,"gamma1"] > 0)

#########################################
# Figure caterpiller plot for covariates

# NDVI
mua1_rowsP <- grepl(colnames(mc), pattern = "^mu_a1")
alpha1_rowsP <- grepl(colnames(mc), pattern = "^alpha1\\[")
alpha1P <- cbind(mc[,alpha1_rowsP], mc[,mua1_rowsP])
alpha1df <- data.frame()
for(i in 1:(nspec+1)){
  PSpos <- mean(alpha1P[,i] > 0)
  PSneg <- mean(alpha1P[,i] < 0)
  temp <- data.frame(Species = spelistCom[i], Median = median(alpha1P[,i]), ci(alpha1P[,i], 0.89), PS = ifelse(mean(alpha1P[,i]) > 0, PSpos, PSneg), Cov = c('NDVI'))
  alpha1df <- rbind(alpha1df,temp)
}

# NDVI2
mua2_rowsP <- grepl(colnames(mc), pattern = "^mu_a2")
alpha2_rowsP <- grepl(colnames(mc), pattern = "^alpha2\\[")
alpha2P <- cbind(mc[,alpha2_rowsP],mc[,mua2_rowsP])
alpha2df <- data.frame()
for(i in 1:(nspec+1)){
  PSpos <- mean(alpha2P[,i] > 0)
  PSneg <- mean(alpha2P[,i] < 0)
  temp <- data.frame(Species = spelistCom[i], Median = median(alpha2P[,i]), ci(alpha2P[,i], 0.89),  PS = ifelse(mean(alpha2P[,i]) > 0, PSpos, PSneg), Cov = c('NDVI^2'))
  alpha2df <- rbind(alpha2df,temp)
}

# Trend
mua3_rowsP <- grepl(colnames(mc), pattern = "^mu_a3")
alpha3_rowsP <- grepl(colnames(mc), pattern = "^alpha3\\[")
alpha3P <- cbind(mc[,alpha3_rowsP],mc[,mua3_rowsP])
alpha3df <- data.frame()
for(i in 1:(nspec+1)){
  PSpos <- mean(alpha3P[,i] > 0)
  PSneg <- mean(alpha3P[,i] < 0)
  temp <- data.frame(Species = spelistCom[i], Median = median(alpha3P[,i]), ci(alpha3P[,i], 0.89), PS = ifelse(mean(alpha3P[,i]) > 0, PSpos, PSneg), Cov = c('Time'))
  alpha3df <- rbind(alpha3df,temp)
}

postAlpha <- rbind(alpha1df, alpha2df, alpha3df)

# visualise beta means & 89% CIs
(pNG <- ggplot(postAlpha, aes(x=factor(Species, levels = spelistCom), y=Median, colour=PS*100)) +
    geom_hline(yintercept=0,linetype = 'dashed', linewidth=1) + coord_flip() + facet_wrap(~Cov, scales='free_x') +
    geom_point(size=3) + geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0, size=1) +
    scale_colour_gradient2(name='Posterior\nprobability (%)', low='#FDE725FF', mid = '#22A884FF', high='#440154FF', midpoint = 80) + ylab('Parameter estimate') + xlab(NULL) +
    theme(axis.text=element_text(size=11), 
          axis.title=element_text(size=12), axis.line=element_line(linewidth=1),
          axis.ticks=element_line(size=1), axis.ticks.length=unit(2, "mm"),
          strip.text=element_text(size=11,face='bold'), strip.background=element_blank(),
          legend.title=element_text(size=12), legend.text=element_text(size=12), 
          plot.title=element_text(hjust=0.5,face='bold',size=12),
          panel.border=element_blank(), panel.grid=element_blank()) + ggtitle('Estimated number of groups'))

ggsave(file = "./Naboisho/Figures/ParametersNG.jpg", plot = pNG, width = 9, height = 4)


#########################################
# NDVI predictions

# Intercepts
# Alphas
alpha0_P <- grepl(colnames(mc), pattern = "^alpha0\\[")
alpha0P <- apply(mc[,alpha0_P], MARGIN = 2, median) 

num <- seq(0,300, 30)
alpha0P2 <-NULL
for(i in 1:nspec){
  temp <- alpha0P[(1+num[i]):(30+num[i])]
  alpha0P2 <- rbind(alpha0P2, temp)
}

# Overall trend across time
mualpha0_P <- grepl(colnames(mc), pattern = "^mu.alpha0\\[")
mualpha0P <- mc[,mualpha0_P]


# Slopes
# NDVI
alpha1_P <- grepl(colnames(mc), pattern = "^alpha1\\[")
alpha1_Post <- mc[,alpha1_P]
alpha1P <- apply(alpha1_Post, MARGIN = 2, median) 

alpha2_P <- grepl(colnames(mc), pattern = "^alpha2\\[")
alpha2_Post <- mc[,alpha2_P]
alpha2P <- apply(alpha2_Post, MARGIN = 2, median)


# Calculate number of groups vs ndvi
load("./Data/NDVI")
ndvi.simscal <- seq(min(NDVI), max(NDVI), 0.1)
ndvi.sim <- (ndvi.simscal * 0.13469) + 0.4855522 # mean 0.4855522 sd 0.13469

ndvicov <- array(NA, dim = c(nreps,length(ndvi.simscal),nspec))
for(s in 1:nspec){
  for(t in 1:nreps){
    ndvicov[t,,s] <- exp(alpha0P2[s,t] + alpha1P[s] * ndvi.simscal + alpha2P[s] * ndvi.simscal * ndvi.simscal)
  }
}

forplotndvi <- NULL
for(i in 1:nspec){
  for(t in 1:nreps){
    temp <- data.frame(x = ndvi.sim, y = ndvicov[t,,i], time = Time[t], Species = speciesList[i])
    forplotndvi <- rbind(forplotndvi, temp)
  }}


# Overall ndvi across time
ndvicov2 <- array(NA, dim = c(nrow(mc),length(ndvi.sim),nspec))
for(i in 1:nrow(mc)){
  for(s in 1:nspec){
    ndvicov2[i,,s] <- exp(mualpha0P[i,s] + alpha1_Post[i,s] * ndvi.simscal + alpha2_Post[i,s] * ndvi.simscal * ndvi.simscal)
  }
}

forplotndvi2 <- NULL
for (s in 1:nspec){
  for(i in 1:length(ndvi.sim)){
    temp <- data.frame(Species = speciesList[s], y = median(ndvicov2[,i,s]), yUpper = ci(ndvicov2[,i,s], 0.89)$CI_high, yLower = ci(ndvicov2[,i,s], 0.89)$CI_low)
    forplotndvi2 <- rbind(forplotndvi2, temp)
  }}
forplotndvi2$x <- rep(ndvi.sim, nspec)

forplotndvi$Species <- factor(forplotndvi$Species, levels=  c("Hartebeest", 'Eland','Giraffe',"Grant's gazelle",'Impala',"Thomson's gazelle",'Topi','Warthog','Wildebeest',"Plains zebra", "Community"))

forplotndvi2$Species <- factor(forplotndvi2$Species, levels= c("Hartebeest", 'Eland','Giraffe',"Grant's gazelle",'Impala',"Thomson's gazelle",'Topi','Warthog','Wildebeest',"Plains zebra", "Community"))

# Add community response
mua0_P <- grepl(colnames(mc), pattern = "^mu_a0")
mua0_P <- mc[,mua0_P]

mua1_P <- grepl(colnames(mc), pattern = "^mu_a1")
mua1_P <- mc[,mua1_P]
median(mua1_P)
ci(mua1_P, ci = 0.89)
mean(mua1_P > 0)

mua2_P <- grepl(colnames(mc), pattern = "^mu_a2")
mua2_P <- mc[,mua2_P]
median(mua2_P)
ci(mua2_P, ci = 0.89)
mean(mua2_P < 0)

commndvi <- array(NA, dim = c(nrow(mc),length(ndvi.sim)))
for(i in 1:nrow(mc)){
  commndvi[i,] <- exp(mua0_P[i] + mua1_P[i] * ndvi.simscal + mua2_P[i] * ndvi.simscal * ndvi.simscal)
}

commndvi2 <- data.frame()
for(i in 1:length(ndvi.sim)){
  temp <- data.frame(Species = "Community", y = median(commndvi[,i]), yUpper = ci(commndvi[,i], 0.89)$CI_high, yLower = ci(commndvi[,i], 0.89)$CI_low)
  commndvi2 <- rbind(commndvi2, temp)
}
commndvi2$x <- ndvi.sim
commndvi2$Species <- factor(commndvi2$Species, levels= c("Hartebeest", 'Eland','Giraffe',"Grant's gazelle",'Impala',"Thomson's gazelle",'Topi','Warthog','Wildebeest',"Plains zebra", "Community"))

# Plot 
ndvi <- ggplot(forplotndvi, aes(x = x, y = y, group = factor(time))) + facet_wrap(~Species, scales="free_y") + geom_line(colour = 'grey') + geom_ribbon(data = forplotndvi2, aes(ymin = yLower, ymax = yUpper, group = Species), fill = "blue", alpha = 0.5) + geom_line(data = forplotndvi2, aes(x = x, y = y, group = Species), size=1, colour = 'blue') +  geom_ribbon(data = commndvi2, aes(ymin = yLower, ymax = yUpper, group = Species), fill = 'black', alpha = 0.5) + geom_line(data = commndvi2, aes(x= x, y = y, group = Species), size=1, colour = 'black') +
  labs(y=expression("Nr. of groups"), x = "NDVI") +
  theme(axis.text=element_text(size=11), 
        axis.title=element_text(size=12), 
        axis.line=element_line(size=0.5),
        axis.ticks=element_line(size=0.5), 
        strip.text=element_text(size=11), 
        strip.background=element_blank(), 
        panel.border=element_blank(), 
        panel.grid=element_blank(), legend.position = 'none') +
  scale_color_manual(values = rep('darkgray', 34))
ndvi
ggsave(file = "./Figures/AbundanceNDVI.jpg", plot = ndvi, width = 8, height = 6)


#########################################
# Figure Abundance

ab_cols <- grepl(colnames(mc), pattern = "^AByearsp\\[")
abpost <- mc[,ab_cols]/9.6 #transform into density by dividing total abundance across 8 transects by total area for the 8 transects (calculated from the polygons after buffering transects)

sppnames <-  c("Hartebeest", 'Eland','Giraffe',"Grant's gazelle",'Impala',"Thomson's gazelle",'Topi','Warthog','Wildebeest',"Plains zebra")
sppnames <- rep(sppnames, each = nreps)
season <- rep(Time, nspec)

forplot <- NULL
for (i in 1:length(sppnames)){
  temp <- data.frame(Species = sppnames[i], time = season[i], MedianAb = median(abpost[,i]), ci(abpost[,i], 0.89,'HDI'), mean = mean(abpost[,i]), sd = sd(abpost[,i]))
  forplot <- rbind(temp, forplot)
}


(Figure3 <- ggplot(forplot, aes(x=time, y=MedianAb)) + 
    facet_wrap(~Species, scales="free_y", nrow = 4, ncol = 3) +
    labs(y=expression("Estimated ind/km"^2), x = "Time") + 
    theme(axis.text=element_text(size=8), 
          axis.title=element_text(size=12), 
          axis.line=element_line(size=0.5),
          axis.ticks=element_line(size=0.5), 
          strip.text=element_text(size=11), 
          strip.background=element_blank(), 
          panel.border=element_blank(), 
          panel.grid=element_blank(),
          legend.position = c(0.85, 0.1), legend.title = element_blank()) +
    geom_point(size=2) + geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0, size=0.5, color = 'darkgrey') + scale_x_date(date_breaks = "6 months", date_labels = "%Y %b")) 


ggsave(file = "./Figures/AbundanceTrend.jpg", plot = Figure3, width = 11, height = 6)



### Calculate precision in abundance estimates

CV = data.frame(Species = forplot$Species, Date = forplot$time, cv = forplot$sd/forplot$mean*100)
CVmedians <- CV %>% group_by(Species) %>% dplyr::summarise(median(cv)) %>% as.data.frame()
median(CV$cv)

(CVboxplot <- ggplot(CV, aes(x=reorder(Species,-cv), y=cv)) + geom_boxplot() +  geom_jitter(alpha = 0.5) + labs(x = "", y = "CV (%)") + 
  theme(strip.text = element_blank(),strip.background=element_blank(), 
        legend.position = 'none'))

ggsave(file = "./Figures/BoxplotCVDens.jpg", plot = CVboxplot, width = 15, height = 6)

#Proportion of surveys with CV < 20, <30, <40
sum(CV$cv < 20)/300
sum(CV$cv < 30)/300
sum(CV$cv < 40)/300

##############################################
# Calculate wet season differences in density

density <- NULL
for (i in 1:length(sppnames)){
  temp <- data.frame(Species = sppnames[i], time = season[i], Density = abpost[,i])
  density <- rbind(temp, density)
}

density2018 <- subset(density, time %in% c("2018-01-01", "2018-02-01", "2018-03-01"))
density2019 <- subset(density, time %in% c("2019-01-01", "2019-02-01", "2019-03-01"))
density2020 <- subset(density, time %in% c("2020-01-01", "2020-02-01", "2020-03-01"))

summary2018 <- density2018 %>% group_by(Species) %>% dplyr::summarise(median(Density), ci(Density, 0.89)) %>% mutate(Year = 2018)
summary2019 <- density2019 %>% group_by(Species) %>% dplyr::summarise(median(Density), ci(Density, 0.89)) %>% mutate(Year = 2019)
summary2020 <- density2020 %>% group_by(Species) %>% dplyr::summarise(median(Density), ci(Density, 0.89)) %>% mutate(Year = 2020)

fordensplot <- rbind(summary2018,summary2019,summary2020)

(FigureDens <- ggplot(fordensplot, aes(x=Year, y=`median(Density)`)) + 
    facet_wrap(~Species, scales="free_y", nrow = 4, ncol = 3) +
    labs(y=expression("Estimated ind/km"^2), x = "Year") + 
    theme(axis.text=element_text(size=8), 
          axis.title=element_text(size=12), 
          axis.line=element_line(size=0.5),
          axis.ticks=element_line(size=0.5), 
          strip.text=element_text(size=11), 
          strip.background=element_blank(), 
          panel.border=element_blank(), 
          panel.grid=element_blank()) +
    geom_point(size=2) + geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0, size=0.5, color = 'darkgrey') + scale_x_discrete(limits=c(2018, 2019, 2020))) 
ggsave(file = "./Figures/WetSeasonDensity.jpg", plot = FigureDens, width = 8, height = 6)
