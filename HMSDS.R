###########################################################################
# Monitoring spatiotemporal dynamics of large herbivores across an African
# rangeland using hierarchical multi-species distance sampling

# Ramiro D. Crego, Harry B.M. Wells, Grant Connette, Jared A. Stabach,
# Naitareu Soit & Stewart Thompson

# R and Nimble model code 

# Set seed
set.seed(1984)

# Load Libraries
library(nimble)
library(coda)
library(mcmcOutput)

# Load Data
#Abundance
load("./Data/MultiSpDSdata")
#Covariates
load("./Data/NDVI")

attach(DSdata)

# BUGS Model
HMSDSModel <- nimbleCode({
  # PRIORS
  
  #Mu gamma0
  mu_g ~ dnorm(0,0.01)          #Mean
  sig_g ~ dunif(0,10)         #SD
  tau_g <- pow(sig_g,-2)    #Precision
  
  #gamma1
  gamma1 ~ dnorm(0, 0.01) 
  
  #Alpha0
  mu_a0 ~ dnorm(0, 0.01)        #Mean
  sig_a0 ~ dunif(0, 10)      #SD
  tau_a0 <- pow(sig_a0, -2)    #Precision
  
  #Alpha1
  mu_a1 ~ dnorm(0, 0.01)        #Mean
  sig_a1 ~ dunif(0, 10)      #SD
  tau_a1 <- pow(sig_a1, -2)    #Precision
  
  #Alpha2
  mu_a2 ~ dnorm(0, 0.01)        #Mean
  sig_a2 ~ dunif(0, 10)      #SD
  tau_a2 <- pow(sig_a2, -2)    #Precision
  
  #Alpha3
  mu_a3 ~ dnorm(0, 0.01)        #Mean
  sig_a3 ~ dunif(0, 10)      #SD
  tau_a3 <- pow(sig_a3, -2)    #Precision
  
  #Overdispersion
  r.N ~ dunif(0,50)            #Number of groups
  r.G ~ dunif(0,50)            #Group size
  
  for(s in 1:nspec){
    
    #Gammas0
    gamma0[s] ~ dnorm(mu_g, tau_g)  #Intercept parameter
    
    
    #Expected Number of Groups
    #Random factor for time
    mu_alpha0[s] ~ dnorm(mu_a0, tau_a0)    #Intercept parameter
    sig_alpha0[s] ~ dunif(0, 10)      #SD
    tau_alpha0[s] <- pow(sig_alpha0[s], -2)    #Precision
    
    for (t in 1:nreps) {
      alpha0[t,s] ~ dnorm(mu_alpha0[s], tau_alpha0[s])
    }
    
    # Covariables
    alpha1[s] ~ dnorm(mu_a1, tau_a1)    #Effect parameter NDVI
    alpha2[s] ~ dnorm(mu_a2, tau_a2)    #Effect parameter NDVI2
    alpha3[s] ~ dnorm(mu_a3, tau_a3) # trend
    
    #Average Group Size
    #Random factor for time
    mu_beta0[s] ~ dnorm(0, 0.01)    #Intercept parameter
    sig_beta0[s] ~ dunif(0, 10)      #SD
    tau_beta0[s] <- pow(sig_beta0[s], -2)    #Precision
    
    for (t in 1:nreps) {
      beta0[t,s] ~ dnorm(mu_beta0[s], tau_beta0[s])
    }
    
    #Scale parameter
    sigma[s] <- exp(gamma0[s] + gamma1 * BMass[s])
    
    for(j in 1:nsites){
      # LIKELIHOOD
      
      for(t in 1:nreps){
        
        #Construct cell probabilities for nG cells using numerical integration
        #Sum of the area (rectangles) under the detection function
        
        for(k in 1:nD){
          
          #Half-normal detection function at midpt (length of rectangle)
          g[k,t,j,s] <- exp(-mdpt[k]*mdpt[k]/(2*sigma[s]*sigma[s]))
          
          #Proportion of each interval (width of rectangle) for both sides of the transect
          pi[k,t,j,s] <- v/B
          
          #Detection probability for each distance class k (area of each rectangle)
          f[k,t,j,s] <- g[k,t,j,s] * pi[k,t,j,s]
          
          #Conditional detection probability (scale to 1)
          fc[k,t,j,s] <- f[k,t,j,s]/pcap[t,j,s]
          
        }#end k loop
        
        #Detection probability at each transect (sum of rectangles)
        pcap[t,j,s] <- sum(f[1:nD,t,j,s])
        
        #Observed population @ each t,j,s (N-mixture)
        y[t,j,s] ~ dbin(pcap[t,j,s], N[t,j,s])
        
        #Latent Number of Groups @ each t,j,s (negative binomial)
        N[t,j,s] ~ dpois(lambda.star[t,j,s])
        
        #Expected Number of Groups
        lambda.star[t,j,s] <- rho[t,j,s] * lambda[t,j,s]
        
        #Overdispersion parameter for Expected Number of Groups
        rho[t,j,s] ~ dgamma(r.N, r.N)
        
        #Linear predictor for Expected Number of Groups
        lambda[t,j,s] <- exp(alpha0[t,s] + alpha1[s] * NDVI[t,j] + alpha2[s] * NDVI2[t,j] + alpha3[s] * (t - 15.5)) 

      }#end t loop
      
      #Mean detection probability @ each j,s
      #psite[j,s] <- mean(pcap[1:nreps, j, s])
      
    }#end j loop
    
  }#end s loop
  
  for(s in 1:nspec){
      for(t in 1:nreps){
        
       #Expected Group Size
        gs.lam.star[t,s] <- gs.lam[t,s] * gs.rho[t,s]
        
        #Overdispersion parameter for Expected Group Size
        gs.rho[t,s] ~ dgamma(r.G, r.G)
        
        #Linear predictor for Expected Group Size
        gs.lam[t,s] <- exp(beta0[t,s])

      }#end t loop
  }#end s loop
  
  for(i in 1:nobs){
    #Observed distance classes
    dclass[i] ~ dcat(fc[1:nD, rep[i], site[i], spec[i]])
    #Observed Group Size (zero truncated negative binomial)
    gs[i] ~ T(dpois(gs.lam.star[rep[i], spec[i]]),1,)
  }#end i loop

  for(s in 1:nspec){
    for(j in 1:nsites){
      for(t in 1:nreps){
        #Estimated number of individuals per transect
        ABrep[t,j,s] <- lambda.star[t,j,s] * gs.lam.star[t,s]
      }#end t loop
    }#end j loop
  }#end s loop
  
  for(s in 1:nspec){
    for(t in 1:nreps){
      # Total estimated number of individuals per month
      AByearsp[t,s] <- sum(ABrep[t,1:nsites,s])
    }#end t loop
  }#end s loop
  
})

#Compile data and constants

data <- list(y = y, dclass = dclass,  NDVI = NDVI, NDVI2 = pow(NDVI,2), BMass = as.vector(BMsc), gs = gs) 

const <- list(spec = spec, 
              site = site, 
              rep = rep, 
              nspec = nspec, 
              nD = nD, 
              v = v, 
              B = B, 
              mdpt = mdpt,
              nobs = nobs,
              nsites = nsites, 
              nreps = nreps
              )

# Inital values

Nst <- y + 1

inits <- list(mu_g = runif(1, 3, 6), 
              sig_g=runif(1, 0, 0.5),
              gamma1 = rnorm(1),
              mu_a0 = runif(1, -4, -2), 
              sig_a0 = runif(1, 0, 1),
              mu_a1 = rnorm(1), 
              sig_a1 = runif(1, 0, 1),
              mu_a2 = rnorm(1), 
              sig_a2 = runif(1, 0, 1),
              mu_a3 = rnorm(1), 
              sig_a3 = runif(1, 0, 1),
              mu_beta0 = runif(nspec, 3, 4), 
              sig_beta0 = runif(nspec, 1.5 ,2.5),
              r.N = runif(1, 0, 5), 
              r.G = runif(1, 0, 5), 
              N = Nst
)

# Parameters to save

params <- c('mu_g', 
             'mu_a0','mu_a1','mu_a2', 'mu_a3',
             'gamma0','gamma1',
             'mu_alpha0','alpha0','alpha1','alpha2','alpha3',
             'mu_beta0','beta0',
             'sig_g', 
             'sig_a0','sig_a1','sig_a2','sig_a3','sig_alpha0',
             'sig_beta0','AByearsp',"N", "lambda.star", "gs.lam.star","pcap")

# MCMC settings

modelnim <- nimbleModel(code = HMSDSModel, name = 'MultSpDist', constants = const, data = data, inits = inits)

modolconf <- configureMCMC(modelnim, monitors = params) # for parameters
modMCMC <- buildMCMC(modolconf)
CdistModel <- compileNimble(modelnim)
CdistMCMC <- compileNimble(modMCMC, project = modelnim)
runMCMC(CdistMCMC, niter = 2) #Check for a few iterations
runMCMC_samples <- runMCMC(CdistMCMC, 
                           nchains = 3,
                           nburnin = 100000,  
                           niter = 200000, 
                           thin = 100,
                           samples = TRUE,
                           samplesAsCodaMCMC = T,
                           summary = F,
                           WAIC = FALSE,
                           perChainWAIC = FALSE)

mc <- mcmcOutput(runMCMC_samples)
#save(mc, file = "./ModelOutput/HMSDS_MC.Rdata")

summary <- summary(mc, n.eff=TRUE)
summary <- summary[c(2701:3341,6042:6066,8467:8491),]
View(summary)


# Appendix: Table S1
library(bayestestR)

parnames <- rownames(summary)
table1 <- data.frame()
for (i in c(2701:3341,6042:6066,8467:8491)){
  t <- data.frame(Mean = mean(mc[,i]), Median = median(mc[,i]), ci(mc[,i], 0.89))
  table1 <- rbind(table1, t)
}
table1$Gelman <- summary$Rhat
table1$n.eff <- summary$n.eff
table1$Param <- parnames
table1
#write.csv(table1, file = "./ModelOutput/AppendixTable1.csv")


#### Posterior checks
library(MCMCvis)

diagPlot(mc,c('mu_g', 
              'mu_a0','mu_a1','mu_a2', 'mu_a3',
              'gamma0','gamma1',
              'mu_alpha0','alpha0','alpha1','alpha2','alpha3',
              'mu_beta0', 'beta0',
              'sig_g', 
              'sig_a0','sig_a1','sig_a2','sig_a3',
              'sig_alpha0','sig_beta0'))

#### Check for prior influence on posterior estimates

# Simulate data from the prior used in the model
# Number of iterations should equal the number of draws times the number of chains 
niter <- nrow(mc)
PR <- rnorm(niter, 0, 10) 

MCMCtrace(mc, params = 'mu_g', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'mu_a0', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'mu_a1', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'mu_a2', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'mu_a3', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'gamma0', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'gamma1', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'mu_alpha0', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'alpha0', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'alpha1', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'alpha2', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'alpha3', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'mu_beta0', priors = PR, pdf = FALSE)
MCMCtrace(mc, params = 'beta0', priors = PR, pdf = FALSE)

PR2 <- runif(niter, 0, 10)

MCMCtrace(mc, params = 'sig_g', priors = PR2, pdf = FALSE)
MCMCtrace(mc, params = 'sig_a0', priors = PR2, pdf = FALSE)
MCMCtrace(mc, params = 'sig_a1', priors = PR2, pdf = FALSE)
MCMCtrace(mc, params = 'sig_a2', priors = PR2, pdf = FALSE)
MCMCtrace(mc, params = 'sig_a3', priors = PR2, pdf = FALSE)
MCMCtrace(mc, params = 'sig_alpha0', priors = PR2, pdf = FALSE)
MCMCtrace(mc, params = 'sig_beta0', priors = PR2, pdf = FALSE)
# Following Gimenez et al. (2009), overlap should be < 0.35


