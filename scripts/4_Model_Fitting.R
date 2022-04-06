#----------------------------------------------------------------------------------------------------------------------------------#
# Script by: Lucien Fitzpatrick
# Project: Bur Oak Forecasting (Quercus Quest)
# Purpose: This script takes our filtered NPN data and fits our different model paramters.
# Inputs: Full_Bur_Obs.csv which contains all npn observations nearby source populations matched with weather metrics
#         Daymet_clean_data.csv which contains yearly weather data for each site of observation
# Outputs: Parameter distributions for models
# Notes: Current models are YDAY null model, GDD5, PTTGDD5, linear GTmean
#-----------------------------------------------------------------------------------------------------------------------------------#
library(rjags)
library(coda)
library(dplyr)
library(tidyr)

dat.burst <- read.csv("../data_processed/Full_Bur_Obs.csv")
dat.burst$pln <- as.numeric(factor(dat.burst$individual_id))
dat.burst$loc <- as.numeric(factor(dat.burst$site_id))
dat.burst$sp <- as.numeric(factor(dat.burst$species))
summary(dat.burst)

lat.calc <- read.csv("../data_processed/Daymet_clean_data.csv")

#---------------------------------------------------#
#Here we start with the GDD5 thermal time model.
#This model is fit using calculated GDD5.cum values (y = GDD5.cum on date of known observation)
#This model fits an accumulation threshold (THRESH) for each source location
#---------------------------------------------------#

YDAY_model <- "
model{
  #This is the data model where the y=YDAY is used to fit THRESH
  for(k in 1:n){
      mu[k] <- THRESH[st[k]]  
      y[k] ~ dnorm(mu[k], sPrec)
  }
  
  #This is the process model and fits the accumulation threshold for each site
  for(j in 1:nSt){
    THRESH[j] ~ dnorm(Tprior, aPrec)T(0, )
  }
  
  #This is where our priors are defined
  aPrec ~ dgamma(.01, .01)
  Tprior ~ dunif(90,120)
  sPrec ~ dgamma(.01, .01)

}
"

bur.list <- list(y = dat.burst$yday, n = length(dat.burst$yday),
                 st = as.numeric(factor(dat.burst$state)), nSt = length(unique(dat.burst$state)))


#Setting the number of MCMC chains and their parameters
nchain = 3
inits <- list()
rngs <- c(737, 874, 869)
for(i in 1:nchain){
  inits[[i]] <- list(aPrec = runif(1,1/5000,1/100), THRESH = runif(3,100,300), sPrec = runif(1,1/100,1/20),
                     .RNG.name = "base::Super-Duper",
                     .RNG.seed = as.integer(rngs[i])
  )
}

#This section actually runs the model and then provides ways to check the output and clean it
bur.model   <- jags.model (file = textConnection(YDAY_model),
                           data = bur.list,
                           inits = inits,
                           n.chains = 3)


#Converting the ooutput into a workable format
bur.out   <- coda.samples (model = bur.model,
                           variable.names = c("THRESH", "aPrec"),
                           n.iter = 700000)


#Checking that convergence happened
bud.con <- gelman.diag(bur.out)
bud.con

#Removing burnin before convergence occurred
burnin = 690000                                ## determine convergence from GBR output
bur.burn <- window(bur.out,start=burnin)  ## remove burn-in
summary(bur.burn)

#converting them into a dataframe 
bur.df <- as.data.frame(as.matrix(bur.burn))

#calculating sd from the precison
bur.df$sd <- 1/sqrt(bur.df[,"aPrec"])

summary(bur.df)

bur.df$Species <- unique(dat.burst$species)

names(bur.df)[grep("THRESH", names(bur.df))] <- unique(dat.burst$state)

write.csv(bur.df, file.path("../data_processed/model_output", paste0("Quercus_Quest_YDAY_model_bur.csv")), row.names=F)


#---------------------------------------------------#
#Here we start with the GDD5 thermal time model.
#This model is fit using calculated GDD5.cum values (y = GDD5.cum on date of known observation)
#This model fits an accumulation threshold (THRESH) for each source location
#---------------------------------------------------#

GDD5_model <- "
model{
  #This is the data model where the y=GDD5.cum is used to fit THRESH
  for(k in 1:n){
      mu[k] <- THRESH[st[k]]  
      y[k] ~ dnorm(mu[k], sPrec)
  }
  
  #This is the process model and fits the accumulation threshold for each site
  for(j in 1:nSt){
    THRESH[j] ~ dnorm(Tprior, aPrec)T(0, )
  }
  
  #This is where our priors are defined
  aPrec ~ dgamma(.01, .01)
  Tprior ~ dnorm(180, 50)
  sPrec ~ dgamma(.01, .01)

}
"

bur.list <- list(y = dat.burst$GDD5.cum, n = length(dat.burst$GDD5.cum),
                 st = as.numeric(factor(dat.burst$state)), nSt = length(unique(dat.burst$state)))


#Setting the number of MCMC chains and their parameters
nchain = 3
inits <- list()
rngs <- c(737, 874, 869)
for(i in 1:nchain){
  inits[[i]] <- list(aPrec = runif(1,1/5000,1/100), THRESH = runif(3,100,300), sPrec = runif(1,1/100,1/20),
                     .RNG.name = "base::Super-Duper",
                     .RNG.seed = as.integer(rngs[i])
  )
}

#This section actually runs the model and then provides ways to check the output and clean it
bur.model   <- jags.model (file = textConnection(GDD5_model),
                           data = bur.list,
                           inits = inits,
                           n.chains = 3)


#Converting the ooutput into a workable format
bur.out   <- coda.samples (model = bur.model,
                           variable.names = c("THRESH", "aPrec"),
                           n.iter = 700000)


#Checking that convergence happened
bud.con <- gelman.diag(bur.out)
bud.con

#Removing burnin before convergence occurred
burnin = 690000                                ## determine convergence from GBR output
bur.burn <- window(bur.out,start=burnin)  ## remove burn-in
summary(bur.burn)

#converting them into a dataframe 
bur.df <- as.data.frame(as.matrix(bur.burn))

#calculating sd from the precison
bur.df$sd <- 1/sqrt(bur.df[,"aPrec"])

summary(bur.df)

bur.df$Species <- unique(dat.burst$species)

names(bur.df)[grep("THRESH", names(bur.df))] <- unique(dat.burst$state)

write.csv(bur.df, file.path("../data_processed/model_output", paste0("Quercus_Quest_GDD5_model_bur.csv")), row.names=F)


#---------------------------------------------------#
#Next we have the PTTGDD thermal time model.
#This model is fit using calculated PTTGDD5.cum values (y = PTTGDD5.cum on date of known observation)
#PTTGDD is GDD5 * (day length in hours/24)
#This model fits an accumulation threshold (THRESH) for each source location
#---------------------------------------------------#
PTTGDD_model  <- "
model{
 #This is the data model where the y=PTTGDD5.cum is used to fit THRESH
  for(k in 1:n){
      mu[k] <- THRESH[st[k]]  
      y[k] ~ dnorm(mu[k], sPrec)
  }
  
  #This is the process model and fits the accumulation threshold for each site
  for(j in 1:nSt){
    THRESH[j] ~ dnorm(Tprior, aPrec)T(0, )
  }
  
  #This is where our priors are defined
  aPrec ~ dgamma(.01, .01)
  Tprior ~ dnorm(80, 20)
  sPrec ~ dgamma(.01, .01)

}
"

bur.list <- list(y = dat.burst$PTTGDD.cum, n = length(dat.burst$PTTGDD.cum),
                 st = as.numeric(factor(dat.burst$state)), nSt = length(unique(dat.burst$state)))


#Setting the number of MCMC chains and their parameters
nchain = 3
inits <- list()
rngs <- c(737, 874, 869)
for(i in 1:nchain){
  inits[[i]] <- list(aPrec = runif(1,1/5000,1/100), THRESH = runif(3,10,150), sPrec = runif(1,1/100,1/20),
                     .RNG.name = "base::Super-Duper",
                     .RNG.seed = as.integer(rngs[i])
  )
}

#This section actually runs the model and then provides ways to check the output and clean it
bur.model   <- jags.model (file = textConnection(PTTGDD_model),
                           data = bur.list,
                           inits = inits,
                           n.chains = 3)


#Converting the ooutput into a workable format
bur.out   <- coda.samples (model = bur.model,
                           variable.names = c("THRESH", "aPrec"),
                           n.iter = 700000)


#Checking that convergence happened
bud.con <- gelman.diag(bur.out)
bud.con

#Removing burnin before convergence occurred
burnin = 690000                                ## determine convergence from GBR output
bur.burn <- window(bur.out,start=burnin)  ## remove burn-in
summary(bur.burn)

#converting them into a dataframe 
bur.df <- as.data.frame(as.matrix(bur.burn))

#calculating sd from the precison
bur.df$sd <- 1/sqrt(bur.df[,"aPrec"])

summary(bur.df)

bur.df$Species <- unique(dat.burst$species)

names(bur.df)[grep("THRESH", names(bur.df))] <- unique(dat.burst$state)

write.csv(bur.df, file.path("../data_processed/model_output", paste0("Quercus_Quest_PTTGDD_model_bur.csv")), row.names=F)


#---------------------------------------------------#
#Here we start with the GDD5 thermal time model.
#This model is fit using calculated GDD5.cum values (y = GDD5.cum on date of known observation)
#This model fits an accumulation threshold (THRESH) for each source location
#---------------------------------------------------#

GTmean_model <- "
model{
  #This is the data model where the y=Yday is used to fit an intercept and slope
  for(k in 1:n){
      mu[k] <- int[st[k]] + slope[st[k]]*GTmean[k]  
      y[k] ~ dnorm(mu[k], sPrec)
  }
  
  #This is the process model and fits the accumulation threshold for each site
  for(j in 1:nSt){
    int[j] ~ dnorm(Tprior,iPrec)
    slope[j] ~ dnorm(0,xPrec)
  }
  
  #This is where our priors are defined
  Tprior ~ dunif(90,121)
  sPrec ~ dgamma(.01, .01)
  iPrec ~ dgamma(.01, .01)
  xPrec ~ dgamma(.01, .01)
  
}
"

bur.list <- list(y = dat.burst$yday, n = length(dat.burst$yday), 
                 st = as.numeric(factor(dat.burst$state)), nSt = length(unique(dat.burst$state)),
                 GTmean = dat.burst$GTmean)


#Setting the number of MCMC chains and their parameters
nchain = 3
inits <- list()
rngs <- c(737, 874, 869)
for(i in 1:nchain){
  inits[[i]] <- list(sPrec = runif(1,1/100,1/20),
                     .RNG.name = "base::Super-Duper",
                     .RNG.seed = as.integer(rngs[i])
  )
}

#This section actually runs the model and then provides ways to check the output and clean it
bur.model   <- jags.model (file = textConnection(GTmean_model),
                           data = bur.list,
                           inits = inits,
                           n.chains = 3)


#Converting the ooutput into a workable format
bur.out   <- coda.samples (model = bur.model,
                           variable.names = c("int", "slope", "iPrec", "xPrec"),
                           n.iter = 700000)


#Checking that convergence happened
bud.con <- gelman.diag(bur.out)
bud.con

#Removing burnin before convergence occurred
burnin = 690000                                ## determine convergence from GBR output
bur.burn <- window(bur.out,start=burnin)  ## remove burn-in
summary(bur.burn)

#converting them into a dataframe 
bur.df <- as.data.frame(as.matrix(bur.burn))

names(bur.df)[grep("int", names(bur.df))] <- paste0("int", unique(dat.burst$state))
names(bur.df)[grep("slope", names(bur.df))] <- paste0("slope", unique(dat.burst$state))

#calculating sd from the precison
bur.df$int.sd <- 1/sqrt(bur.df[,"iPrec"])
bur.df$slope.sd <- 1/sqrt(bur.df[,"xPrec"])

summary(bur.df)

bur.df$Species <- unique(dat.burst$species)

write.csv(bur.df, file.path("../data_processed/model_output", paste0("Quercus_Quest_GTmean_model_bur.csv")), row.names=F)

