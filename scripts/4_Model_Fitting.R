#----------------------------------------------------------------------------------------------------------------------------------#
# Script by: Lucien Fitzpatrick
# Project: Bur Oak Forecasting (Quercus Quest)
# Purpose: This script downloads the NPN data and cleans them
# Inputs: 
# Outputs:
# Notes: 
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

#Creating a data frame to check the model output
Check <- data.frame()
l <- 1
#We are using bloom to determine the species becasue it is missing two species: fusiformis and laurifolia


#Creating the row for this species information
Check[l, "species"] <- "Quercus macrocarpa"
Check[l, "bur.nObs"] <- nrow(dat.burst)


SP.burst <- as.data.frame(table(dat.burst[, "state"]))
colnames(SP.burst) <- c("State", "Freq")

#Checking the mean
Check[l, "mean.bur"] <- mean(dat.burst$GDD5.cum, na.rm = T)
Check[l, "sd.bur"] <- sd(dat.burst$GDD5.cum, na.rm = T)

hierarchical_regression <- "
model{
  for(k in 1:n){
      mu[k] <- THRESH[st[k]]  
      y[k] ~ dnorm(mu[k], sPrec)
  }
  
  for(j in 1:nSt){
    THRESH[j] ~ dnorm(Tprior, aPrec)
  }
  aPrec ~ dgamma(1, 1)
  Tprior ~ dunif(1, 600)
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
  inits[[i]] <- list(.RNG.name = "base::Super-Duper",
                     .RNG.seed = as.integer(rngs[i])
  )
}

#---------------------------------------------------------#
#This section actually runs the model and then provides ways to check the output and clean it
#---------------------------------------------------------#
bur.model   <- jags.model (file = textConnection(hierarchical_regression),
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

names(bur.df)[grep("THRESH", names(bur.df))] <- paste(SP.burst$State)

write.csv(bur.df, file.path("../data_processed/model_output", paste0("Quercus_Quest_TT_model_bur.csv")), row.names=F)



#---------------------------------------------------#
#Here we start with the variable thermal time model. THis is like GDDX instead of GDD5 or GDD0
#---------------------------------------------------#
dat.burst <- read.csv("../data_processed/Full_Bur_Obs.csv")
dat.burst$pln <- as.numeric(factor(dat.burst$individual_id))
dat.burst$loc <- as.numeric(factor(dat.burst$site_id))
dat.burst$sp <- as.numeric(factor(dat.burst$species))
summary(dat.burst)

dat.burst <- dat.burst[dat.burst$state == "IL",]
dat.burst <- dat.burst[1:5,]

lat.calc <- read.csv("../data_processed/Daymet_clean_data.csv")

burst.list <- as.list(dat.burst)
for(i in 1:max(lengths(burst.list))){
  yr <- burst.list[["year"]][i]
  yday <- burst.list[["yday"]][i]
  site <- burst.list[["site_id"]][i]
  temp <- lat.calc[lat.calc$year == yr & lat.calc$site == site,]
  burst.list[["temp.vec"]][i] <- list(temp$TMEAN)
}

start <- c(2)
end <- c(365)

for(i in 2:nrow(dat.burst)){
  start[i] <- ((i-1)*365) + 2
  end[i] <- (i*365)
}

hierarchical_regression <- "
model{
  for(k in 1:n){
  
      cs[k,1] <- ifelse(temp[(k*365)-364] < 5, 0, (temp[(k*365)-364]-5))
      
      for (t in start[k]:end[k]) {
        new.tmp[k,t] <- ifelse(temp[t] < 5, 0, (temp[t]-5))

        cs[k,(t-((k-1)*365))] <- cs[k, ((t-((k-1)*365)) - 1)] + new.tmp[k,t]
      }
      
      mu[k] <- cs[k,] >= 200
      y[k] ~ dnorm(mu[k], sPrec)
  }
  
 #for(j in 1:nSt){
 #   THRESH[j] ~ dnorm(Tprior, aPrec)
 #}

  #aPrec ~ dgamma(.01, .01)
  #Tprior ~ dunif(1, 1000)
  sPrec ~ dgamma(.01, .01)

}
"

bur.list <- list(y = burst.list[["yday"]], n = length(burst.list[["yday"]]), 
                 temp = as.numeric(unlist(burst.list[["temp.vec"]])), start = as.numeric(start), end = as.numeric(end),
                 st = as.numeric(factor(burst.list[["state"]])), nSt = length(unique(burst.list[["state"]])))


#Setting the number of MCMC chains and their parameters
nchain = 3
inits <- list()
rngs <- c(737, 874, 869)
for(i in 1:nchain){
  inits[[i]] <- list(.RNG.name = "base::Super-Duper",
                     .RNG.seed = as.integer(rngs[i])
  )
}

#---------------------------------------------------------#
#This section actually runs the model and then provides ways to check the output and clean it
#---------------------------------------------------------#
bur.model   <- jags.model (file = textConnection(hierarchical_regression),
                           data = bur.list,
                           inits = inits,
                           n.chains = 3)


#Converting the ooutput into a workable format
bur.out   <- coda.samples (model = bur.model,
                           variable.names = c("cs"),
                           n.iter = 7000)


#Checking that convergence happened
bud.con <- gelman.diag(bur.out)
bud.con

#Removing burnin before convergence occurred
burnin = 6900                                ## determine convergence from GBR output
bur.burn <- window(bur.out,start=burnin)  ## remove burn-in
summary(bur.burn)

#converting them into a dataframe 
bur.df <- as.data.frame(as.matrix(bur.burn))

#calculating sd from the precison
bur.df$sd <- 1/sqrt(bur.df[,"aPrec"])

summary(bur.df)

bur.df$Species <- unique(dat.burst$species)

names(bur.df)[grep("THRESH", names(bur.df))] <- paste(SP.burst$State)

write.csv(bur.df, file.path("../data_processed/model_output", paste0("Quercus_Quest_VarTT_model_bur.csv")), row.names=F)





#---------------------------------------------------#
#Here we start with the GDD0 because it's the easiest to convert and Lucien was curious
#Does not work yet
#---------------------------------------------------#
hierarchical_regression <- "
model{
  for(k in 1:n){
      mu[k] <- THRESH[st[k]]  
      y[k] ~ dnorm(mu[k], sPrec)
  }
  
  for(j in 1:nSt){
    THRESH[j] ~ dnorm(Tprior, aPrec)
  }
  aPrec ~ dgamma(1, 1)
  Tprior ~ dunif(1, 600)
  sPrec ~ dgamma(.01, .01)

}
"

bur.list <- list(y = dat.burst$GDD0.cum, n = length(dat.burst$GDD0.cum),
                 st = as.numeric(factor(dat.burst$state)), nSt = length(unique(dat.burst$state)))


#Setting the number of MCMC chains and their parameters
nchain = 3
inits <- list()
rngs <- c(737, 874, 869)
for(i in 1:nchain){
  inits[[i]] <- list(.RNG.name = "base::Super-Duper",
                     .RNG.seed = as.integer(rngs[i])
  )
}

#---------------------------------------------------------#
#This section actually runs the model and then provides ways to check the output and clean it
#---------------------------------------------------------#
bur.model   <- jags.model (file = textConnection(hierarchical_regression),
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

names(bur.df)[grep("THRESH", names(bur.df))] <- paste(SP.burst$State)

write.csv(bur.df, file.path("../data_processed/model_output", paste0("Quercus_Quest_GDD0_model_bur.csv")), row.names=F)
