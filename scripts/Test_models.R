
#---------------------------------------------------#
#Next we have the CDD0 thermal time model.
#This model is fit using calculated CDD0.cum values (y = CDD0.cum on date of known observation)
#This model fits an accumulation threshold (THRESH) for each source location
#Currently doesn't work since IL won't reach the temp for MN and OK predicts any sign of chilling will cause budburst
#This failure seems more an issue of the model being poor for our species. May doing a cdd2 would help?
#---------------------------------------------------#
CDD0_model  <- "
model{
 #This is the data model where the y=CDD0.cum is used to fit THRESH
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

bur.list <- list(y = dat.burst$CDD0.cum, n = length(dat.burst$CDD0.cum),
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
bur.model   <- jags.model (file = textConnection(CDD0_model),
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

write.csv(bur.df, file.path("../data_processed/model_output", paste0("Quercus_Quest_CDD0_model_bur.csv")), row.names=F)


#---------------------------------------------------#
#Next we have our GDDX thermal time model. This is curerntly our most complext
#This model is fit using the date of observed earliest budburst (y = day of year of known budburst observation)
#This model fits an accumulation threshold (THRESH) and the base temperature used for accumulation (T_base) for each source location
#CURRENTLY DOESN'T WORK
#---------------------------------------------------#
#Data is converted into a list so we can match the sites yearly temperature values to each observation

dat.burst <- read.csv("../data_processed/Test_Full_Bur_Obs.csv")
dat.burst$pln <- as.numeric(factor(dat.burst$individual_id))
dat.burst$loc <- as.numeric(factor(dat.burst$site_id))
dat.burst$sp <- as.numeric(factor(dat.burst$species))
summary(dat.burst)

lat.calc <- read.csv("../data_processed/Test_Daymet_clean_data.csv")

dat.burst <- dat.burst[dat.burst$state == "MA",]

burst.list <- as.list(dat.burst)
for(i in 1:max(lengths(burst.list))){
  yr <- burst.list[["year"]][i]
  yday <- burst.list[["yday"]][i]
  site <- burst.list[["site_id"]][i]
  temp <- lat.calc[lat.calc$year == yr & lat.calc$site == site,]
  burst.list[["temp.vec"]][i] <- list(temp$TMEAN)
}

rm(lat.calc)

#Creating out start and end terms used to make sure the model loops through each year/site individually
start <- c(2)
end <- c(365)

for(i in 2:nrow(dat.burst)){
  start[i] <- ((i-1)*365) + 2
  end[i] <- (i*365)
}

GDDX_model <- "
model{
  
  #This is the data model where y = yday is used to fir out THRESH and T_Base
  for(k in 1:n){
      
      #Below we calculate the cumulative sum of thermal time for each observation
      #This is tricky because JAGS is not procedural which means objects can only be defined once and everything is defined at the same time
      
      cs[k,1] <- ifelse(temp[1] < T_base, 0, (temp[1]-T_base))
      #cs[k,1] <- ifelse(temp[1] < 5, 0, (temp[1]-5))
      
      for (t in start[k]:end[k]) {
        new.tmp[k,t] <- ifelse(temp[t] < T_base, 0, (temp[t]-T_base))
        #new.tmp[k,t] <- ifelse(temp[t] < 5, 0, (temp[t]-5))

        cs[k,(t-((k-1)*365))] <- (cs[k, ((t-((k-1)*365)) - 1)] + new.tmp[k,t])
      }
      
      #In JAGS a comparision operator returns a 1 for true and a 0 for false
      #We use this to calculate the first yday the Threshold is crossed using some math
      
      check[k, 1:365] <- cs[k,] >= THRESH
      mu[k] <- (366-sum(check[k,]))
      y[k] ~ dnorm(mu[k], sPrec)
      
  }
  
 #This is the process model where the accumlation threshold (THRESH) and base temperature of accumulation (T_base) are fit for each site  

  THRESH ~ dnorm(Tprior, aPrec)T(0, )
  T_base ~ dnorm(TmpPrior, tPrec)

 
  #This is where our priors are defined
  Tprior ~ dnorm(180, 50)
  aPrec ~ dgamma(.01, .01)
  
  TmpPrior ~dunif(-10, 10)
  tPrec ~dgamma(.01, .01)
  
  sPrec ~ dgamma(.01, .01)

}
"

bur.list <- list(y = burst.list[["yday"]], n = length(burst.list[["yday"]]), 
                 temp = as.numeric(unlist(burst.list[["temp.vec"]])), start = as.numeric(start), end = as.numeric(end))


#Setting the number of MCMC chains and their parameters
nchain = 3
inits <- list()
rngs <- c(737, 874, 869)
for(i in 1:nchain){
  inits[[i]] <- list(THRESH = runif(1,100,300), T_base = runif(1, 0, 5),
                     .RNG.name = "base::Super-Duper",
                     .RNG.seed = as.integer(rngs[i])
  )
}

#This section actually runs the model and then provides ways to check the output and clean it
bur.model   <- jags.model (file = textConnection(GDDX_model),
                           data = bur.list,
                           inits = inits,
                           n.chains = 3)


#Converting the ooutput into a workable format
bur.out   <- coda.samples (model = bur.model,
                           variable.names = c("THRESH", "aPrec", "T_base", "tPrec"),
                           n.iter = 70000)


#Checking that convergence happened
bud.con <- gelman.diag(bur.out)
bud.con

#Removing burnin before convergence occurred
burnin = 69000                                ## determine convergence from GBR output
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

