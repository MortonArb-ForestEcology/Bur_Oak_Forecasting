#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick
# Project: Phenology forecasting app
# Purpose: This script is for the creation of the sql database that is pulled for the app
# Inputs: Thermal Time model output for Oak species
#         Quercuscollection googlesheet
#         Oak_collection_budburst.csv
#         Weather_ArbCOOP_historical_latest.csv created by M1_Meterology_download.R
#         Weather_Arb_forecast_ensemble_latest.csv created by M1_Meterology_download.R
# Outputs: Arb_Pheno.db created
# Notes: #Currently this establishes the database on a local device and then loads it with our information of interest
#-----------------------------------------------------------------------------------------------------------------------------------#
library(dplyr)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(tidyr)
path.weath <- "data_raw/meteorology/"
path.ghcn=c("data_raw/meteorology/GHCN_extracted/")
dir.met <- "data_raw/meteorology"

path.temp <- "MortonArb_Buroak_Forecast/data/"
if(!dir.exists(path.temp)) dir.create(path.temp, recursive=T, showWarnings = F)

path.met <- "MortonArb_Buroak_Forecast/data/meteorology"
if(!dir.exists(path.met)) dir.create(path.met, recursive=T, showWarnings = F)

path.burst <- file.path("MortonArb_Buroak_Forecast/data/budburst")
if(!dir.exists(path.burst)) dir.create(path.burst, recursive=T, showWarnings = F)

set.seed(903)
#Reading in our different models and formatting them so they can share a dataframe
path.mod <- "../../data_processed/model_output/"

GDD5.mod <- read.csv(paste0(path.mod, "Quercus_Quest_GDD5_model_bur.csv"))

GDD5.model <- gather(GDD5.mod, state, THRESH, IL:OK, factor_key=TRUE)
GDD5.model <- GDD5.model[c("Species","state","THRESH")]
GDD5.model$model <- "GDD5"
GDD5.model$int <- NA
GDD5.model$slope <- NA

GDD5.model <- do.call(rbind,lapply(split(GDD5.model, GDD5.model$state),function(x) x[sample(nrow(x), 1000), ]))

PTTGDD5.mod <- read.csv(paste0(path.mod, "Quercus_Quest_PTTGDD_model_bur.csv"))

PTTGDD5.model <- gather(PTTGDD5.mod, state, THRESH, IL:OK, factor_key=TRUE)
PTTGDD5.model <- PTTGDD5.model[c("Species","state","THRESH")]
PTTGDD5.model$model <- "PTTGDD5"
PTTGDD5.model$int <- NA
PTTGDD5.model$slope <- NA

PTTGDD5.model <- do.call(rbind,lapply(split(PTTGDD5.model, PTTGDD5.model$state),function(x) x[sample(nrow(x), 1000), ]))


YDAY.mod <- read.csv(paste0(path.mod, "Quercus_Quest_YDAY_model_bur.csv"))

YDAY.model <- gather(YDAY.mod, state, THRESH, IL:OK, factor_key=TRUE)
YDAY.model <- YDAY.model[c("Species","state","THRESH")]
YDAY.model$model <- "YDAY"
YDAY.model$int <- NA
YDAY.model$slope <- NA

YDAY.model <- do.call(rbind,lapply(split(YDAY.model, YDAY.model$state),function(x) x[sample(nrow(x), 1000), ]))


GTmean.mod <- read.csv(paste0(path.mod, "Quercus_Quest_GTmean_model_bur.csv"))

GTmean.model <- gather(GTmean.mod, intstate, int, intIL:intOK, factor_key=TRUE)
GTmean.model <- gather(GTmean.model, slopestate, slope, slopeIL:slopeOK, factor_key=TRUE)

GTmean.model$state <- car::recode(GTmean.model$intstate, "'intIL'='IL'; 'intMN'='MN'; 'intOK'='OK'")

GTmean.model <- GTmean.model[c("Species","state", "int", "slope")]
GTmean.model$model <- "GTmean"
GTmean.model$THRESH <- NA

GTmean.model <- do.call(rbind,lapply(split(GTmean.model, GTmean.model$state),function(x) x[sample(nrow(x), 1000), ]))


#Combining all the models into one dataframe despite different parameters
comb.model <- rbind(GDD5.model, PTTGDD5.model, GTmean.model, YDAY.model)

comb.model$state <- factor(comb.model$state, levels = c("MN", "IL", "OK"))
comb.model$model <- factor(comb.model$model, levels = c("YDAY", "GDD5", "GTmean", "PTTGDD5"))
rownames(comb.model) <- NULL
write.csv(comb.model, file.path(path.burst,"Model_params_distribution.csv"),row.names = F)


#rownames(b.model) <- NULL

#Reading in the oak observations
dat.b <- read.csv(file.path(path.mod, "../Full_Bur_Obs.csv"))

#Reading in the historical weather
dat.ghcn <- read.csv(file.path(dir.met, "Weather_ArbCOOP_historical_latest.csv"))

dat.ghcn$DATE <- as.Date(dat.ghcn$DATE)

#Reading in our latest forecast
dat.forecast <- read.csv(file.path(dir.met, paste0("Mortonarb_daily_FORECAST-READY-LONGRANGE_", Sys.Date(),".csv")))
vars.agg <- c("TMEAN", "GDD0.cum", "GDD5.cum", "CDD0.cum", "CDD2.cum", "PTTGDD5.cum")
ens.forecast <- list()

ens.good <- unique(dat.forecast[!is.na(dat.forecast$TMAX) & dat.forecast$DATE == "2022-04-20", "ENS"])

dat.forecast <- dat.forecast[dat.forecast$ENS %in% ens.good,]

#Populating our uncertainty
for(VAR in vars.agg){
  ens.forecast[[VAR]] <- aggregate(dat.forecast[,VAR],
                                   by=dat.forecast[,c("DATE", "YDAY", "TYPE")],
                                   FUN=mean, na.rm=T)
  names(ens.forecast[[VAR]])[names(ens.forecast[[VAR]])=="x"] <- "mean"
  ens.forecast[[VAR]]$min <- aggregate(dat.forecast[,VAR],
                                       by=dat.forecast[,c("DATE", "YDAY", "TYPE")],
                                       FUN=min, na.rm=T)$x
  ens.forecast[[VAR]]$max <- aggregate(dat.forecast[,VAR],
                                       by=dat.forecast[,c("DATE", "YDAY", "TYPE")],
                                       FUN=max, na.rm=T)$x
}

if(Sys.Date()<=as.Date(paste0(lubridate::year(Sys.Date()), "-06-20"))){
  dat.ghcn$threshB <- dat.ghcn$GDD5.cum
  ens.forecast$threshB <- ens.forecast$GDD5.cum
} else {
  dat.ghcn$threshB <- dat.ghcn$CDD2.cum 
  ens.forecast$threshB <- ens.forecast$CDD0.cum
}

fc.sp <- as.data.frame(do.call(rbind, ens.forecast))
fc.sp$VAR <- gsub("\\..*","",(row.names(fc.sp)))
rownames(fc.sp) <- NULL
colnames(fc.sp) <- c("PRED.DATE", "YDAY", "TYPE", "mean", "min", "max", "VAR")

write.csv(fc.sp, file.path(path.temp, "meteorology", paste0("Forecast_data_", Sys.Date(),".csv")), row.names = F)

calc.thresh.bud <- function(dat, VAR, THRESH){
  min(dat[which(dat[,VAR] >= THRESH),"YDAY"])
}

for(MOD in unique(comb.model$model)){
  #subsetting to the particular model
  b.model <- comb.model[comb.model$model == MOD,]
  #Need to optomize this later
  pred.st <- data.frame()
  
  lim.comb <- data.frame(matrix(ncol = 6, nrow = length(unique(b.model$state))))
  colnames(lim.comb) <- c("state", "lb", "ub", "mean", "min", "max")
  
  count <- 1
  
  for(ST in unique(b.model$state)){
    
    #GTmean is a linear model so it needs to be handled differently
    if(MOD == "GTmean"){
      int <- b.model[b.model$state == ST, "int"]
      slope <- b.model[b.model$state == ST, "slope"]
      
      lin <- data.frame(int,slope)
      
      pred.array <- array(dim=c(length(int), length(unique(dat.forecast$ENS))))
      
      #Here we predict using the intercept and the slope to predict the yday
      ens <- unique(dat.forecast$ENS)
      for(i in 1:length(ens)){
        GT <- unique(dat.forecast[dat.forecast$ENS==ens[i],"GTmean"])[1]
        pred.array[,i] <- lin$int + (lin$slope* GT) 
      }
    
    }else if(MOD == "YDAY"){ #Now our basic yday model which uses a basic yday prediction. This is our null model
      thresh <- b.model[b.model$state == ST, "THRESH"]
      
      pred.array <- array(dim=c(length(thresh), length(unique(dat.forecast$ENS))))
      #We want to pull from the Budburst_Model table for out GDD predictions
      ens <- unique(dat.forecast$ENS)
      for(i in 1:length(ens)){
        pred.array[,i] <- thresh
      }
    } else { #Here we calculate GDD5 and PTTGDD5 using our calc.thresh.bud function
      Var <- paste0(MOD, ".cum")
      
      thresh <- b.model[b.model$state == ST, "THRESH"]
      
      pred.array <- array(dim=c(length(thresh), length(unique(dat.forecast$ENS))))
      #We want to pull from the Budburst_Model table for out GDD predictions
      ens <- unique(dat.forecast$ENS)
      for(i in 1:length(ens)){
        pred.array[,i] <- unlist(lapply(dat=dat.forecast[dat.forecast$ENS==ens[i],], FUN=calc.thresh.bud, VAR=Var, thresh))
      } 
    }
    pred.df <- data.frame(yday=pred.array)
    pred.long <- tidyr::gather(pred.df)
    colnames(pred.long) <- c("ens.ref", "yday")
    pred.long$state <- ST
    
    for(i in 1:length(unique(pred.long$ens.ref))){
      ref <- paste0("yday.",i)
      pred.long[pred.long$ens.ref == ref, "ens.num"] <- ens[as.numeric(i)]
    }
    
    prop.df <- aggregate(state~yday+ens.num, data=pred.long, FUN = length)
    colnames(prop.df) <- c("yday", "ens", "count")
    prop.df$ens.group <- as.numeric(gsub("\\..*","", prop.df$ens))
    prop.df$Proportion <- prop.df$count/nrow(pred.df)
    prop.df$state <- ST
    
    pred.st <- rbind(pred.st, prop.df)
    
    quant <- quantile(pred.array, c(0.025, 0.975))
    
    # Create some useful indices and labels
    dat.lim <- data.frame(ST, quant[1], quant[2],
                          mean(pred.array), min(pred.array), max(pred.array))
    
    colnames(dat.lim) <- c("state", "lb", "ub", "mean", "min", "max")
    row.names(dat.lim) <- NULL
    
    lim.comb[count, ] <- dat.lim
    count <- count +1
  

  }
  write.csv(pred.st, file.path(path.burst, paste0(MOD,"_Prop_Oak_Budburst_Prediction_",Sys.Date() ,".csv")), row.names = F)
  write.csv(lim.comb, file.path(path.burst, paste0(MOD, "_Oak_Prediciton_Summary_",Sys.Date() ,".csv")), row.names = F)
}


#Checking the dates we have done a forecast for so we can pick using the slider
past.fc <- list.files(path = file.path(path.met), full.names = F)
past.dates <- data.frame((gsub("\\..*", "", gsub(".*_", "", past.fc))), past.fc)
colnames(past.dates) <- c("Date", "File")


Model.name <- unique(comb.model$model)

model.df <- data.frame(Model.name)
model.df$Model.name <- factor(model.df$Model.name, levels = c("YDAY", "GDD5", "GTmean", "PTTGDD5"))

write.csv(model.df, file.path(path.temp, "Model_List.csv"), row.names = F)
write.csv(past.dates, file.path(path.temp, "Old_Forecast_List.csv"), row.names = F)
write.csv(dat.ghcn, file.path(path.temp, "Historical_Weather.csv"), row.names = F)
write.csv(dat.b, file.path(path.temp, "Full_Bur_Obs.csv"), row.names = F)


print("Data Organizaiton & Prediction Workflow Complete!")
