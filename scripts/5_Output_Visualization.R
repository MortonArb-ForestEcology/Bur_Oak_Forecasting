#----------------------------------------------------------------------------------------------------------------------------------#
# Script by: Lucien Fitzpatrick
# Project: Bur Oak Forecasting (Quercus Quest)
# Purpose: This script downloads the NPN data and cleans them
# Inputs: 
# Outputs:
# Notes: 
#-----------------------------------------------------------------------------------------------------------------------------------#
library(ggplot2)
library(cowplot)
library(tidyr)

#Reading in budburst model
path.mod <- "../data_processed/model_output/"

wide.mod <- read.csv(paste0(path.mod, "Quercus_Quest_TT_model_bur.csv"))

b.model <- gather(wide.mod, state, THRESH, IL:OK, factor_key=TRUE)

convg.df <- read.csv(file.path(path.mod, paste0("Quercus_Quest_convergence.csv")))
good.state <- convg.df[convg.df$burst.converge < 1.05, "state"]

safe_colorblind_palette <- c("#88CCEE",  "#DDCC77", "#CC6677", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

b.model$state <- factor(b.model$state, levels = c("MN", "IL", "OK"))

#set.seed(903)
#Taking a  random sample of 1000 pulls
#b.model <- do.call(rbind,lapply(split(b.model, b.model$state),function(x) x[sample(nrow(x), 1000), ]))

#rownames(b.model) <- NULL

png(width=9, height=6, units="in", res=600, filename= paste0("../figures/TT_Bur_Oak_GDD5_dist.png"))
ggplot() + 
  geom_density(data=b.model, aes(x=THRESH, fill = state, color = state), adjust=3.5, alpha=0.5) +
  scale_x_continuous(name="GDD5.cum", expand=c(0,0))  +
  scale_y_continuous(name="Probability of Bud Burst" ,expand=c(0,0)) +
  scale_color_manual(values = safe_colorblind_palette)+
  scale_fill_manual(values = safe_colorblind_palette)+
  theme_bw() +
  ggtitle(paste0("Budburst GDD5 Distribution"))
dev.off()

dir.met <- c("../../Phenology_Forecasting/scripts/Morton_Bloom_Forecast/data_raw/meteorology/")

#Reading in our latest forecast
dat.forecast <- read.csv(file.path(dir.met, paste0("Mortonarb_daily_FORECAST-READY-LONGRANGE_2022-03-28.csv")))
vars.agg <- c("TMEAN", "GDD0.cum", "GDD5.cum", "CDD0.cum", "CDD2.cum")
ens.forecast <- list()

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

fc.sp <- as.data.frame(do.call(rbind, ens.forecast))
fc.sp$VAR <- gsub("\\..*","",(row.names(fc.sp)))
rownames(fc.sp) <- NULL
colnames(fc.sp) <- c("PRED.DATE", "YDAY", "TYPE", "mean", "min", "max", "VAR")

#write.csv(fc.sp, file.path(path.temp, "meteorology", paste0("Forecast_data_", Sys.Date(),".csv")), row.names = F)

calc.bud <- function(dat, VAR, THRESH){
  min(dat[which(dat[,VAR] >= THRESH),"YDAY"])
}

#Need to optomize this later
pred.st <- data.frame()

lim.comb <- data.frame(matrix(ncol = 6, nrow = length(unique(b.model$state))))
colnames(lim.comb) <- c("state", "lb", "ub", "mean", "min", "max")

count <- 1
for(ST in unique(b.model$state)){
  #Pulling out the gdd5.cum vlaues
  
  thresh <- b.model[b.model$state == ST, "THRESH"]
  thresh <- bur.mod$THRESH
  
  pred.array <- array(dim=c(length(thresh), length(unique(dat.forecast$ENS))))
  #We want to pull from the Budburst_Model table for out GDD predictions
  ens <- unique(dat.forecast$ENS)
  for(i in 1:length(ens)){
    pred.array[,i] <- unlist(lapply(dat=dat.forecast[dat.forecast$ENS==ens[i],], FUN=calc.bud, VAR="GDD5.cum", thresh))
  }
  #pred.df <- data.frame(yday=as.vector(pred.array))
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

path.burst <- "../data_processed/"
write.csv(pred.st, file.path(path.burst, paste0("GDD5_Prop_Oak_Budburst_Prediction.csv")), row.names = F)
write.csv(lim.comb, file.path(path.burst, paste0("GDD5_Oak_Prediciton_Summary.csv")), row.names = F)

#pred.st <- read.csv(file.path(path.burst, paste0("GDD5_Prop_Oak_Budburst_Prediction.csv")))
#lim.comb <- read.csv(file.path(path.burst, paste0("GDD5_Oak_Prediciton_Summary.csv")))

# Creating some day and axis labels
day.labels <- data.frame(Date=seq.Date(as.Date("2021-01-01"), as.Date("2021-12-31"), by="month"))
day.labels$yday <- lubridate::yday(day.labels$Date)
day.labels$Text <- paste(lubridate::month(day.labels$Date, label=T), lubridate::day(day.labels$Date))
summary(day.labels)

day.labels2 <- data.frame(Date=seq.Date(as.Date("2021-01-03"), as.Date("2021-7-01"), by="day"))
day.labels2$yday <- lubridate::yday(day.labels2$Date)
day.labels2$Text <- paste(lubridate::month(day.labels2$Date, label=T), lubridate::day(day.labels2$Date))
summary(day.labels2)

pred.st$state <- factor(pred.st$state, levels = c("MN", "IL", "OK"))
  
png(width=9, height=6, units="in", res=600, filename= paste0("../figures/TT_Bur_Oak_YDAY_dist.png"))
ggplot() + 
  geom_density(data=pred.st, aes(x=yday, color = state, fill = state), adjust=3.5, alpha=0.5) +
  scale_x_continuous(name="Day of Year", expand=c(0,0), breaks=day.labels2$yday[seq(8, nrow(day.labels2), by=7)], labels=day.labels2$Text[seq(8, nrow(day.labels2), by=7)])  +
  scale_y_continuous(name="Probability of Bud Burst" ,expand=c(0,0), breaks=day.labels2$yday[seq(8, nrow(day.labels2), by=7)], labels=day.labels2$Text[seq(8, nrow(day.labels2), by=7)]) +
  scale_color_manual(values = safe_colorblind_palette)+
  scale_fill_manual(values = safe_colorblind_palette)+
  theme_bw() +
  ggtitle(paste0("Budburst DOY Distribution"))
dev.off()
