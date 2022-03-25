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

#Reading in budburst model
path.mod <- "../data_processed/model_output/"

wide.mod <- read.csv(paste0(path.mod, "Quercus_Quest_TT_model_bur.csv"))

b.model <- gather(wide.mod, state, THRESH, IL:OK, factor_key=TRUE)

convg.df <- read.csv(file.path(path.mod, paste0("Quercus_Quest_convergence.csv")))
good.state <- convg.df[convg.df$burst.converge < 1.05, "state"]

#set.seed(903)
#Taking a  random sample of 1000 pulls
#b.model <- do.call(rbind,lapply(split(Budburst_Model, Budburst_Model$state),function(x) x[sample(nrow(x), 1000), ]))

#rownames(b.model) <- NULL

pdf("../figures/TT_Bur_Oak_GDD5_dist.pdf")
for(ST in unique(b.model$state)){
  
  quant <- quantile(b.model[b.model$state == ST,"THRESH"], c(0.025, 0.975))
  
  plot.gdd <- ggplot() + 
    geom_density(data=b.model[b.model$state == ST,], aes(x=THRESH), fill = "darkgreen", adjust=3.5, alpha=0.5) +
    geom_vline(aes(xintercept=quant[1]), linetype="dashed") +
    geom_vline(aes(xintercept=quant[2]), linetype="dashed") +
    scale_x_continuous(name="GDD5.cum", expand=c(0,0))  +
    scale_y_continuous(name="Probability of Bud Burst" ,expand=c(0,0)) +
    theme_bw() +
    guides(fill="none") +
    ggtitle(paste0(ST, ": Budburst DOY 95% C.I."))+
    theme(legend.position = c(0.5, 0.2),
          legend.title=element_blank(),
          legend.background = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
  
  title <- ggdraw() + 
    draw_label(
      paste0(ST, " sourced trees predicted budburst timing at TMA 95% CI GDD5.cum range ", round(quant[1]),
             " to ", round(quant[2])),
      fontface = 'bold',
      x = 0,
      hjust = 0, 
      size = 10,
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  final.plot <- plot_grid(
    title, plot.gdd,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  print(final.plot)
}
dev.off()    


dir.met <- c("../../Phenology_Forecasting/scripts/Morton_Bloom_Forecast/data_raw/meteorology/")

#Reading in our latest forecast
dat.forecast <- read.csv(file.path(dir.met, paste0("Mortonarb_daily_FORECAST-READY-LONGRANGE_2022-03-18.csv")))
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
write.csv(pred.st, file.path(path.burst, paste0("Prop_Oak_Budburst_Prediction.csv")), row.names = F)
write.csv(lim.comb, file.path(path.burst, paste0("Oak_Prediciton_Summary.csv")), row.names = F)

# Creating some day and axis labels
day.labels <- data.frame(Date=seq.Date(as.Date("2021-01-01"), as.Date("2021-12-31"), by="month"))
day.labels$yday <- lubridate::yday(day.labels$Date)
day.labels$Text <- paste(lubridate::month(day.labels$Date, label=T), lubridate::day(day.labels$Date))
summary(day.labels)

calc.bud <- function(dat, VAR, THRESH){
  min(dat[which(dat[,VAR] >= THRESH),"YDAY"])
}

day.labels2 <- data.frame(Date=seq.Date(as.Date("2021-01-03"), as.Date("2021-7-01"), by="day"))
day.labels2$yday <- lubridate::yday(day.labels2$Date)
day.labels2$Text <- paste(lubridate::month(day.labels2$Date, label=T), lubridate::day(day.labels2$Date))
summary(day.labels2)

pdf("../figures/TT_Bur_Oak_Yday_dist.pdf")
for(ST in unique(pred.st$state)){
  plot.yday <- ggplot() + 
    geom_density(data=pred.st[pred.st$state == ST,], aes(x=yday), fill = "darkgreen", adjust=3.5, alpha=0.5) +
    geom_vline(data=lim.comb[lim.comb$state == ST,], aes(xintercept=lb), linetype="dashed") +
    geom_vline(data=lim.comb[lim.comb$state == ST,], aes(xintercept=ub), linetype="dashed") +
    scale_x_continuous(name="Day of Year", expand=c(0,0), breaks=day.labels2$yday[seq(8, nrow(day.labels2), by=7)], labels=day.labels2$Text[seq(8, nrow(day.labels2), by=7)])  +
    scale_y_continuous(name="Probability of Bud Burst" ,expand=c(0,0), breaks=day.labels2$yday[seq(8, nrow(day.labels2), by=7)], labels=day.labels2$Text[seq(8, nrow(day.labels2), by=7)]) +
    theme_bw() +
    guides(fill="none") +
    ggtitle(paste0(ST, ": Budburst DOY 95% C.I."))+
    theme(legend.position = c(0.5, 0.2),
          legend.title=element_blank(),
          legend.background = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
  
  title <- ggdraw() + 
    draw_label(
      paste0(ST, " sourced trees predicted budburst timing at TMA 95% CI date range ", as.Date((lim.comb[lim.comb$state == ST,"lb"]-1), origin = paste0(lubridate::year(Sys.Date()),"-01-01")) ,
             " to ", as.Date((lim.comb[lim.comb$state == ST,"ub"]-1), origin = paste0(lubridate::year(Sys.Date()),"-01-01"))),
      fontface = 'bold',
      x = 0,
      hjust = 0, 
      size = 10,
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  final.plot <- plot_grid(
    title, plot.yday,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  print(final.plot)
}
dev.off()     
