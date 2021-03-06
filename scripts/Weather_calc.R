#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick
# Project: Collections phenology vulnerability
# Purpose: This script serves to calculate weather covariates of interest. 
#         Currently Growing degree days at 5C, Growing degree days at 0C, Number of chill days, and Growing season mean temperature
# Inputs: Lubridate package
# Outputs:This function will take a data frame of daily weather data and produce the following summary statistics
#         GDD5 = Growing degree days at 5 degrees C 
#         GDD0 = Growing degree days at 0 degrees C
#         NCD = Number of chilling days 
#         GTmean = Growing season mean temperature
# Notes: The defaults for this funcion are
#       Starting year of interest                       y_start = 1975
#       Ending year of interest                         y_end = 2019
#       Julian yday for start of growing season         g_start = 1
#       Julian yday for end of growing season           g_end = 60
#
#-----------------------------------------------------------------------------------------------------------------------------------#
# Calculating the Tmean for the growing season of that year
weather_calc <- function(met.all, g_start=1, g_end=90){
  
  
  #Calculating mean temperature, growing degree days using 5C, and gorwing degree days using 0c
  met.all$TMEAN <- (met.all$tmax..deg.c. + met.all$tmin..deg.c.)/2
  met.all$GDD5 <- ifelse(met.all$TMEAN>5, met.all$TMEAN-5, 0)
  met.all$GDD0 <- ifelse(met.all$TMEAN>0, met.all$TMEAN, 0)
  met.all$CDD0 <- ifelse(met.all$tmin..deg.c.<0, 0-met.all$tmin..deg.c., 0)
  
  #Creating empty columns to fill in the loop
  met.all$GDD5.cum <- NA
  met.all$GDD0.cum <- NA
  met.all$NCD <- NA
  met.all$PTTGDD.cum <- NA
  met.all$CDD0.cum <- NA
  
  #Setting the beginning and end (using julian yday) of the growing season for our Growing season mean
  g_start <- 1
  g_end <- 90
  
  
  # Calculate the cumulative growing degree days for each day/year
  for(YR in unique(met.all$year)){
    dat.tmp <- met.all[met.all$year==YR, ]
    dat.gtmean <- dat.tmp[(dat.tmp$yday>=g_start & dat.tmp$yday<=g_end), ]
    gtmean <- mean(dat.gtmean$TMEAN, na.rm = TRUE)
    dat.tmp$Date <- as.Date(paste(dat.tmp$year, dat.tmp$yday, sep="-"), format="%Y-%j")
    gdd5.cum=0; gdd0.cum=0; pttgdd.cum=0; cdd0.cum=0;
    d5.miss = 0; d0.miss=0; ptt.miss=0; c0.miss = 0;
    ncd = 0
    for(i in 1:nrow(dat.tmp)){
      if(is.na(dat.tmp$GDD5[i]) & d5.miss<=3){ 
        d5.miss <- d5.miss+1 # Let us miss up to 3 consecutive days
        gdd5.cum <- gdd5.cum+0
      } else {
        d5.miss = 0 # reset to 0
        gdd5.cum <- gdd5.cum+dat.tmp$GDD5[i] 
      }
      
      if(is.na(dat.tmp$GDD0[i]) & d0.miss<=3){
        d0.miss <- d0.miss+1 # Let us miss up to 3 consecutive days
        gdd0.cum <- gdd0.cum+0
      } else {
        d0.miss = 0 # reset to 0
        gdd0.cum <- gdd0.cum+dat.tmp$GDD0[i] 
      }
      
      if(is.na(dat.tmp$CDD0[i]) & c0.miss<=3){
        c0.miss <- c0.miss+1 # Let us miss up to 3 consecutive days
        cdd0.cum <- cdd0.cum+0
      } else {
        c0.miss = 0 # reset to 0
        cdd0.cum <- cdd0.cum+dat.tmp$CDD0[i] 
      }
      
      if(!is.na(dat.tmp$TMEAN[i]) & dat.tmp$TMEAN[i] < 0){
        ncd <- ncd + 1
      }
      
      if(is.na(dat.tmp$GDD5[i]) & ptt.miss<=3){ 
        ptt.miss <- ptt.miss+1 # Let us miss up to 3 consecutive days
        pttgdd.cum <- pttgdd.cum+0
      } else {
        ptt.miss = 0 # reset to 0
        pttgdd.cum <- pttgdd.cum+(dat.tmp$GDD5[i] * (dat.tmp$dayl..s.[i]/86400)) #86400 is converting seconds of daylength into proportion of day 
      }
      
      dat.tmp[i,"GDD5.cum"] <- gdd5.cum
      dat.tmp[i,"GDD0.cum"] <- gdd0.cum
      dat.tmp[i,"CDD0.cum"] <- cdd0.cum
      dat.tmp[i, "NCD"] <- ncd
      dat.tmp[i,"PTTGDD.cum"] <- pttgdd.cum
      dat.tmp[i, "GTmean"] <- gtmean
    }
    
    met.all[met.all$year==YR, "GDD5.cum"] <- dat.tmp$GDD5.cum
    met.all[met.all$year==YR, "GDD0.cum"] <- dat.tmp$GDD0.cum
    met.all[met.all$year==YR, "CDD0.cum"] <- dat.tmp$CDD0.cum
    met.all[met.all$year==YR, "NCD"] <- dat.tmp$NCD
    met.all[met.all$year==YR, "GTmean"] <- dat.tmp$GTmean
    met.all[met.all$year==YR, "PTTGDD.cum"] <- dat.tmp$PTTGDD.cum
    met.all[met.all$year==YR, "Date"] <- as.Date(dat.tmp$Date)
  }
  return(met.all)
}