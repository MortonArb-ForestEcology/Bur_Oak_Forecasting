#----------------------------------------------------------------------------------------------------------------------------------#
# Script by: Lucien Fitzpatrick
# Project: Bur Oak Forecasting (Quercus Quest)
# Purpose: This script downloads the NPN data and cleans them
# Inputs: 
# Outputs:
# Notes: 
#-----------------------------------------------------------------------------------------------------------------------------------#
npn.zones <- read.csv("../data_processed/Filtered_NPN.csv")
npn.bur <- npn.zones[npn.zones$close.check == 1,]

path.daymet <- "../data_raw/DAYMET"
if(!dir.exists(path.daymet)) dir.create(path.daymet, recursive=T)
if(!dir.exists("../data_processed")) dir.create("../data_processed", recursive=T)

#Creating combined data frame so that we can donwload from all sites with one go considering the amount of overlap
NPN.pts <- aggregate(year~site_id+latitude+longitude, data=npn.bur, 
                     FUN=min)
names(NPN.pts)[names(NPN.pts=="year")] <- "yr.start"
NPN.pts$yr.end <- aggregate(year~site_id+latitude+longitude, data=npn.bur, 
                            FUN=max)[,"year"]
summary(NPN.pts)


#Writing the csv file of lat and longs because daymetr batch function needs to read a file instead of a dataframe
write.csv(NPN.pts, file.path(path.daymet, "NPN_points.csv"), row.names=FALSE)


#Downloading all of the damet data for each point. Internal =TRUE means it creates a nested list. Set false to actually download a file
#At its current size (991 sites) +.this can take about an hour to run
lat.list <- daymetr::download_daymet_batch(file_location = file.path(path.daymet, "NPN_points.csv"),
                                           start = 1980,
                                           end = 2019,
                                           internal = T)


#removing failed downloads 
lat.list <- lat.list[sapply(lat.list, function(x) is.list(x))]


# This gives us a list with one layer per site (I think)
length(lat.list)
names(lat.list) <- NPN.pts$site # Giving the different layers of the list the site names they correspond to

#Lets look at the structure of what we are given
summary(lat.list)

list.met <- list()
for(i in seq_along(lat.list)){
  list.met[[i]] <- data.frame(site=NPN.pts$site_id[i], latitude=NPN.pts$latitude[i], longitude=NPN.pts$longitude[i], lat.list[[i]]$data)
}
names(list.met) <-  NPN.pts$site

rm(lat.list) # Removing lat.list to save memory
#Reading in our function for calculating weather statistics of interest
source(file.path("Weather_calc.R"))

#Running our function to calculate weather statistics. Default year range is 1975-2019. Growing seaosn is yday 1 to 120
list.met<- lapply(list.met, weather_calc)

lat.calc <- dplyr::bind_rows(list.met)
summary(lat.calc)

write.csv(lat.calc, "../data_processed/Daymet_clean_data.csv", row.names=F)


library(dplyr)
library(tidyr)
dir.create("../data_processed/model_output", recursive = T, showWarnings = F)

npn.bur$date <- as.Date(paste(npn.bur$year, npn.bur$yday), format="%Y %j")
npn.bur$species <- as.character(npn.bur$species)

#Reading in the weather data
#lat.calc <- read.csv("../data_processed/Daymet_clean_data.csv")
lat.calc$Freeze <- ifelse(lat.calc$tmin..deg.c.<0, 1, 0)
lat.calc$site <- as.character(lat.calc$site)
# lat.calc$date <- as.Date(lat.calc$date)

# Merging met into the phenology data
dat.npnbud <- merge(npn.bur, lat.calc[,c("year", "yday", "site", "TMEAN", "GDD5.cum", "GDD0.cum", "CDD5.cum","GTmean")], by.x = c("year", "yday", "site_id"), by.y = c("year", "yday", "site"))
summary(dat.npnbud)
summary(dat.npnbud[!dat.npnbud$flag.3sig,])

# names(arb.burst)[names(arb.burst) %in% names(dat.npnbud)]
cols.keep <- c("species", "site_id", "state", "latitude", "longitude", "individual_id", "date", "year", "yday", "TMEAN" ,"GDD5.cum", "GDD0.cum", "CDD5.cum", "GTmean", "flag.3sig", "flag.4sig")

bud.all <- dat.npnbud[, cols.keep]

write.csv(bud.all, "../data_processed/Full_Bur_Obs.csv", row.names = F)

