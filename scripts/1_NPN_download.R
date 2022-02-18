#----------------------------------------------------------------------------------------------------------------------------------#
# Script by: Lucien Fitzpatrick
# Project: Bur Oak Forecasting (Quercus Quest)
# Purpose: This script downloads the NPN data and cleans them
# Inputs: 
# Outputs:
# Notes: 
#-----------------------------------------------------------------------------------------------------------------------------------#
# # NOTE: requires the github version of rnpn; CRAN won't work

#devtools::install_github("bluegreen-labs/phenor", force = TRUE)
library(phenor)

#Download all the bur oak budburst phenometrics
#Including 2021 observations seems to break the formatting script. For now I'm just ignoring that year but we should come back
raw.npn <- pr_dl_npn(species = 101, phenophase = 371, start = "2000-01-01", end = "2021-01-01")
raw.npn[raw.npn==-9999] <- NA

raw.npn$species <- as.factor(raw.npn$species)
raw.npn$species_id <- as.factor(raw.npn$species_id)
raw.npn$individual_id <- as.factor(raw.npn$individual_id)
raw.npn$phenophase_id <- as.factor(raw.npn$phenophase_id)
raw.npn$phenophase_description <- as.factor(raw.npn$phenophase_description)

#Removing any observations that don't have a NO within 10 days
dat.npn <- raw.npn[!is.na(raw.npn$numdays_since_prior_no) & raw.npn$numdays_since_prior_no<=10,]

#Only using the first budburst event for any individual and year. This is to remove repeat budburst observations
npn.first <- aggregate(first_yes_doy ~ site_id + latitude + longitude + elevation_in_meters + species_id + species + individual_id + phenophase_id + phenophase_description + first_yes_year, data=dat.npn, FUN=min)

# Looking at jsut sp
# Filter by the summer equinox: June 21 (day 172)
npn.filter <- npn.first[npn.first$first_yes_doy<172,]
npn.filter$year <- npn.filter$first_yes_year
npn.filter$yday <- npn.filter$first_yes_doy
summary(npn.filter)

#Checking for outliers to flag
summary(npn.filter$species)
npn.filter$flag.3sig <- NA
npn.filter$flag.4sig <- NA
for(SPP in unique(npn.filter$species)){
    row.spp <- which(npn.filter$species==SPP)
    
    spp.mean <- mean(npn.filter$yday[row.spp])
    spp.sd <- sd(npn.filter$yday[row.spp])
    
    npn.filter$flag.3sig[row.spp] <- ifelse(npn.filter$yday[row.spp]<spp.mean-3*spp.sd | npn.filter$yday[row.spp]>spp.mean+3*spp.sd, T, F)
    npn.filter$flag.4sig[row.spp] <- ifelse(npn.filter$yday[row.spp]<spp.mean-4*spp.sd | npn.filter$yday[row.spp]>spp.mean+4*spp.sd, T, F)
    
    
}

summary(npn.filter)

#This is where you would remove any flagged observations. Currently we don't have any flagged
npn.filter <- npn.filter[npn.filter$flag.4sig == F, ]

#Converting into phenoR's preferred format
form.npn <- pr_fm_npn(npn.filter)

#Useful for identifying the parameters of different models
path <- sprintf("%s/extdata/parameter_ranges.csv",path.package("phenor"))
par_ranges <- read.table(path,
                         header = TRUE,
                         sep = ",")

#Way to automatically calculate ideal parameters.
par.fit <- pr_fit_parameters(par = NULL, data = form.npn, cost = rmse, model = "TT", method = "GenSA", lower = c(1,-5,0), upper = c(365,10,2000), control = NULL)

#Manually chosen parameters for testing
par = as.vector(c(150, -5, 400))

#Converting into the "flat" format needed for modelling 
dat.mod <- pr_flatten(form.list)

mod.out <- TT(par = par.fit$par, data= dat.mod)

