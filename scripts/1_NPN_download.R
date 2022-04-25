#----------------------------------------------------------------------------------------------------------------------------------#
# Script by: Lucien Fitzpatrick
# Project: Bur Oak Forecasting (Quercus Quest)
# Purpose: This script downloads bur oak NPN data, cleans them, and subsets by proximity to source populations
# Inputs: Mother Tree Locations csv which contains the location of all mother trees
# Outputs: Filtered_NPN.csv which contains NPN bur oak observations with a column to flag if they are in the source area or not
#          Filtered_NPN_Sum_Stats.csv which is a summary table of the NPN bur oak observations
#          Bur_Oak_Source_Area_NPN.pdf which contains an image of all NPN bur oak observations and their proximity to source
# Notes: Creating the spatial hull for OK requires a manual adjustment. It currently works but is an oddity and could be an issue down the line
#-----------------------------------------------------------------------------------------------------------------------------------#
library(rgdal)
library(ggplot2)
library(sf)
library(rnpn)

list <- npn_species()
list <- list[list$genus == "Quercus",]

#path.goog <- "G:/My Drive/Phenology - NPN/"
path.goog <- "/Volumes/GoogleDrive/My Drive/Phenology - Oaks - Hipp 2022/Phenology - NPN/"


#Loading in mother tree locations.
#Currently this was manually added but can pivot to pulling from a google drive
#dat.moth <- read.csv("../data_raw/Mother_Tree_Locations.csv")
dat.moth <- read.csv(paste0(path.goog, "Mother_Tree_Locations.csv"))

path.raw <- "../data_raw/"
if(!dir.exists(path.raw)) dir.create(path.raw, recursive=T, showWarnings = F)
write.csv(dat.moth, "../data_raw/Mother_Tree_Locations.csv", row.names = F)


dat.moth$state <- substr(dat.moth$Unique.Code,1,2)
dat.moth <- dat.moth[!is.na(dat.moth$Latitude),]

#Download all the bur oak budburst phenometrics
raw.npn <- npn_download_individual_phenometrics(phenophase_ids =c(371),  years=2000:2021, species_ids = 101, request_source="The Morton Arboretum")
raw.npn[raw.npn==-9999] <- NA

raw.npn$species <- as.factor(raw.npn$species)
raw.npn$species_id <- as.factor(raw.npn$species_id)
raw.npn$individual_id <- as.factor(raw.npn$individual_id)
raw.npn$phenophase_id <- as.factor(raw.npn$phenophase_id)
raw.npn$phenophase_description <- as.factor(raw.npn$phenophase_description)

#Removing any observations that don't have a NO within 10 days
dat.npn <- raw.npn[!is.na(raw.npn$numdays_since_prior_no) & raw.npn$numdays_since_prior_no<=10,]

#Only using the first budburst event for any individual and year. This is to remove repeat budburst observations
npn.first <- aggregate(first_yes_doy ~ site_id + latitude + longitude + elevation_in_meters + species_id + species + individual_id + phenophase_id + phenophase_description + first_yes_year + state, data=dat.npn, FUN=min)

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

#This is where you remove any flagged observations. Currently we don't have any flagged
npn.filter <- npn.filter[npn.filter$flag.4sig == F, ]

#Setting the coordinate reference system (CRS) for the maps we are going to make
proj4string <-CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#Creating the minimum convex polygon that captures all the points within it
Moth.mcp <- dat.moth[,c("Longitude", "Latitude", "state")]
coordinates(Moth.mcp) <- c("Longitude", "Latitude")
proj4string(Moth.mcp) <- CRS( "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0" )

#There is likley a way to do this without this package but I couldn't figure it out
#Maybe consult Shiven and Christy because of their PCA space hulls.
full.mcp <- adehabitatHR::mcp(Moth.mcp)

#Converting the mcp back into something ggplot can work with
full.mcp@data$id <- rownames(full.mcp@data)# create a data.frame from our spatial object
full.mcpdata <- fortify(full.mcp, region = "id")# merge the "fortified" data with the data from our spatial object
full.mcpdf <- merge(full.mcpdata, full.mcp@data,
                    by = "id")

#Defining the data frame that will contain the observations near source populations
npn.zones <- data.frame()
pdf("../figures/Bur_Oak_Source_Area_NPN.pdf")
for(ST in unique(full.mcpdf$id)){
    #Subsetting by state
    st.mcpdf <- full.mcpdf[full.mcpdf$id == ST,]
    #Converting to a shapefile
    mcp.shp <- st_as_sf(st.mcpdf, coords = c("long", "lat"),crs =proj4string )
    
    npn.plot <- npn.filter[npn.filter$state == ST,]
    #Converting to a shapefile
    npn.sf <- st_as_sf(npn.plot, coords = c("longitude", "latitude"),crs =proj4string )
    
    
    #For some reason Oklahoma is weird and not calculating properly. This manual fix helps for visuals but is suspcious
    #While this is messy, this part of the workflow isn't likely to change so I'm ok with it for now.
    if(ST == "OK"){
        st.mcpdf[4,] <- c("a", -97.52167, 34.98056, 6, FALSE, 1, "OK.1",3.84365860650178e-06)
        st.mcpdf$long <- as.numeric(st.mcpdf$long)
        st.mcpdf$lat <- as.numeric(st.mcpdf$lat)
        st.mcpdf$order <- as.integer(st.mcpdf$order)
    }
    
    #Again could be automated but didn't seem worht the time investment for data that won't change
    polygon_list = list(rbind(c(st.mcpdf[1,"long"] , st.mcpdf[1,"lat"]), c(st.mcpdf[2,"long"] , st.mcpdf[2,"lat"]), 
                              c(st.mcpdf[3,"long"] , st.mcpdf[3,"lat"]), c(st.mcpdf[4,"long"] , st.mcpdf[4,"lat"]),
                              c(st.mcpdf[5,"long"] , st.mcpdf[5,"lat"]), c(st.mcpdf[6,"long"] , st.mcpdf[6,"lat"]),
                              c(st.mcpdf[7,"long"] , st.mcpdf[7,"lat"]), c(st.mcpdf[8,"long"] , st.mcpdf[8,"lat"])))
    
    #Make the polygon a shapefile
    sf.poly <- st_polygon(polygon_list)
    
    #Create a buffer around the point. 1.5 means 1.5 decimal degrees which is ~160km or ~100m
    pbuf <-  st_buffer(sf.poly, 1.5)
    
    #This is a series of conversions necessary to get the buffer to a format where we can check the contained npn points
    pbuf.sfc <- st_sfc(pbuf)
    pbuf.sf <- st_sf(pbuf.sfc)
    pbuf.sf <- st_set_crs(pbuf.sf, proj4string)
    
    #Flagging the points in the buffer with a 1
    npn.plot$close.check <- sapply(st_intersects(npn.sf,pbuf.sf), function(z) if (length(z)==0) FALSE else z[TRUE])
    
    npn.zones <- rbind(npn.zones, npn.plot)
}

path.out <- "../data_processed/"

if(!dir.exists(path.out)) dir.create(path.out, recursive=T, showWarnings = F)

#Creating a csv that contains a column flagging which obs are near enough to a mother tree
write.csv(npn.zones, paste0(path.out, "Filtered_NPN.csv"), row.names = F)

#Filtering out the obs that are too far
npn.filter <- npn.zones[npn.zones$close.check == 1,]

#Manual filling for a table
State <- c("IL", "MN", "OK")
n_sites <- c(length(unique(npn.filter[npn.filter$state == "IL", "site_id"])), length(unique(npn.filter[npn.filter$state == "MN", "site_id"])), length(unique(npn.filter[npn.filter$state == "OK", "site_id"])))
n_trees <- c(length(unique(npn.filter[npn.filter$state == "IL", "individual_id"])), length(unique(npn.filter[npn.filter$state == "MN", "individual_id"])), length(unique(npn.filter[npn.filter$state == "OK", "individual_id"])))
n_obs <- c(nrow(npn.filter[npn.filter$state == "IL", ]), nrow(npn.filter[npn.filter$state == "MN", ]), nrow(npn.filter[npn.filter$state == "OK", ]))

sum.tab <- data.frame(State, n_sites, n_trees, n_obs) 

#Writing the summary stats of the npn obs used for modeling
write.csv(sum.tab, paste0(path.out, "Filtered_NPN_Sum_Stats.csv"), row.names = F)
