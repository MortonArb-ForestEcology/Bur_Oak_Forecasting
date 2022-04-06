#----------------------------------------------------------------------------------------------------------------------------------#
# Script by: Lucien Fitzpatrick
# Project: Bur Oak Forecasting (Quercus Quest)
# Purpose: This script creates figures showing the mother tree source populations locations and NPN observations
# Inputs: Filtered_NPN.csv which contains NPN bur oak observations with a column to flag if they are in the source area or not
#         Mother Tree Locations csv which contains the location of all mother trees
# Outputs: NPN_Bur_Oak_Observation_Locations.png which shows all bur oak observations in NPN
#          Bur_Oak_Source_Area_NPN_points.pdf which shows the convex hull around mother trees and nearby NPN points
#          Bur_Oak_Source_Area_Expanded_NPN_Range.pdf whcih shows the extended range of all the mother tree hulls
# Notes: This script only produces figures and isn't required in the overall workflow
#-----------------------------------------------------------------------------------------------------------------------------------#
library(ggplot2)
library(dplyr)
library(sf)
library(raster)

#Bur oak npn observations with the POTENTIAL to be filter. Not currently subset for figure creation
npn.filter <- read.csv("../data_processed/Filtered_NPN.csv")

#Loading in mother tree locations.
#Currently this was manually added but can pivot to pulling from a google drive
dat.moth <- read.csv("../data_raw/Mother_Tree_Locations.csv")
dat.moth$state <- substr(dat.moth$Unique.Code,1,2)
dat.moth <- dat.moth[!is.na(dat.moth$Latitude),]

#Setting the coordinate reference system (CRS) for the maps we are going to make
proj4string <-CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#Downloading maps for all the states and general us map
il_counties <- map_data("county", "illinois")
il_counties$state <- "IL"

mn_counties <- map_data("county", "minnesota") 
mn_counties$state <- "MN"

ok_counties <- map_data("county", "oklahoma") 
ok_counties$state <- "OK"

us_map <- map_data("usa")

#Converting the bur oak obs into a shapefile for graphing
npn.sf <- st_as_sf(npn.filter, coords = c("longitude", "latitude"),crs =proj4string )

#Figure of all bur oak observations
png(width=9, height=6, units="in", res=600, filename= paste0("../figures/NPN_Bur_Oak_Observation_Locations.png"))
ggplot()+
  geom_polygon(data= us_map, aes(long, lat, group=region), fill = "white", colour = "grey50")+
  geom_sf(data = npn.sf, color = "magenta", shape = 18, size = 3)+
  ggtitle(paste0("Bur Oak NPN Observation Locations"))+
  coord_sf()
dev.off()

all_counties <- rbind(il_counties,mn_counties,ok_counties)

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

#Pdf containing the mother tree hulls at their normal size
pdf("../figures/Bur_Oak_Source_Area_NPN_points.pdf")
for(ST in unique(all_counties$state)){
  dat.poly <- dat.moth[dat.moth$state == ST,]
  Moth.tree <- st_as_sf(dat.poly, coords = c("Longitude", "Latitude"),crs =proj4string )
  
  counties <- all_counties[all_counties$state == ST,]
  
  st.mcpdf <- full.mcpdf[full.mcpdf$id == ST,]
  mcp.shp <- st_as_sf(st.mcpdf, coords = c("long", "lat"),crs =proj4string )
  
  npn.plot <- npn.filter[npn.filter$state == ST,]
  npn.sf <- st_as_sf(npn.plot, coords = c("longitude", "latitude"),crs =proj4string )
  
  
  #For some reason Oklahoma is weird and not calculating properly. This manual fix helps for visuals but is suspcious
  if(ST == "OK"){
    st.mcpdf[4,] <- c("a", -97.52167, 34.98056, 6, FALSE, 1, "OK.1",3.84365860650178e-06)
    st.mcpdf$long <- as.numeric(st.mcpdf$long)
    st.mcpdf$lat <- as.numeric(st.mcpdf$lat)
    st.mcpdf$order <- as.integer(st.mcpdf$order)
  }
  
  plot_full <- ggplot()+
    geom_polygon(data= counties, aes(long, lat, group=subregion), fill = "white", colour = "grey50")+
    geom_sf(data = Moth.tree)+
    #geom_sf(data = pbuf.sf ,fill = "blue", color = "blue", alpha = 0.5)+
    geom_polygon(data = st.mcpdf, aes(x= long, y =lat, group = group),fill = "blue", color = "blue", alpha = 0.5)+
    geom_sf(data = npn.sf, color = "magenta", shape = 18, size = 3)+
    ggtitle(paste0("Bur Oak Mother Tree Locations (", ST, ")"))+
    coord_sf()
  print(plot_full)
  
  sub_counties <- counties[counties$long >= (min(dat.poly$Longitude)-1) & counties$long <= (max(dat.poly$Longitude)+1) & counties$lat >= (min(dat.poly$Latitude)-1) & counties$lat <= (max(dat.poly$Latitude)+1) ,  ]
  plot_zoom <- ggplot()+
    geom_polygon(data= sub_counties, aes(long, lat, group=subregion), fill = "white", colour = "grey50")+
    geom_sf(data = Moth.tree)+
    #geom_sf(data = pbuf.sf ,fill = "blue", color = "blue", alpha = 0.5)+
    geom_sf(data = npn.sf, color = "magenta", shape = 18, size = 3)+
    geom_polygon(data = st.mcpdf, aes(x= long, y =lat, group = group), fill = "blue", color = "blue", alpha = 0.5)+
    ggtitle(paste0("Bur Oak Mother Tree Locations (Zoomed in ", ST, ")"))+
    coord_sf()
  print(plot_zoom)
  
}
dev.off()


#Pdf containing the mother tree hulls at their expanded range
pdf("../figures/Bur_Oak_Source_Area_Expanded_NPN_Range.pdf")
for(ST in unique(all_counties$state)){
  
  #mother tree locations
  dat.poly <- dat.moth[dat.moth$state == ST,]
  Moth.tree <- st_as_sf(dat.poly, coords = c("Longitude", "Latitude"),crs =proj4string )
  
  counties <- all_counties[all_counties$state == ST,]
  
  #mother tree hulls
  st.mcpdf <- full.mcpdf[full.mcpdf$id == ST,]
  mcp.shp <- st_as_sf(st.mcpdf, coords = c("long", "lat"),crs =proj4string )
  
  #npn points
  npn.plot <- npn.filter[npn.filter$state == ST,]
  npn.sf <- st_as_sf(npn.plot, coords = c("longitude", "latitude"),crs =proj4string )

  
  #For some reason Oklahoma is weird and not calculating properly. This manual fix helps for visuals but is suspcious
  if(ST == "OK"){
    st.mcpdf[4,] <- c("a", -97.52167, 34.98056, 6, FALSE, 1, "OK.1",3.84365860650178e-06)
    st.mcpdf$long <- as.numeric(st.mcpdf$long)
    st.mcpdf$lat <- as.numeric(st.mcpdf$lat)
    st.mcpdf$order <- as.integer(st.mcpdf$order)
  }
  
  polygon_list = list(rbind(c(st.mcpdf[1,"long"] , st.mcpdf[1,"lat"]), c(st.mcpdf[2,"long"] , st.mcpdf[2,"lat"]), 
                            c(st.mcpdf[3,"long"] , st.mcpdf[3,"lat"]), c(st.mcpdf[4,"long"] , st.mcpdf[4,"lat"]),
                            c(st.mcpdf[5,"long"] , st.mcpdf[5,"lat"]), c(st.mcpdf[6,"long"] , st.mcpdf[6,"lat"]),
                            c(st.mcpdf[7,"long"] , st.mcpdf[7,"lat"]), c(st.mcpdf[8,"long"] , st.mcpdf[8,"lat"])))
  
  sf.poly <- st_polygon(polygon_list)
  
  #Create a buffer around the point. 1.5 means 1.5 decimal degrees which is ~160km or ~100m
  pbuf <-  st_buffer(sf.poly, 1.5)
  
  #This is a series of conversions necessary to get the buffer to a format where we can check the contained npn points
  pbuf.sfc <- st_sfc(pbuf)
  pbuf.sf <- st_sf(pbuf.sfc)
  pbuf.sf <- st_set_crs(pbuf.sf, proj4string)
  
  plot_full <- ggplot()+
    geom_polygon(data= counties, aes(long, lat, group=subregion), fill = "white", colour = "grey50")+
    geom_sf(data = Moth.tree)+
    geom_sf(data = pbuf.sf ,fill = "blue", color = "blue", alpha = 0.5)+
    geom_sf(data = npn.sf, color = "magenta", shape = 18, size = 3)+
    ggtitle(paste0("Bur Oak Mother Tree Locations (", ST, ")"))+
    coord_sf()
  print(plot_full)
  
}
dev.off()
