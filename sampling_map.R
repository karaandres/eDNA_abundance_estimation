
library(ggplot2)
library(dplyr)
library(stringr)
library(ggmap)
library(maps)
library(mapdata)
library(rgdal)
library(maptools)
library(mapproj)
library(sp)
library(raster)
library(grDevices)
library(RColorBrewer)

sampling_locations <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/goby_abundance_eDNA_sample_sheet_ForPete.csv", header=TRUE)
lakes <- readOGR("/Users/kbja10/Downloads/ne_10m_lakes") # file path to lake shapefile
great_lakes <- subset(lakes, name_alt=="Great Lakes") # make object for Great Lakes
vol.25 <- sampling_locations[sampling_locations$Expected_Volume_L==0.25,]
vol100 <- sampling_locations[sampling_locations$Expected_Volume_L==100,]
map(database = "state", regions = c("Wisconsin","Illinois","Indiana","Michigan","Ohio","Pennsylvania","New York"), 
    col = c("white"), fill = TRUE, border = "lightgray")
plot(great_lakes, col = "blue4", add = TRUE) # add Great Lakes 

points(x = vol.25$Longitude, y = vol.25$Latitude, col = "black", pch=21, bg="lightgray")
points(x = vol100$Longitude, y = vol100$Latitude, col = "black", pch=21, bg="cadetblue", cex=2)

map.scale(y=41, x=-92.5, ratio=FALSE, relwidth = 0.3) # add a scale bar
box() # add a box around the map

