## code to create a shapefile around radar stations
## may use the resulting polygon to select motus towers in range of radar

library(sf)

radar.data <- read.csv ("./radarFiles/RadarStationLocations.csv")

radar.ON <- subset(radar.data, Province == "ON")

radar_sf <- st_as_sf(radar.ON, 
                     coords = c("Longitude", "Latitude"), 
                     crs = 4326)

#reproject to lambert conformal conic
radar_projected <- st_transform(radar_sf, crs = 3348)

#create buffers around radar to make the polygon
# using three sizes to choose from
# 1) 240 km, the doppler range
# 2) 100 km, closer increases optimal detections
# 3) 50 km, closest but still allows some range
radar_buffer_240 <- st_buffer(radar_projected, 
                              dist = 240000)
radar_buffer_100 <- st_buffer(radar_projected, 
                              dist = 100000)
radar_buffer_50 <- st_buffer(radar_projected, 
                             dist = 50000)
#

#dissolve buffers into one polygon
buffer_union_240 <- st_union(radar_buffer_240)
buffer_union_100 <- st_union(radar_buffer_100)
buffer_union_50 <- st_union(radar_buffer_50)

#view the polygon and the radar locations
library(mapview)

mapview(radar_projected) + 
  mapview(buffer_union_240) +
  mapview(buffer_union_100) +
  mapview(buffer_union_50)

#save dissolved buffers
st_write(buffer_union_240, "radarFiles/buffer_union_240.shp")
st_write(buffer_union_100, "radarFiles/buffer_union_100.shp")
st_write(buffer_union_50, "radarFiles/buffer_union_50.shp")

