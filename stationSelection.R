#code for tower/species selection from the WebDataSummary

library(lubridate)
library(dplyr)
library(tidyr)
library (sf)

Sys.setenv(TZ = "UTC")

station_pairs <- read.csv("station_pairs_summary.csv") #24969 obs

#====================================================================
# station pair filtering column creation for later filtering
#===================================================================

#filter out same station pairs and station unknown

filtered_pairs <- station_pairs %>% 
  filter(
    previousStationName != stationName,
    previousStationName != "Unknown",
    previousStationName != "Unknown station"
  )
#21940 obs

#loading great lakes watershed polygon, contains subbasins for each lake
GLWatershed <- st_read("./greatlakes_subbasins/greatlakes_subbasins.shp")
#checking CRS
st_crs(GLWatershed) #ESPG 6269

#transform both the polygon and station pairs to lambert conformal conic

#need to run current and previous stations separate;y
current_sf <- st_as_sf(filtered_pairs, 
                       coords = c("lon", "lat"), 
                       crs = 4326) %>%
  st_transform(crs = 3348)

previous_sf <- st_as_sf(filtered_pairs, 
                        coords = c("lon_previous", "lat_previous"), 
                        crs = 4326) %>%
  st_transform(crs = 3348)

GLWS_proj <- st_transform(GLWatershed, crs = 3348)

# check if either of the stations are within the GL watershed, 
filtered_pairs$current_in_GLWS  <- lengths(st_intersects(current_sf, 
                                                         GLWS_proj)) > 0
filtered_pairs$previous_in_GLWS <- lengths(st_intersects(previous_sf, 
                                                         GLWS_proj)) > 0

# calculate distance between points
filtered_pairs$distance_km <- as.numeric(st_distance(current_sf, 
                                                     previous_sf, 
                                                     by_element = TRUE)
                                         ) / 1000 #convert to km

# Mark TRUE if distance is 30 - 965 km, otherwise FALSE
# < 30 km could be in range of another tower
# > 965 km is beyond what a bird could fly in a day during migration

filtered_pairs$valid_distance <- filtered_pairs$distance_km >= 30 & 
                                  filtered_pairs$distance_km <= 965

#==================================================================
# station pair filtering by criteria
#==================================================================

filtered_df <- filtered_pairs %>%
  filter(
    # Keep if either the current OR the previous station is in GLWS
    (current_in_GLWS | previous_in_GLWS),
    
    # Keep only if paired points falls within 30km to 965km distance
    valid_distance == TRUE
  )

str(filtered_df) #9214 obs

write.csv(filtered_df, "StationPairsFilterd.csv")

### creating a map of the points

#separate out spring migration

spring <- subset(filtered_df, season == "Spring Migration") 
#2910 obs

current_map_GLWS <- st_as_sf(spring, 
                             coords = c("lon", "lat"), 
                             crs = 4326) %>% 
  st_transform(crs = 3348)

prev_map_GLWS <- st_as_sf(spring, 
                          coords = c("lon_previous", "lat_previous"), 
                          crs = 4326) %>% 
  st_transform(crs = 3348)

mapview(GLWS_proj, 
        col.regions = "gray", 
        alpha.regions = 0.3, 
        layer.name = "GLWS") +
  mapview(current_map_GLWS, 
          col.regions = "blue", 
          layer.name = "Current Station") +
  mapview(prev_map_GLWS, 
          col.regions = "orange",
          layer.name = "Previous Station")
