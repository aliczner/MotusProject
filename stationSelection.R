#code for tower/species selection from the WebDataSummary

library(lubridate)
library(dplyr)
library(tidyr)
library (sf)

Sys.setenv(TZ = "UTC")

station_pairs <- read.csv("observationSummary.csv") #78591 obs

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
#55283 obs

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

# Mark TRUE if distance is 30 - 965 km, otherwise
#FALSE_LESS if distance is less than 30 km 
# FALSE_MORE if distance is more than 965 km

filtered_pairs$valid_distance <- ifelse(filtered_pairs$distance_km < 30, 
                                        "FALSE_LESS",
                                        ifelse(filtered_pairs$distance_km > 965, 
                                               "FALSE_MORE", 
                                               "TRUE"))

#==================================================================
# station pair filtering by criteria
#==================================================================

filtered_df <- filtered_pairs %>%
  filter(
    # Keep if either the current OR the previous station is in GLWS
    (current_in_GLWS | previous_in_GLWS))

str(filtered_df) #15117 obs

table(filtered_pairs$valid_distance)
#Valid_distance
#TRUE 18623
#FALSE_LESS 27922
#FALSE_MORE 3646

write.csv(filtered_df, "StationPairsFiltered.csv")

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

#================================================================
# seeing the number of unique station pairs and unique stations
#================================================================

pair_counts_by_season <- filtered_df %>%
  group_by(season) %>%
  summarise(
    unique_pairs_count = n_distinct(previousStationID, stationID)
  )

(pair_counts_by_season)
#season           unique_pairs_count
#<chr>                         <int>
#  1 Fall Migration               3064
#2 Spring Migration               2313
#3 Summer                         1387
#4 Winter                          719

unique_stations_by_season <- filtered_df %>%
  select(season, previousStationID, stationID) %>%
  pivot_longer(
    cols = c(previousStationID, stationID), 
    values_to = "all_stations"
  ) %>%
  group_by(season) %>%
  summarise(
    unique_stations_count = n_distinct(all_stations)
  )
unique_stations_by_season
#season           unique_stations_count
#<chr>                            <int>
#  1 Fall Migration                   637
#2 Spring Migration                   543
#3 Summer                             378
#4 Winter                             300

#=================================================================
# getting less refined data
#================================================================

station_pairs_collapsed <- filtered_df %>%
#Group by the columns that define the paired stations
  group_by(
    previousStationID, 
    stationID, 
    season,
    lat, 
    lon, 
    lat_previous, 
    lon_previous,
    distance_km, 
  ) %>%
  
# Collapse the rows and count the number of unique species for each pair
  summarise(
    stationName = first(stationName),
    previousStationName = first(previousStationName),
    number_of_species = n_distinct(species),
    
    #Sum the total movements across all species for this pair of stations
    total_movements   = sum(number_of_movements, na.rm = TRUE),
    .groups = "drop"
  )
station_pairs_collapsed

write.csv(station_pairs_collapsed, "StationPairsCollapsed.csv")

# counting movements for individual stations

# Isolate the previous station data
prev_station <- filtered_df %>%
  select(
    station_ID = previousStationID, 
    stationName = previousStationName, 
    season, 
    species,
    number_of_movements
  )

# Isolate the current tower data
current_station <- filtered_df %>%
  select(
    station_ID = stationID, 
    stationName = stationName,
    season, 
    species, 
    number_of_movements
  )

# Combine them, group by station/season, and calculate summaries
individual_station_summary <- bind_rows(prev_station, 
                                       current_station) %>%
  group_by(station_ID, 
           season) %>%
  summarise(
    # Count how many unique species visited this specific station
    stationName = first(stationName),
    number_of_species = n_distinct(species),
    
    # Sum up all movements passing through this specific station
    total_movements   = sum(number_of_movements)
    )
write.csv(individual_station_summary, "StationIndividualSummary.csv")
