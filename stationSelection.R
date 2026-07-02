#code for tower/species selection from the WebDataSummary

library(lubridate)
library(dplyr)
library(tidyr)
library (sf)

Sys.setenv(TZ = "UTC")


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
