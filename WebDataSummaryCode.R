#obtaining data from motus towers downloading from motus website

library(motus)
library(lubridate)
library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(tidyr)

Sys.setenv(TZ = "UTC") 

# ================================================================
# summarize station overview workflow
# ===============================================================
# year range, missing years, species #, detection #, tag #

data_active <- read.csv("./StationDownloads/Motus-data-project-1_detections_downloaded-2026-06-03.csv")
#at first just did active towers to prevent crashing
data_inactive <- read.csv ("./StationDownloads/Motus-data-multiple-stations_detections_downloaded-2026-06-09.csv")

#the inactive data had frequency as num, convert so it matches
data_inactive <- data_inactive %>%  
  mutate(tagFrequency = as.character(tagFrequency))

#join active and inactive stations together

data_raw <- bind_rows(data_active, data_inactive)

# Format the date columns using lubridate functions
data_cleaned <- data_raw %>%
  mutate(
    tsStart_dt = dmy_hms(tsStart, tz = "UTC"),
    tsEnd_dt   = dmy_hms(tsEnd, tz = "UTC"),
    year_start = year(tsStart_dt),
    year_end   = year(tsEnd_dt)
    )

#================================================================
# Station metadata json file
#================================================================

# Station metadata json file with GPS coordinates
stationsJSON <- jsonlite::read_json("./JSONFiles/Stations.json")

stationsJSON <- jsonlite::fromJSON("./JSONFiles/Stations.json", 
                                   simplifyDataFrame = TRUE)

readr::write_csv(stationsJSON$results, "StationsJSON.csv")

#getting coordinates from stationsJSON to make a lookup table for 
# previous tower coordinates
id_coords_lookup <- stationsJSON$results %>% 
  select(stationID, latitude, longitude) %>% 
  filter(!is.na(latitude) & !is.na(longitude)) %>% 
  distinct(stationID, .keep_all = TRUE)

# Add previous tower coordinates using the ID columns
summary_coords <- data_cleaned %>% 
  # convert to integer, from character
  mutate(previousStationID = as.integer(previousStationID)) %>% 
  
  # Join the lookup table to the previousStationID column
  left_join(
    id_coords_lookup, 
    by = c("previousStationID" = "stationID")
  ) %>% 
  
  # Rename the new columns oordinates
  rename(
    lat_previous = latitude,
    lon_previous = longitude
  )

#===================================================================
# data cleaning and binning by time periods for station overview
# ==================================================================

station_overview <- summary_coords %>%
#Classify migration periods and seasons based on tsStart_dt
  mutate(
#using month and day as a number (MMDD)
    md = as.numeric(format(tsStart_dt, "%m%d")),
    season = case_when(
      md >= 0415 & md < 0615  ~ "Spring Migration", #mid April to mid June
      md >= 0615 & md < 0815  ~ "Summer", #mid april to mid August
      md >= 0815 & md < 1115  ~ "Fall Migration", #mid Aug to mid Nov
      TRUE                     ~ "Winter" # mid-Nov to mid-Apr 
    )
  ) %>%

  # Group by station and season
  group_by(stationID, 
           previousStationID,
           season, 
           lat, 
           lon, 
           lat_previous, 
           lon_previous) %>%
  summarise(
    stationName = first(stationName), #adds station names without grouping
    previousStationName = first(previousStationName),
    earliest_detection_raw = min(tsStart_dt, na.rm = TRUE),
    latest_detection_raw = max(tsEnd_dt, na.rm = TRUE),
    
    min_y = min(c(year_start, year_end), na.rm = TRUE),
    max_y = max(c(year_start, year_end), na.rm = TRUE),
    
    years_present_list = list(unique(c(year_start, year_end))),
    
    number_of_species = n_distinct(species),
    number_of_detections = n(),
    number_of_tags = n_distinct(tagDeployID),
    
    .groups = "drop"
  ) %>%
  
  rowwise() %>%
  mutate(
    earliest_detection = format(earliest_detection_raw, "%d %b %Y %H:%M:%S UTC"),
    latest_detection   = format(latest_detection_raw, "%d %b %Y %H:%M:%S UTC"),
    
    year_range = if_else(min_y == max_y, 
                         as.character(min_y), 
                         paste(min_y, max_y, sep = "-")),
    number_of_years = max_y - min_y,
    
    # Find the missing operational years in the sequence range
    full_sequence = list(min_y:max_y),
    missing_years_vec = list(setdiff(unlist(full_sequence), 
                                     unlist(years_present_list))),
    
    # Save these explicitly as normal, flat data columns (No list nesting here!)
    num_years_no_detection = length(unlist(missing_years_vec)),
    years_missing_detections = if_else(
      num_years_no_detection > 0,
      paste(unlist(missing_years_vec), collapse = ", "),
      ""
    )
  ) %>%
  ungroup() %>%
  select(
    stationID,
    previousStationID,
    stationName, 
    previousStationName,
    season, 
    earliest_detection,
    latest_detection,
    year_range, 
    number_of_years, 
    num_years_no_detection,
    years_missing_detections,
    number_of_species, 
    number_of_detections, 
    number_of_tags,
    lat,
    lon,
    lat_previous,
    lon_previous
  )

write_csv(station_overview, "combined_station_overview.csv")

#======================================================================
# station network yearly by species summary
#=====================================================================
#comparing which towers are most often connected for each sp and year

station_pairs_summary <- summary_coords %>%

  #Classify Seasons 
  mutate(
    md = as.numeric(format(tsStart_dt, "%m%d")),
    season = case_when(
      md >= 0415 & md < 0615  ~ "Spring Migration",
      md >= 0615 & md < 0815  ~ "Summer",
      md >= 0815 & md < 1115  ~ "Fall Migration",
      TRUE                     ~ "Winter"
    )
  ) %>%
  
  # Group by station connections species, and season
  group_by(stationID,
           previousStationID,
           species, 
           season,
           lat,
           lon,
           lat_previous,
           lon_previous) %>%
  
  # Aggregate tracking and stopover metrics for the pair
  summarise(
    stationName = first(stationName),
    previousStationName = first(previousStationName),
    number_of_movements = n(),
    number_of_individual_tags = n_distinct(tagDeployID), 
    
    # Year range for this specific connection
    min_year = min(c(year_start, year_end)),
    max_year = max(c(year_start, year_end)),
    
    
    .groups = "drop"
  ) %>%
  
  # Create year range text column
  mutate(
    years_active = if_else(min_year == max_year, 
                           as.character(min_year), 
                           paste(min_year, max_year, sep = "-"))
  ) %>%
  
  # Order from where, to where, what species, and when
  arrange(number_of_movements)

station_pairs_summary

write.csv(station_pairs_summary, "station_pairs_summary.csv")


#===================================================================
# data organization for observation level
# ==================================================================

data_obs <- summary_coords %>%
  #Classify migration periods and seasons based on tsStart_dt
  mutate(
    #using month and day as a number (MMDD)
    md = as.numeric(format(tsStart_dt, "%m%d")),
    season = case_when(
      md >= 0415 & md < 0615  ~ "Spring Migration", #mid April to mid June
      md >= 0615 & md < 0815  ~ "Summer", #mid april to mid August
      md >= 0815 & md < 1115  ~ "Fall Migration", #mid Aug to mid Nov
      TRUE                     ~ "Winter" # mid-Nov to mid-Apr 
    ))

## adding columns for sunrise and sunset, and movement duration

library(suntools)

data_obs <- data_obs %>%
  mutate(
    movement_duration_hours = as.numeric(difftime(tsEnd_dt, 
                                                  tsStart_dt, 
                                                  units = "hours")),
    # Calculate sunrise and sunset in UTC
    sunrise_utc = sunriset(cbind(lon, lat), 
                           tsStart_dt, 
                           direction = "sunrise", 
                           POSIXct.out = TRUE)$time,
    sunset_utc  = sunriset(cbind(lon, lat), 
                           tsStart_dt, 
                           direction = "sunset", 
                           POSIXct.out = TRUE)$time,
    
    #Convert to  Eastern Time so they make  sense
    sunrise_local = with_tz(sunrise_utc, tzone = "America/Toronto"),
    sunset_local  = with_tz(sunset_utc, tzone = "America/Toronto")
  )
    
write.csv(data_obs, "observationSummary.csv")

#saving a file for species natural history to be uses as lookup

speciesList <- unique(data_obs$species)
write.csv (speciesList, "data_obs_speciesList.csv")


