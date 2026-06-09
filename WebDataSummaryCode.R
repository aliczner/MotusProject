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
data_inactive <- data_inactive |> 
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

station_overview <- data_cleaned %>%
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
  #Binning the stopover duration
  mutate(
    stopover_bin = case_when(
      stopoverDurationHours == 0                  ~ "stopover_0",
      stopoverDurationHours > 0  & stopoverDurationHours < 1   ~ "stopover_under_1",
      stopoverDurationHours >= 1  & stopoverDurationHours <= 12  ~ "stopover_1_12",
      stopoverDurationHours > 12 & stopoverDurationHours <= 24  ~ "stopover_12_24",
      stopoverDurationHours > 24 & stopoverDurationHours <= 48  ~ "stopover_24_48",
      stopoverDurationHours > 48 & stopoverDurationHours <= 72  ~ "stopover_48_72",
      stopoverDurationHours > 72 & stopoverDurationHours <= 96  ~ "stopover_72_96",
      stopoverDurationHours > 96 & stopoverDurationHours <= 168 ~ "stopover_96_168",
      stopoverDurationHours > 168                 ~ "stopover_over_168",
      TRUE                                        ~ "stopover_unknown"
    )
  ) %>%
  # Group by station AND season
  group_by(stationName, season) %>%
  summarise(
    earliest_detection_raw = min(tsStart_dt, na.rm = TRUE),
    latest_detection_raw   = max(tsEnd_dt, na.rm = TRUE),
    
    min_y = min(c(year_start, year_end), na.rm = TRUE),
    max_y = max(c(year_start, year_end), na.rm = TRUE),
    
    years_present_list = list(unique(c(year_start, year_end))),
    
    number_of_species = n_distinct(species),
    number_of_detections = n(),
    number_of_tags = n_distinct(tagDeployID),
    
    #count the number of detections in each stopover bins and make summaries
    stopover_0       = sum(stopover_bin == "stopover_0", na.rm = TRUE),
    stopover_under_1  = sum(stopover_bin == "stopover_under_1", na.rm = TRUE),
    stopover_1_12     = sum(stopover_bin == "stopover_1_12", na.rm = TRUE),
    stopover_12_24    = sum(stopover_bin == "stopover_12_24", na.rm = TRUE),
    stopover_24_48   = sum(stopover_bin == "stopover_24_48", na.rm = TRUE),
    stopover_48_72   = sum(stopover_bin == "stopover_48_72", na.rm = TRUE),
    stopover_72_96    = sum(stopover_bin == "stopover_72_96", na.rm = TRUE),
    stopover_96_168  = sum(stopover_bin == "stopover_96_168", na.rm = TRUE),
    stopover_over_168 = sum(stopover_bin == "stopover_over_168", na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  
  rowwise() %>%
  mutate(
    earliest_detection = format(earliest_detection_raw, "%d %b %Y %H:%M:%S UTC"),
    latest_detection   = format(latest_detection_raw, "%d %b %Y %H:%M:%S UTC"),,
    
    year_range = if_else(min_y == max_y, 
                         as.character(min_y), 
                         paste(min_y, max_y, sep = "-")),
    number_of_years = max_y - min_y,
    
    # Find the missing operational years in the sequence range
    full_sequence = list(min_y:max_y),
    missing_years_vec = list(setdiff(unlist(full_sequence), unlist(years_present_list))),
    
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
    stationName, 
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
    stopover_0,
    stopover_under_1,
    stopover_1_12,
    stopover_12_24,
    stopover_24_48,
    stopover_48_72,
    stopover_72_96,
    stopover_96_168,
    stopover_over_168
  )

write_csv(station_overview, "combined_station_overview.csv")

#======================================================================
# station network yearly by species summary
#=====================================================================
#comparing which towers are most often connected for each sp and year

station_pairs_summary <- data_cleaned %>%

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
  
  # Bin the stopover duration at the current station
  mutate(
    stopover_bin = case_when(
      stopoverDurationHours == 0                  ~ "stopover_0",
      stopoverDurationHours > 0  & stopoverDurationHours < 1   ~ "stopover_under_1",
      stopoverDurationHours >= 1  & stopoverDurationHours <= 12  ~ "stopover_1_12",
      stopoverDurationHours > 12 & stopoverDurationHours <= 24  ~ "stopover_12_24",
      stopoverDurationHours > 24 & stopoverDurationHours <= 48  ~ "stopover_24_48",
      stopoverDurationHours > 48 & stopoverDurationHours <= 72  ~ "stopover_48_72",
      stopoverDurationHours > 72 & stopoverDurationHours <= 96  ~ "stopover_72_96",
      stopoverDurationHours > 96 & stopoverDurationHours <= 168 ~ "stopover_96_168",
      stopoverDurationHours > 168                 ~ "stopover_over_168",
      TRUE                                        ~ "stopover_unknown"
    )
  ) %>%
  
  # Group by station connections species, and season
  group_by(previousStationName, stationName, species, season) %>%
  
  # Aggregate tracking and stopover metrics for this specific link
  summarise(
    number_of_movements = n(),
    number_of_individual_tags = n_distinct(tagDeployID), 
    
    # Year range for this specific connection
    min_year = min(c(year_start, year_end)),
    max_year = max(c(year_start, year_end)),
    
    # Stopover duration breakdowns upon arrival at the current station
    stopover_0        = sum(stopover_bin == "stopover_0"),
    stopover_under_1  = sum(stopover_bin == "stopover_under_1"),
    stopover_1_12    = sum(stopover_bin == "stopover_1_12"),
    stopover_12_24    = sum(stopover_bin == "stopover_12_24"),
    stopover_24_48   = sum(stopover_bin == "stopover_24_48"),
    stopover_48_72   = sum(stopover_bin == "stopover_48_72"),
    stopover_72_96    = sum(stopover_bin == "stopover_72_96"),
    stopover_96_168   = sum(stopover_bin == "stopover_96_168"),
    stopover_over_168 = sum(stopover_bin == "stopover_over_168"),
    
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


#================================================================
# Station metadata json file
#================================================================


stationsJSON <- jsonlite::read_json("./JSONFiles/Stations.json")

stationsJSON <- jsonlite::fromJSON("./JSONFiles/Stations.json", 
                                   simplifyDataFrame = TRUE)

readr::write_csv(stationsJSON$results, "StationsJSON.csv")

