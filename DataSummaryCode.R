#obtaining data from motus towers

library(motus)
library(lubridate)
library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(tidyr)

Sys.setenv(TZ = "UTC") #sets the working environment to be in UTC to match motus

#tower at Western UWO - SSC (#11645)
#receiver CTT-7292E3C5E490

# ================================================================
# summarize station overview workflow
# ===============================================================
# year range, missing years, species #, detection #, tag #

# Define the folder path with the station data
folder_path <- "StationDownloads"

# List all CSV files in the folder
file_list <- list.files(path = folder_path, 
                        pattern = "\\.csv$", 
                        full.names = TRUE)

#apply function to the file list
all_detections <- map_df(file_list, function(file) {
  
# Extract Station ID from the filename
  filename_clean <- basename(file) #basename removes filepath from name
  #looks for "station-" then captures all digits until hits a number
  extracted_station_id <- str_extract(filename_clean, "(?<=station-)\\d+")
  
# read_csv (readr) to import and read each .csv file in the folder
  data_raw <- read_csv(file, show_col_types = FALSE)
  
# Format the date columns using lubridate functions
  data_cleaned <- data_raw %>%
    mutate(
      CurrentStationID = extracted_station_id,
      tsStart_dt = dmy_hms(tsStart, tz = "UTC"),
      tsEnd_dt   = dmy_hms(tsEnd, tz = "UTC"),
      year_start = year(tsStart_dt),
      year_end   = year(tsEnd_dt)
    )
  
  return(data_cleaned)
})
  
# overview summary 
station_overview <- all_detections %>%
group_by(stationID = CurrentStationID) %>%
  summarise(
    earliest_detection_raw = min(tsStart_dt, na.rm = TRUE),#date range
    latest_detection_raw   = max(tsEnd_dt, na.rm = TRUE),
    
    min_y = min(c(year_start, year_end), na.rm = TRUE), #year range
    max_y = max(c(year_start, year_end), na.rm = TRUE),
    
    years_present = list(unique(c(year_start, year_end))), #list years
    
    number_of_species = n_distinct(species),
    number_of_detections = n(),
    number_of_tags = n_distinct(tagDeployID),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    earliest_detection = format(earliest_detection_raw, 
                                "%d %b %Y %H:%M:%S UTC"),
    latest_detection   = format(latest_detection_raw, 
                                "%d %b %Y %H:%M:%S UTC"),
    
    full_sequence = list(min_y:max_y), #list years with data
    #find if there are missing data years by comparing two lists of years
    missing_years_vec = list(setdiff(full_sequence, years_present)),
    #puts the number of missing years, or 0 if there are none
    num_years_no_detection = length(missing_years_vec),
    years_missing_detections = if_else(
      num_years_no_detection > 0, 
      paste(missing_years_vec, collapse = ", "), 
      ""
    ),
    
    year_range = if_else(min_y == max_y, 
                         as.character(min_y), #if only one year prints that year
                         paste(min_y, 
                               max_y, 
                               sep = "-") 
                         ),
    number_of_years = max_y - min_y
  ) %>%
  ungroup() %>%
  select(
    stationID, 
    earliest_detection,
    latest_detection,
    year_range, 
    number_of_years, 
    num_years_no_detection,
    years_missing_detections,
    number_of_species, 
    number_of_detections, 
    number_of_tags
  )
print(station_overview)
write.csv(station_overview, "combined_station_overview.csv")

#======================================================================
# summary data by species by station
# =====================================================================
species_station_overview <- all_detections %>% 
  group_by(stationID = CurrentStationID, species) %>% 
  summarize (
    number_of_detections = n (),
    number_of_tags_detected = n_distinct(tagDeployID),
    .groups = "drop"
  ) %>% 
  arrange(species)

print (species_station_overview)

write.csv(species_station_overview, "species_station_overview.csv")

species_station_yearly <- all_detections %>%
  group_by(stationID = CurrentStationID, species, year = year_start) %>%
  summarise(
    number_of_detections = n(),
    number_of_tags_detected = n_distinct(tagDeployID),
    .groups = "drop"
  ) %>% 
  arrange(species, year)

print(species_station_yearly)

write.csv(species_station_yearly, "species_station_yearly.csv")

# ============================================================
# species by day by year by station
# ============================================================

detailed_data <- all_detections %>%
  mutate(
    date_only = as.Date(tsStart_dt), #add column of just date no time
    day_of_year = yday(tsStart_dt) # add column of Julian date, might be useful
  ) %>%
  group_by(stationID = CurrentStationID, 
           species, 
           year = year_start, 
           date_only, 
           day_of_year) %>%
  summarise(
    number_of_detections = n(),
    number_of_tags_detected = n_distinct(tagDeployID),
    .groups = "drop"
  ) %>%
  arrange(stationID, date_only, species) #order the detections

print(detailed_data)

write.csv(detailed_data, "species_station_daily.csv")

#===================================================================
# species summary
# =================================================================

species_overview <- all_detections %>% 
  group_by(species = species, 
           year = year_start) %>% 
  summarise (
    number_of_stations = n_distinct(CurrentStationID), 
    number_of_detections = n (), 
    number_of_tags_detected = n_distinct(tagDeployID), 
    .groups = "drop"
    ) %>% 
    arrange (species, year)

write.csv(species_overview, "species_overview.csv")

