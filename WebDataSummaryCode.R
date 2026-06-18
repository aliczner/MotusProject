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

data_raw <- bind_rows(data_active, data_inactive) %>%
  select(-tagFrequency) %>%  # some rows have it some do not
  distinct() 

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
  select(stationID, 
         latitude, 
         longitude) %>% 
  filter(!is.na(latitude) & !is.na(longitude)) %>% 
  # Convert stationID to integer
  mutate(stationID = as.integer(stationID)) %>% 
  # Group by ID and take the first coordinate 
  group_by(stationID) %>% 
  summarise(
    latitude = first(latitude),
    longitude = first(longitude),
    .groups = "drop"
  )

# Add previous station coordinates using the ID columns
summary_coords <- data_cleaned %>% 
  # Convert both ID columns to integers
  mutate(
    previousStationID = as.integer(previousStationID),
    stationID         = as.integer(stationID)
  ) %>% 
  
  # Join the lookup table to the previousStationID column
  left_join(
    id_coords_lookup, 
    by = c("previousStationID" = "stationID")
  ) %>% 
  
  # Rename the new columns coordinates
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


station_pairs <- read.csv("observationSummary.csv") #53515 obs
str(station_pairs)

#====================================================================
# station pair filtering column creation for later filtering
#===================================================================

#filter out same station pairs and station unknown

filtered_pairs <- station_pairs %>% 
  filter(
    previousStationID != stationID,
    previousStationName != "Unknown",
    previousStationName != "Unknown station"
  )
#55285 obs

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

str(filtered_df) #27680 obs

table(filtered_df$valid_distance)
#Valid_distance
#TRUE 11945
#FALSE_LESS 14376
#FALSE_MORE 1359

filtered_df$flight <- ifelse(filtered_df$movement_duration_hours < 12,
                             "flight",
                             "incidence") 

write.csv(filtered_df, "StationPairsFiltered.csv")

#==================================================
# summary data by individual
#=================================================

length(unique(filtered_df$species)) #124
length(unique(filtered_df$tagDeployID)) #4494

library(ggplot2)

# Calculate the number of detections per tag
tag_counts <- filtered_df %>% 
  group_by(tagDeployID) %>% 
  summarise(total_detections = n(), .groups = "drop")
# max = 473
# median = 3
# mean = 6.2

# Plot the histogram
ggplot(tag_counts, aes(x = total_detections)) +
  geom_histogram(binwidth = 5, fill = "purple4", color = "white") +
  labs(
    title = "Distribution of Detections per Tag",
    x = "Number of Detections",
    y = "Count of Tags"
  ) +
  theme_minimal()

#zoomed in histogram without such a long tail
ggplot(tag_counts, aes(x = total_detections)) +
  geom_histogram(binwidth = 5, 
                 boundary = 0, 
                 fill = "purple4", 
                 colour = "white") +
  coord_cartesian(xlim = c(1, 50)) + #zooms in
  labs(
    title = "Distribution of Detections per Tag",
    x = "Number of Detections",
    y = "Count of Tags"
  ) +
  theme_minimal()

#====================================================================
# checking for, and counting duplicate tsStart/End
#====================================================================

# adding in a buffer for potential processing lag between stations
buffer_seconds <- 3

duplicate_counts <- filtered_df %>%
#sort by tagID and date so that observations are chronological
  arrange(tagDeployID, tsStart_dt) %>%
  group_by(tagDeployID) %>%
  
  # calculating time difference within row and across rows (next step in flight)
  mutate(
    # within row time difference 
    within_row = as.numeric(difftime(tsEnd_dt, 
                                     tsStart_dt, 
                                     units = "secs")),
    # across row time difference behind
    cross_row_lag = as.numeric(difftime(lag(tsEnd_dt), 
                                    tsStart_dt, 
                                    units = "secs")),
    #across row time difference forward to ensure all pairs are flagged
    cross_row_lead = as.numeric(difftime(tsEnd_dt, 
                                         lead(tsStart_dt), 
                                         units = "secs"))
  ) %>%
  ungroup() %>%
  
  # Count the occurrences
  summarise(
    total_rows_in_dataset = n(),
    # Flag 1: Start and end are the same
    Flag_1_duration = sum(within_row >= 0 & within_row <= buffer_seconds,
                          na.rm = TRUE),
    # Flag 2: Overlaps with the previous row OR overlaps with the next row
    Flag_2_duration = sum(cross_row_lag > -buffer_seconds | cross_row_lead > -buffer_seconds,
                          na.rm = TRUE)
  )
duplicate_counts 
#  total_rows_in_dataset Flag_1_duration Flag_2_duration
#               27680              26           3548


# adding a column that classifies the detections by flags
# stationDuplicateFlag:
  # none means no issues
  # flag 1 the start and end time are the same
  # flag 2 multiple towers detected at once, differing start times between rows

filtered_df <- filtered_df %>%
  arrange(tagDeployID, tsStart_dt) %>%
  group_by(tagDeployID) %>%
  mutate(
    stationDuplicateFlag = case_when(
      # Label flag 1
      as.numeric(difftime(tsEnd_dt, 
                          tsStart_dt, 
                          units = "secs")) >= 0 & 
        as.numeric(difftime(tsEnd_dt, 
                            tsStart_dt, 
                            units = "secs")) <= buffer_seconds ~ "flag_1",
      # Label flag 2
      as.numeric(difftime(lag(tsEnd_dt), 
                          tsStart_dt, 
                          units = "secs")) > -buffer_seconds |
        as.numeric(difftime(tsEnd_dt, 
                            lead(tsStart_dt), 
                            units = "secs")) > -buffer_seconds ~ "flag_2",
      
      # Label everything else
      TRUE ~ "none"
    ),
      rolling_max_end = accumulate(as.POSIXct(tsEnd_dt),
                                   pmax),
      
      # If a row's start time is > rolling max end time of all previous rows, 
      # it is a clean break
      new_event_break = as.POSIXct(tsStart_dt) > lag(rolling_max_end),
      
      # Replace initial NA from the lag step with FALSE
      new_event_break = ifelse(is.na(new_event_break), 
                               FALSE,
                               new_event_break),
      
      # Summing the breaks to give a unique cluster ID
      eventClusterID = cumsum(new_event_break) + 1
    ) %>%
      ungroup() %>%
      
      select(-rolling_max_end, -new_event_break) #remove temp columns

head(filtered_df)

write.csv(filtered_df, "StationPairsFiltered.csv")


#=================================================================
# to fix flag_1s
#=================================================================

library(geosphere)

# Update flag_1 rows and calculate the midpoint
df_with_midpoints <- filtered_df %>%
  mutate(
    # Create the matrix inputs for the geosphere function
    p1 = cbind(lon, lat),
    p2 = cbind(lon_previous, lat_previous),
    
    # Calculate the true great-circle midpoint
    mid_coords = geosphere::midPoint(p1, p2),
    
    # Update current lat/lon and station name ONLY for flag_1 rows
    lon = if_else(stationDuplicateFlag == "flag_1", 
                  mid_coords[, 1], 
                  lon),
    lat = if_else(stationDuplicateFlag == "flag_1", 
                  mid_coords[, 2], 
                  lat),
    stationName = if_else(stationDuplicateFlag == "flag_1", 
                          paste0(previousStationName, " & ", stationName), 
                          stationName)
  ) %>% 
  select(-p1, -p2, -mid_coords)

# remove the flag_1 row from the path
final_cleaned_df <- df_with_midpoints %>%
  arrange(tagDeployID, tsStart_dt) %>%
  group_by(tagDeployID) %>%
  mutate(
    # Check if the row above was a flag_1 row 
    follows_flag1 = lag(stationDuplicateFlag == "flag_1",
                        default = FALSE),
    
    # If row follows a flag_1, add midpoint details
    previousStationName = case_when(
      follows_flag1 == TRUE ~ lag(stationName),
      TRUE                  ~ previousStationName
    ),
    previousStationID = case_when(
      follows_flag1 == TRUE ~ lag(stationID),
      TRUE                  ~ previousStationID
    ),
    lon_previous = case_when(
      follows_flag1 == TRUE ~ lag(lon),
      TRUE                  ~ lon_previous
    ),
    lat_previous = case_when(
      follows_flag1 == TRUE ~ lag(lat),
      TRUE                  ~ lat_previous
    ),
    
    #rename from flag_1 to flag_1_resolved
    stationDuplicateFlag = case_when(
      follows_flag1 == TRUE ~ "flag_1_resolved",
      TRUE                  ~ stationDuplicateFlag
    )
  ) %>%
  # drop the flag_1 rows 
  filter(stationDuplicateFlag != "flag_1") %>%
  
  # Recalculate distances with midpoints
  mutate(
    distance_km = case_when(
      follows_flag1 == TRUE ~ geosphere::distHaversine(cbind(lon_previous, 
                                                             lat_previous), 
                                                       cbind(lon, lat)) / 1000,
      TRUE                  ~ distance_km
    )
  ) %>%
  select(-follows_flag1) %>% # drops temporary columns
  ungroup()

#===========================================
# fix flag 2s
#===========================================

# all flag_2s have two events
# one event is a glitch (going back in time, the other is true)

# Keep only the forward-moving data
cleaned_df <- final_cleaned_df %>%
  filter(
    # If it's a flag_2, keep it if it moves forward/synchronously in time
    (stationDuplicateFlag == "flag_2" & tsStart_dt <= tsEnd_dt) |
      
      # Keep all other rows 
      (stationDuplicateFlag != "flag_2")
  )


write.csv(cleaned_df, "StationPairsFiltered.csv") 

#===============================================================
# adding additional flight information
#===============================================================

#=====================================================================
# creating flight paths
# ====================================================================

library(dplyr)
library(tidyr)
library(lubridate)

prepared_df <- filtered_df %>% 
  arrange(tagDeployID, tsStart_dt) %>%
  group_by(tagDeployID) %>%
  mutate( #calculate the lag time between steps in hours
    step_gap = as.numeric(difftime(tsStart_dt, 
                                   lag(tsEnd_dt), 
                                   units = "hours")),
    #order can get messy for stations within range ending up with negative
    #values. converting these to zero
    step_gap = ifelse(step_gap < 0, 0, step_gap),
    #classify as separate flight if lag time is > 12 hrs
    new_flight = is.na(step_gap) | (step_gap) > 12,
    #numbering the flights
    flight_ID = paste0("Flight_", cumsum(new_flight))
  ) %>%
  ungroup() %>%
  select(
    tagDeployID, flight_ID, season, species, 
    stationID, stationName, lon, lat, 
    previousStationID, previousStationName, lon_previous, lat_previous, 
    tsStart_dt, tsEnd_dt, movement_duration_hours, distance_km, 
    sunrise_local, sunset_local
  )

#now split into nested list
flights_list <- prepared_df %>% 
  # Split the data into a list by tagID
  split(.$tagDeployID) %>% 
  
  # loop through each tag and split by Flight ID
  lapply(function(individual_df) {
    split(individual_df, individual_df$flight_ID)
  })

str(flights_list, max.level = 2) #set to level 2 to see just tag and flight list
#======================================================================
# creating a map of the points
#======================================================================

library(sf)
library(mapview)

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


