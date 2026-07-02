#obtaining data from motus towers downloading from motus website

library(motus)
library(lubridate)
library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(tidyr)
library(purrr)
library(sf)
library(geosphere)

Sys.setenv(TZ = "UTC") 

# ================================================================
# summarize station overview workflow
# ===============================================================
# year range, missing years, species #, detection #, tag #

# downloaded data for all states/ON in GLWS but in batches
# the full region would not download
data_ON <- read.csv("./StationDownloads/Motus-data-region-CA-ON_detections_downloaded-2026-06-23.csv")
data_NY <- read.csv("./StationDownloads/Motus-data-region-US-NY_detections_downloaded-2026-06-23.csv")
data_PA <- read.csv("./StationDownloads/Motus-data-region-US-PA_detections_downloaded-2026-06-23.csv")
data_OH <- read.csv("./StationDownloads/Motus-data-region-US-OH_detections_downloaded-2026-06-23.csv")
data_IN <- read.csv("./StationDownloads/Motus-data-region-US-IN_detections_downloaded-2026-06-23.csv")
data_IL <- read.csv("./StationDownloads/Motus-data-region-US-IL_detections_downloaded-2026-06-23.csv")
data_WI <- read.csv("./StationDownloads/Motus-data-region-US-WI_detections_downloaded-2026-06-23.csv")

# join all the station data 
data_raw <- bind_rows(data_ON,
                      data_NY,
                      data_PA,
                      data_OH,
                      data_IN,
                      data_IL,
                      data_WI) %>%  # some rows have it some do not
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
# data cleaning and binning 
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
    
write_csv(data_obs, "observationSummary.csv")

#====================================================================
# station pair filtering column creation for later filtering
#===================================================================

station_pairs <- read.csv("observationSummary.csv") #142303 obs
str(station_pairs)

#filter out same station pairs and station unknown

filtered_pairs <- station_pairs %>% 
  filter(
    previousStationID != stationID,
    previousStationName != "Unknown",
    previousStationName != "Unknown station"
  )
str(filtered_pairs)
#94643 obs

#loading great lakes watershed polygon, contains subbasins for each lake
GLWatershed <- st_read("./greatlakes_subbasins/greatlakes_subbasins.shp")
#checking CRS
st_crs(GLWatershed) #ESPG 6269

#transform both the polygon and station pairs to lambert conformal conic

#need to run current and previous stations separately
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

str(filtered_df) #35865  obs

table(filtered_df$valid_distance)
#Valid_distance
#TRUE 13706
#FALSE_LESS 20377
#FALSE_MORE 1782

filtered_df$flight <- ifelse(filtered_df$movement_duration_hours < 12,
                             "flight",
                             "incidence") 

write_csv(filtered_df, "StationPairsFiltered.csv")
#35865 obs

#====================================================================
# Removing observations where tsStart is after tsEnd
#====================================================================
cleaned_backwards <- filtered_df %>%
  mutate(
    tsStart_dt = ymd_hms(tsStart_dt, tz = "UTC"),
    tsEnd_dt = ymd_hms(tsEnd_dt, tz = "UTC")
  ) %>%
  filter(tsEnd_dt >= tsStart_dt)

str(cleaned_backwards) #27309

#================================================================
# calculating centroids when detections overlap more stations
#================================================================
#will also look for false positives so midpoints are not calculated
# for impossible distances
#none of the detections occurred at the same time

# flag false positives 
final_cleaned_df <- cleaned_backwards %>%
  mutate(
    track_status = case_when(
      # flag as local overlapping detections with high speed and low distance
      travelDistanceKMs <= 30 & speedKmH > 120 ~ "LOCAL_OVERLAP",
      # flag as false positive for very high speeds >30 km
      speedKmH > 120                           ~ "FALSE_POSITIVE",
      # all other flights flagged as valid
      TRUE                                     ~ "VALID_FLIGHT"
    )
  ) %>%
# arrange chronologically
  arrange(tagDeployID, tsStart_dt) %>%
  group_by(tagDeployID) %>%
  
  # Calculate midpoints for local overlaps mathematically (keeps numeric data intact)
  mutate(
    lon = if_else(track_status == "LOCAL_OVERLAP" & !is.na(lon_previous), (lon + lon_previous) / 2, lon),
    lat = if_else(track_status == "LOCAL_OVERLAP" & !is.na(lat_previous), (lat + lat_previous) / 2, lat),
    stationName = if_else(
      track_status == "LOCAL_OVERLAP", 
      paste0(previousStationName, " & ", stationName), 
      stationName
    )
  ) %>%
  
  # arrange chronologically
  arrange(tagDeployID, tsStart_dt) %>%
  group_by(tagDeployID) %>%
  
  # calculate centroids using geosphere
  mutate(
    centroid = if_else(
      track_status == "LOCAL_OVERLAP" & !is.na(lon_previous),
      #centroid requires a matrix
      geosphere::midPoint(cbind(lon_previous, 
                                lat_previous), 
                          cbind(lon, lat)),
      cbind(lon, lat)
    ),
    lon = centroid[, 1],
    lat = centroid[, 2],
    
    stationName = if_else(
      track_status == "LOCAL_OVERLAP", 
      paste0(previousStationName, " & ", stationName), 
      stationName
    )
  ) %>%
  select(-centroid) %>% # remove temporary matrix column
  
  # Put centroid coordinates to the next sequential row
  mutate(
    follows_overlap = lag(track_status == "LOCAL_OVERLAP", 
                          default = FALSE),
    
    previousStationName = if_else(follows_overlap, 
                                  lag(stationName), 
                                  previousStationName),
    previousStationID   = if_else(follows_overlap, 
                                  lag(stationID), 
                                  previousStationID),
    lon_previous        = if_else(follows_overlap, 
                                  lag(lon), 
                                  lon_previous),
    lat_previous        = if_else(follows_overlap, 
                                  lag(lat), 
                                  lat_previous),
    
    # Update track_status directly
    track_status = if_else(follows_overlap & 
                             track_status == "VALID_FLIGHT", 
                           "LOCAL_OVERLAP_RESOLVED", 
                           track_status)
  ) %>%
  
  # Remove the duplicate row, unless it is the only one
  filter(track_status != "LOCAL_OVERLAP" | lead(follows_overlap, 
                                                default = FALSE) == FALSE) %>%
  
  # Recalculate distances for the centroids
  mutate(
    distance_km = if_else(
      follows_overlap,
      geosphere::distHaversine(cbind(lon_previous, 
                                     lat_previous), 
                               cbind(lon, lat)) / 1000,
      distance_km
    )
  ) %>%
  select(-follows_overlap) %>% 
  ungroup()

table(final_cleaned_df$track_status)

# =============================================
# adding additional flight information
#===============================================================

flightInfo_df <- final_cleaned_df %>%
  # sort by individual and in chronological order
  arrange(tagDeployID, tsStart_dt) %>%
  group_by(tagDeployID) %>%
  mutate(
      # Mark a TRUE every time a new flight starts (first row or > 12 hrs)
      new_flight = row_number() == 1 | movement_duration_hours > 12,
      
      # add flight number per individual
      flight_number = cumsum(new_flight),
      
      # Create the unique text identifier needed for splitting/nesting later
      flight_ID = paste0("Flight_", flight_number)
    ) %>%
  # Regroup by animal AND flight number for step and angle math
  group_by(tagDeployID, 
           flight_number) %>%
  mutate(
    # Speed of each step in a flight (km divided by travel duration hours)
    # Using if_else to prevent division-by-zero errors if duration is 0
    step_speed_kmh = if_else(travelDurationHours > 0, 
                             distance_km / travelDurationHours, 
                             0),
    
    # Turning angle of each step
    # calculate absolute headings for current step
    current_bearing = geosphere::bearing(cbind(lon_previous, lat_previous), 
                                         cbind(lon, lat)),
    # Look at the bearing of the next step within the same flight path
    next_bearing = lead(current_bearing),
    
    # Calculate the change in direction (-180 to +180 degrees)
    turning_angle = next_bearing - current_bearing,
    turning_angle = (turning_angle + 180) %% 360 - 180
  ) %>%
  ungroup() %>%
  
  mutate(
    # Convert character columns to formal date-time objects
    tsStart_posix   = ymd_hms(tsStart_dt),
    sunrise_posix   = ymd_hms(sunrise_local),
    sunset_posix    = ymd_hms(sunset_local),
    
    # Add a column for whether the flight was during daylight or at night
    diel_period = if_else(tsStart_posix >= sunrise_posix & 
                            tsStart_posix <= sunset_posix, 
                          "daylight", 
                          "night"),
    
    # Calculate  hour gaps from twilight transitions
    hours_from_sunset  = as.numeric(abs(difftime(tsStart_posix, 
                                                 sunset_posix, 
                                                 units = "hours"))),
    hours_from_sunrise = as.numeric(abs(difftime(tsStart_posix, 
                                                 sunrise_posix, 
                                                 units = "hours"))),
    
    # Add nearSun column based on the 1-hour window proximity
    nearSun = case_when(
      hours_from_sunset <= 1  ~ "sunset",
      hours_from_sunrise <= 1 ~ "sunrise",
      TRUE                    ~ "none"
    )
  ) %>%
  # Clean up temporary helper columns to keep the data tidy
  select(-current_bearing, 
         -next_bearing,
         -tsStart_posix, 
         -sunrise_posix, 
         -sunset_posix, 
         -hours_from_sunset,
         -hours_from_sunrise)

names(flightInfo_df)

str(flightInfo_df) #27181 obs

# adding the subbasin information to the dataframe
#needs to be done separately for current vs previous station

# current station
flightInfo_sf_current <- st_as_sf(
  flightInfo_df, 
  coords = c("lon", "lat"), 
  crs = st_crs(GLWatershed), 
  remove = FALSE
)

# joining to the data for current station
flightInfo_geo <- st_join(flightInfo_sf_current, 
                          GLWatershed, 
                          join = st_intersects) %>%
  rename(subbasin = merge) # rename current subbasin column

# previous station, the geometry from current needs to be removed
flightInfo_df_temp <- st_drop_geometry(flightInfo_geo)

# Re-convert to sf using the previous coordinates
flightInfo_sf_prev <- st_as_sf(
  flightInfo_df_temp,
  coords = c("lon_previous", "lat_previous"), 
  crs = st_crs(GLWatershed), 
  remove = FALSE
)

# joinging to the data for the previous station
flightInfo_geo <- st_join(flightInfo_sf_prev, 
                            GLWatershed, 
                            join = st_intersects) %>%
  rename(subbasin_previous = merge) %>% # rename previous subbasin column
  st_drop_geometry()# Drop the final geometry column

table(flightInfo_geo$subbasin, useNA= "ifany")
table(flightInfo_geo$subbasin_previous, useNA = "ifany")

write.csv(flightInfo_geo, "StationPairsFiltered.csv", row.names = FALSE)

#==============================================================
# adding animal information (tag site, sex, age)
#==============================================================

#the species names do not match, us the species metadata to pull species
#name by the species number

speciesMeta <- read.csv("./StationDownloads/GLWS_Species_Metadata.csv")

animalMeta <- read.csv("./StationDownloads/GLWS_Animal_Metadata.csv")

# Add the species name to animalMeta using a left join
animalMeta <- animalMeta %>%
  left_join(
    speciesMeta %>% 
      select(speciesID, speciesName = english) %>% 
      distinct(speciesID, .keep_all = TRUE),
    by = c("species" = "speciesID")
  )

#rename the columns so when it is joined it makes sense
#there are many instances with the same tagID even by species
#it looks like this happens when redeploying a tag that fell off, so will
#sort chronologically, and take the most recent record
animalMeta <- animalMeta %>%
  select(tagID, 
         name, 
         tag_year, 
         lat_tagSite, 
         lon_tagSite, 
         age, 
         sex) %>%
  arrange(tagID, 
          name, 
          tag_year) %>%  # Sorts chronologically by tag year
  group_by(tagID,
           name) %>%
  slice_tail(n = 1) %>% #takes the last row (the most recent)
  ungroup()

#adding in tag site info and sex and age info to the main df

animalInfo <- flightInfo_geo %>%
  left_join(
    animalMeta %>% 
      select(tagID, 
             lat_tagSite, 
             lon_tagSite, 
             age, 
             sex, 
             name),
    by = c("tagDeployID" = "tagID", 
           "species" = "name")
  )
  

write_csv(animalnfo, "StationPairsFiltered.csv", row.names = FALSE)
#==================================================
# summary data by individual
#=================================================

length(unique(flightInfo_df$species)) #140
length(unique(flightInfo_df$tagDeployID)) #5452

library(ggplot2)

# Calculate the number of detections per tag
tag_counts <- flightInfo_df %>% 
  group_by(tagDeployID) %>% 
  summarise(total_detections = n(), .groups = "drop")
# max = 38
# median = 3
# mean = 2.97

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


