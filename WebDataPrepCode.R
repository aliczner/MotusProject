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
data_MI <- read.csv("./StationDownloads/Motus-data-region-US-MI_detections_downloaded-2026-07-23.csv")
data_MN <- read.csv("./StationDownloads/Motus-data-region-US-MN_detections_downloaded-2026-07-23.csv")

# join all the station data 
data_raw <- bind_rows(data_ON,
                      data_NY,
                      data_PA,
                      data_OH,
                      data_IN,
                      data_IL,
                      data_WI,
                      data_MI,
                      data_MN) %>%  
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
    md = as.numeric(format(tsEnd_dt, "%m%d")),
    season = case_when(
      md >= 0301 & md < 0615  ~ "Spring Migration", #March to mid June
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

station_pairs <- read.csv("observationSummary.csv") #147043 obs
str(station_pairs)

#filter out same station pairs and station unknown

filtered_pairs <- station_pairs %>% 
  filter(
    previousStationID != stationID,
    previousStationName != "Unknown",
    previousStationName != "Unknown station"
  )
str(filtered_pairs)
#97441 obs

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



#==================================================================
# station pair filtering by criteria
#==================================================================

filtered_df <- filtered_pairs %>%
  filter(
    # Keep if either the current OR the previous station is in GLWS
    (current_in_GLWS | previous_in_GLWS))

str(filtered_df) #36837  obs

filtered_df$flight <- ifelse(filtered_df$movement_duration_hours < 12,
                             "flight",
                             "incidence") 

write_csv(filtered_df, "StationPairsFiltered.csv")

#====================================================================
# Removing observations where tsStart is after tsEnd
#====================================================================
cleaned_backwards <- filtered_df %>%
  mutate(
    tsStart_dt = ymd_hms(tsStart_dt, tz = "UTC"),
    tsEnd_dt = ymd_hms(tsEnd_dt, tz = "UTC")
  ) %>%
  filter(tsEnd_dt >= tsStart_dt)

str(cleaned_backwards) #28018


# =============================================================
# adding additional flight information
#===============================================================

flightInfo_df <- cleaned_backwards %>%
  arrange(tagDeployID, tsEnd_dt) %>%
  group_by(tagDeployID) %>%
  mutate(
    flight_type = if_else(movement_duration_hours > 12, 
                          "incidence", 
                          "flight"),
    
    # Flag start of new flight or put incidence
    is_new_flight = flight_type == "flight" & 
      (row_number() == 1 | lag(flight_type) == "incidence"),
    
    # Track flight sequence numbers
    flight_number = cumsum(is_new_flight),
    
    flight_ID = if_else(
      flight_type == "incidence",
      "incidence",
      paste0("Flight_", flight_number)
    )
  ) %>%
  group_by(tagDeployID, flight_ID) %>%
  mutate(
    step_speed_kmh = if_else(travelDurationHours > 0, 
                             distance_km / travelDurationHours, 
                             0),
    current_bearing = bearing(cbind(lon_previous, lat_previous),
                              cbind(lon, lat)),
    next_bearing = lead(current_bearing),
    turning_angle = next_bearing - current_bearing,
    turning_angle = (turning_angle + 180) %% 360 - 180
  ) %>%
  ungroup() %>%
  mutate(
    # Convert date-time strings
    tsStart_posix = ymd_hms(tsStart_dt),
    tsEnd_posix   = ymd_hms(tsEnd_dt),
    sunrise_posix = ymd_hms(sunrise_local),
    sunset_posix  = ymd_hms(sunset_local),
    
    # Classify diel period
    diel_period = case_when(
      tsStart_posix >= sunrise_posix & 
        tsEnd_posix <= sunset_posix ~ "daylight",
      tsEnd_posix < sunrise_posix | tsStart_posix > sunset_posix ~ "night",
      TRUE ~ "mixed"
    ),
    
    # Calculate hour gaps
    hours_from_sunset  = as.numeric(abs(difftime(tsEnd_posix, 
                                                 sunset_posix, 
                                                 units = "hours"))),
    hours_from_sunrise = as.numeric(abs(difftime(tsEnd_posix, 
                                                 sunrise_posix, 
                                                 units = "hours"))),
    
    # Near sun classification
    nearSun = case_when(
      hours_from_sunset >= 2  ~ "sunset",
      hours_from_sunrise <= 1 ~ "sunrise",
      TRUE                    ~ "none"
    )
  ) %>%
  # Clean up temporary calculations
  select(
    -current_bearing, 
    -next_bearing,
    -tsEnd_posix,
    -tsStart_posix, 
    -sunrise_posix, 
    -sunset_posix, 
    -hours_from_sunset,
    -hours_from_sunrise
  )
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
      select(speciesID, 
             speciesName = english) %>% 
      distinct(speciesID, 
               .keep_all = TRUE),
    by = c("species" = "speciesID")
  )


#rename the columns so when it is joined it makes sense
#there are many instances with the same tagID even by species
#it looks like this happens when redeploying a tag that fell off, so will
#sort chronologically, and take the most recent record

#first need to make a column for tagging year
# but the dtStart column has weird formatting 
animalMeta <- animalMeta %>%
  mutate(
    tag_year = as.numeric(str_extract(dtStart, 
                                      "\\b\\d{4}\\b")),
    species_clean = str_trim(str_remove(speciesName, 
                                        "\\s*\\(.*\\)"))
  ) %>% 
  rename(
    lat_tagSite = latitude,
    lon_tagSite = longitude
  )


animalMeta <- animalMeta %>%
  select(id, 
         species_clean, 
         tag_year, 
         lat_tagSite, 
         lon_tagSite, 
         age, 
         sex) %>%
  group_by(id,
           tag_year,
           species_clean) %>%
  ungroup()

#adding in tag site info and sex and age info to the main df

animalInfo <- flightInfo_geo %>%
  left_join(
    animalMeta %>% 
      select(id, 
             lat_tagSite, 
             lon_tagSite, 
             age, 
             sex, 
             species_clean,
             tag_year),
    by = c("tagDeployID" = "id", 
           "year_start" = "tag_year",
           "species" = "species_clean")
  )

write.csv(animalInfo, "StationPairsFiltered.csv", row.names = FALSE)

#==================================================================
# adding in animal type
#==================================================================

bats <- c(
  "Silver-haired Bat", "Eastern Red Bat", "Northern Hoary Bat", 
  "Little Brown Myotis", "Indiana Myotis", "Big Brown Bat", 
  "American Hoary Bat", "Northern Myotis", "Eastern Pipistrelle", 
  "Cave Myotis", "Eastern Small-footed Myotis", "Southern Yellow Bat"
)

shorebirds <- c(
  "Red Knot", "Ruddy Turnstone", "Black-bellied Plover", "Semipalmated Plover", 
  "Sanderling", "White-rumped Sandpiper", "Semipalmated Sandpiper", "Dunlin", 
  "Dunlin (hudsonia)", "Hudsonian Whimbrel", "American Woodcock", "Least Sandpiper", 
  "Lesser Yellowlegs", "Pectoral Sandpiper", "American Oystercatcher", 
  "Western Sandpiper", "Short-billed Dowitcher", "Piping Plover", 
  "Short-billed Dowitcher (griseus)", "Solitary Sandpiper", "Greater Yellowlegs"
)

raptors <- c(
  "Merlin", "Northern Saw-whet Owl", "Long-eared Owl", "Sharp-shinned Hawk", 
  "Peregrine Falcon", "Sharp-shinned Hawk (Madrean)", "American Kestrel", 
  "American Kestrel (Northern)"
)

water_birds <- c(
  "Black-crowned Night Heron", "Canada Goose", "Wood Duck", 
  "Least Bittern", "Black Tern", "Common Tern"
)

rails <- c(
  "Clapper Rail", "Virginia Rail", "Common Gallinule", "Sora"
)

insects <- c(
  "Common Green Darner", "Monarch"
)

groundbirds_doves <- c(
  "Sharp-tailed Grouse", "Rock Pigeon (Feral Pigeon)", 
  "Mourning Dove", "White-winged Dove"
)

aerial_insectivores <- c(
  "Common Nighthawk", "Eastern Whip-poor-will", "Chimney Swift", "Common Poorwill"
)

woodpeckers <- c(
  "Yellow-bellied Sapsucker", "Lewis's Woodpecker", "Northern Flicker"
)

cuckoos <- c(
  "Black-billed Cuckoo"
)

# add column for animalType
animalInfo <- animalInfo %>%
  mutate(
    # Detailed category column
    AnimalType = case_when(
      species %in% bats ~ "Bat",
      species %in% shorebirds ~ "Shorebird",
      species %in% raptors ~ "Raptor",
      species %in% water_birds ~ "Water Bird",
      species %in% rails ~ "Rail",
      species %in% insects ~ "Insect",
      species %in% groundbirds_doves ~ "Groundbird & Dove",
      species %in% aerial_insectivores ~ "Aerial Insectivore",
      species %in% woodpeckers ~ "Woodpecker",
      species %in% cuckoos ~ "Cuckoo",
      TRUE ~ "Songbird"
    ),
    
    # Broad category column (Bird, Bat, or Insect)
    Animal = case_when(
      species %in% bats ~ "Bat",
      species %in% insects ~ "Insect",
      TRUE ~ "Bird" 
    )
  )


#===================================================================
#histograms of takeoff/landing
#===================================================================

#spring night flights departure
spring_night_flights_start <- animalInfo %>%
  filter(season == "Spring Migration" & 
           flight_type == "flight" &
           diel_period == "night" &
           Animal == "Bird") %>% 
  group_by(tagDeployID, flight_number) %>% 
  slice_min(tsStart_dt,
            n=1,
            with_ties = FALSE) %>%  #makes sure only 1 row is selected
  mutate(
    sunset_dt = ymd_hms(sunset_local),
    
    # If the flight starts in the morning hours (before noon), 
    # the relevant sunset was on the previous day, this fixes that
    effective_sunset = if_else(
      hour(tsStart_dt) < 12, 
      sunset_dt - days(1), 
      sunset_dt
    ),
    
    # Calculate hours since the actual previous sunset
    hours_since_sunset = as.numeric(difftime(tsStart_dt, 
                                             effective_sunset, 
                                             units = "hours"))
  )

# make the histogram
ggplot(spring_night_flights_start, aes(x = hours_since_sunset)) +
  geom_histogram(binwidth = 0.5, fill = "#202C59", color = "white") +
  labs(
    title = "Spring Night Flight: Start",
    x = "Hours Since Sunset",
    y = "Number of Flights"
  ) +
  theme_minimal()

# spring mixed flights departure
spring_mixed_flights_start <- animalInfo %>%
  filter(season == "Spring Migration" & 
           flight_type == "flight" &
           diel_period == "mixed" &
           Animal == "Bird") %>% 
  group_by(tagDeployID, flight_number) %>% 
  slice_min(tsStart_dt,
            n=1,
            with_ties = FALSE) %>%  
  mutate(
    sunset_dt = ymd_hms(sunset_local),
    
    effective_sunset = if_else(
      hour(tsStart_dt) < 12, 
      sunset_dt - days(1), 
      sunset_dt
    ),
    
    hours_since_sunset = as.numeric(difftime(tsStart_dt, 
                                             effective_sunset, 
                                             units = "hours"))
  )

ggplot(spring_mixed_flights_start, aes(x = hours_since_sunset)) +
  geom_histogram(binwidth = 0.5, fill = "#202C59", color = "white") +
  labs(
    title = "Spring Mixed (Day/Night) Flights: Start",
    x = "Hours Since Sunset",
    y = "Number of Flights"
  ) +
  theme_minimal()


#spring night flight end/landing
spring_night_end <- animalInfo %>%
  filter(season == "Spring Migration" & 
           flight_type == "flight" &
           diel_period %in% c("night") &
           Animal == "Bird") %>%  
  group_by(tagDeployID, flight_number) %>% 
  slice_max(tsEnd_dt,
            n=1,
            with_ties = FALSE) %>% 
  mutate(
    sunrise_dt = ymd_hms(sunrise_local),
    
    effective_sunrise = if_else(
      hour(tsEnd_dt) > 12, 
      sunrise_dt + days(1), 
      sunrise_dt
    ),
    
    hours_relative_to_sunrise = as.numeric(difftime(tsEnd_dt, 
                                                    sunrise_dt, 
                                                    units = "hours"))
  ) %>% 
  filter(hours_relative_to_sunrise <= 5)

ggplot(spring_night_end, aes(x = hours_relative_to_sunrise)) +
  geom_histogram(binwidth = 0.5, fill = "#EF3054", color = "white") +
  labs(
    title = "Spring Night Flights: End Time",
    x = "Hours Relative to Sunrise (0 = Sunrise)",
    y = "Number of Flights"
  ) +
  theme_minimal()

#spring mixed flight end/landing
spring_mixed_end <- animalInfo %>%
  filter(season == "Spring Migration" & 
           flight_type == "flight" &
           diel_period %in% c("mixed") &
           Animal == "Bird") %>%  
  group_by(tagDeployID, flight_number) %>% 
  slice_max(tsEnd_dt,
            n=1,
            with_ties = FALSE) %>% 
  mutate(
    sunrise_dt = ymd_hms(sunrise_local),
    
    effective_sunrise = if_else(
      hour(tsEnd_dt) > 12, 
      sunrise_dt + days(1), 
      sunrise_dt
    ),
    
    hours_relative_to_sunrise = as.numeric(difftime(tsEnd_dt, 
                                                    sunrise_dt, 
                                                    units = "hours"))
  )

ggplot(spring_mixed_end, aes(x = hours_relative_to_sunrise)) +
  geom_histogram(binwidth = 0.5, fill = "#EF3054", color = "white") +
  labs(
    title = "Spring Mixed Flights: End Time",
    x = "Hours Relative to Sunrise (0 = Sunrise)",
    y = "Number of Flights"
  ) +
  theme_minimal()

#=======================================================
# adding column for migratory timing
#=======================================================
#remove subspecies, change to just species level
animalInfo <- animalInfo %>%
  mutate(species = sub(" \\(.*\\)", "", species))
#the space matches the space before the subspecies
#double backslashes are needed to tell R to look for ( )
#.* match any character, so any words in the bracket match

#used Winger et al 2019, Ralph 1981, Dufour et al 2026, Cornell BOW
nocturnal <- c (
  "White-throated Sparrow","Dark-eyed Junco", "Song Sparrow", "Swamp Sparrow",
  "Ovenbird", "Hermit Thrush", "Tennessee Warbler", "Fox Sparrow", 
  "American Tree Sparrow", "Magnolia Warbler", "White-crowned Sparrow",
  "Northern Waterthrush", "Yellow-rumped Warbler", "American Redstart",
  "Blackpoll Warbler", "Veery", "Palm Warbler", "Chestnut-sided Warbler",
  "Savannah Sparrow","Black-throated Green Warbler", "Baltimore Oriole",
  "Red-eyed Vireo", "Connecticut Warbler", "Mourning Warbler", 
  "Canada Warbler","Gray Catbird", "Wood Thrush", "Swainson's Thrush",
  "Gray-cheeked Thrush", "Black-throated Blue Warbler", 
  "Northern Waterthrush", "American Redstart", "Scarlet Tanager",
  "Nashville Warbler", "Lincoln's Sparrow", "Bay-breasted Warbler", 
  "Swainson's Warbler", "Painted Bunting", "Prothonotary Warbler",
  "Kirtland's Warbler", "American Woodcock", "Baird's Sparrow", 
  "Bank Swallow", "Bicknell's Thrush", "Black-bellied Plover",
  "Black-billed Cuckoo", "Black-crowned Night Heron", "Black Tern",
  "Bobolink", "Brown Thrasher", "Cerulean Warbler", "Clapper Rail", 
  "Common Gallinule", "Common Poorwill", "Common Tern", "Eastern Towhee",
  "Eastern Whip-poor-will", "Golden-winged Warbler", "Grasshopper Sparrow",
  "Henslow's Sparrow", "Kentucky Warbler", "Least Bittern", 
  "Least Sandpiper",  "Long-eared Owl", "Northern Saw-whet Owl",
  "Northern Yellow Warbler", "Pectoral Sandpiper", "Piping Plover", 
  "Prairie Warbler","Snowy Plover", "Solitary Sandpiper", "Sora", 
  "Vesper Sparrow", "Virginia Rail", "Willow Flycatcher", "Wood Duck",
  "Worm-eating Warbler","Yellow-bellied Sapsucker"
  
  
)

diurnal <- c(
  "American Kestrel", "American Oystercatcher", "American Pipit",
  "American Robin", "Barn Swallow", "Blue Jay", "Brown-headed Cowbird",
  "Chestnut-collared Longspur", "Chimney Swift", "Cliff Swallow", 
  "Evening Grosbeak", "Horned Lark", "Lewis's Woodpecker",  
  "Loggerhead Shrike", "Merlin", "Mourning Dove", "Peregrine Falcon", 
  "Purple Finch", "Purple Martin","Rusty Blackbird", "Sprague's Pipit",
  "Tree Swallow", "White-winged Dove" 
)

nonmigratory <- c(
  "Black-capped Chickadee", "Puerto Rican Oriole", "Rock Pigeon",
  "Sharp-tailed Grouse" 
)

mixed <- c (
  "Canada Goose", "Cedar Waxwing", "Common Nighthawk", "Dunlin", 
  "Greater Yellowlegs", "Hudsonian Whimbrel", "Lesser Yellowlegs",
  "Northern Flicker", "Pine Siskin", "Red Knot", "Ruddy Turnstone",
  "Sanderling", "Semipalmated Plover", "Semipalmated Sandpiper",
  "Sharp-shinned Hawk", "Snow Bunting","Short-billed Dowitcher",
  "Western Sandpiper","White-rumped Sandpiper" 
  
)

# add column for migration timing
animalInfo <- animalInfo %>%
  mutate(
    # Detailed category column
    MigrateTime = case_when(
      species %in% nocturnal ~ "nocturnal",
      species %in% diurnal ~ "diurnal",
      species %in% mixed ~ "mixed",
      species %in% nonmigratory ~ "nonmigratory",
      TRUE ~ "unclassified"
    ))
          
write.csv(animalInfo, "StationPairsFiltered.csv", row.names = FALSE)
