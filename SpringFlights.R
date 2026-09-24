# code for investigating spring migration phenology

library(dplyr)
library(lubridate)
library(ggplot2)
library(ggridges)
library(sf)
library(sfheaders)
library(spatstat.geom)
library(spatstat.explore)
library(terra)


fullData <- read.csv ("StationPairsFiltered.csv")

springData <- fullData %>% 
  filter(
    season == "Spring Migration",
  )
nrow(springData)

springTable <- springData %>%  
  mutate(
    tsEnd_dt = as_datetime(tsEnd_dt, tz = "GMT"),
    # merging georgian bay with lake huron so they are not separate. 
    subbasin = if_else(subbasin == "geo_bay", "lk_huron", subbasin)
  )
#=========================================================
# creating flight paths
#=========================================================

library(sfheaders)

flight_steps <- springTable %>%
  filter(flight_type != "incidence") %>%
  group_by(tagDeployID, flight_ID) %>%
  arrange(tsEnd_dt, .by_group = TRUE) %>%
  ungroup()

geoms <- lapply(seq_len(nrow(flight_steps)), function(i) {
  st_linestring(matrix(
    c(flight_steps$lon_previous[i], flight_steps$lon[i],
      flight_steps$lat_previous[i], flight_steps$lat[i]),
    ncol = 2
  ))
})
flight_lines <- st_sf(flight_steps, 
                      geometry = st_sfc(geoms, 
                                        crs = 4326))

#===========================================================
# creating a figure of search effort
# ======================================================
library(dplyr)
library(tidyr)
library(ggplot2)

flight_lines %>%
  select(lon, lon_previous, lon_tagSite) %>%
  pivot_longer(
    cols = c(lon, lon_previous, lon_tagSite),
    names_to = "longitude_type",
    values_to = "longitude"
  ) %>%
  mutate(
    longitude_type = case_match(
      longitude_type,
      "lon" ~ "Current Station",
      "lon_previous" ~ "Previous Station",
      "lon_tagSite" ~ "Tag Site"
    )
  ) %>%
  filter(!is.na(longitude) & longitude < 0) %>%
  ggplot(aes(x = longitude, 
             color = longitude_type, 
             fill = longitude_type)) +
  geom_density(alpha = 0.3, 
               linewidth = 0.9) +
  facet_wrap(~ longitude_type, 
             ncol = 1, 
             scales = "free_y") +
  scale_color_manual(values = c(
    "Current Station" = "#2b5c8f", 
    "Previous Station" = "#2d6a4f", 
    "Tag Site" = "#d95f02"
  )) +
  scale_fill_manual(values = c(
    "Current Station" = "#2b5c8f", 
    "Previous Station" = "#2d6a4f", 
    "Tag Site" = "#d95f02"
  )) +
  labs(
    title = "Longitude Distributions",
    x = "Longitude (°W)",
    y = "Density of Flights"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    strip.text = element_text(size = 11)
  )

#================================================
# checking flight timing
#=================================================

flight_time_summary <- flight_lines %>%
  group_by(MigrateTime, species, diel_period) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(MigrateTime, species) %>%
  mutate(
    total_flights = sum(n),
    percentage = (n / total_flights) * 100,
    species_label = paste0(species, " (", total_flights, ")")
  ) %>%
  ungroup()

library(ggplot2)

# diurnal barplot
plot_diurnal <- flight_time_summary %>%
  filter(MigrateTime == "diurnal") %>%
  ggplot(aes(x = reorder(species_label, percentage), 
             y = percentage, 
             fill = diel_period)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_manual (
    values = c (
      "daylight" = "#ECA72C",
      "night" = "#31263E"
    )) +
  labs(
    title = "Diurnal",
    x = "",
    y = "Percentage of Flights (%)",
    fill = "Diel Period"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# Mixed  Plot
plot_mixed <- flight_time_summary %>%
  filter(MigrateTime == "mixed") %>%
  ggplot(aes(x = reorder(species_label, percentage), 
             y = percentage, 
             fill = diel_period)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_manual (
    values = c (
      "daylight" = "#ECA72C",
      "night" = "#31263E"
    )) +
  labs(
    title = "Mixed",
    x = "",
    y = "Percentage of Flights (%)",
    fill = "Diel Period"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# nocturnal
plot_nocturnal <- flight_time_summary %>%
  filter(MigrateTime == "nocturnal") %>%
  ggplot(aes(x = reorder(species_label, percentage), 
             y = percentage, 
             fill = diel_period)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_manual (
    values = c (
      "daylight" = "#ECA72C",
      "night" = "#31263E"
    )) +
  labs(
    title = "Nocturnal",
    x = "",
    y = "Percentage of Flights (%)",
    fill = "Diel Period"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

#==============================================================
# partitioning night flights for nocturnal migrants
#==============================================================
# If the flight starts in the morning  its true starting sunset 
# was actually on the previous calendar day (- 1 day).
# Otherwise, it belongs to the sunset of the current calendar day.

### nocturnal take off
noc.takeoff <- flight_lines %>%  
  filter(MigrateTime == "nocturnal",
         Animal == "Bird") %>%  
  mutate(
    start_time = ymd_hms(tsStart_dt),
    sunset_time_dt = ymd_hms(sunset_utc_previous),
    true_sunset = if_else(
      hour(start_time) < 12, 
      sunset_time_dt - days(1), 
      sunset_time_dt
    ),
    
    # Calculate continuous positive hours
    hours_since_sunset = as.numeric(difftime(start_time, 
                                             true_sunset, 
                                             units = "hours"))
  ) %>%
  # Filter for flights starting between 0 and 2 hours after sunset
  filter(
    between(hours_since_sunset, -0.5, 2)
  )

### nocturnal landing

noc.landing <- flight_lines %>%  
  filter(
    MigrateTime == "nocturnal",
    Animal == "Bird") %>%  
  mutate(
    end_time = ymd_hms(tsEnd_dt),
    sunrise_time_dt = ymd_hms(sunrise_utc_previous),
    
    sunrise_raw_diff = as.numeric(difftime(end_time, 
                                           sunrise_time_dt, 
                                           units = "hours")),
    true_sunrise = case_when(
      sunrise_raw_diff > 12  ~ sunrise_time_dt + days(1),
      sunrise_raw_diff < -12 ~ sunrise_time_dt - days(1),
      TRUE                   ~ sunrise_time_dt
    ),
  
  # Calculate hours relative to sunrise (negative means before sunrise) 
  hours_relative_to_sunrise = as.numeric(difftime(end_time, 
                                                  true_sunrise, 
                                                  units = "hours"))
) %>%

  filter(
    between(hours_relative_to_sunrise, -4, 0.5)
  )



#==================================================================
# Multispecies Line Kernel Density Estimation Overlay All Birds
#===================================================================
library(terra)
library(sf)
library(dplyr)
library(purrr)
library(mapview)
library(stringr)
library(spatstat.geom)
library(spatstat.explore)

# create a blank template of 2.5 km for the region
regionTemplate <- rast(ext(flight_lines.pj), 
                       resolution = 2500, 
                       crs = st_crs(flight_lines.pj)$wkt)


#creating flight lines density
flight_lines.pj <- st_transform(flight_lines, crs = 3978)

Bird_lines.pj <- flight_lines.pj %>%
  filter (
    Animal == "Bird"
  )

st_write(Bird_lines.pj, 
         "Bird_lines.pj.gpkg", 
         delete_layer = TRUE)

# get the list of species to loop through
species_list <- unique(Bird_lines.pj$species)

# loop through the species to create the rasters
species_rasters <- lapply(species_list, 
                          function(sp_name) {
  
  df_sp <- Bird_lines.pj %>% 
    filter(species == sp_name)
  
  flightsRast <- rasterize(vect(df_sp), 
                           regionTemplate, 
                           field = 1, 
                           fun = "sum", 
                           background = 0)
  
  #this is for the smoothing window
  weightMatrix <- focalMat(flightsRast, 
                           d = 5000, 
                           type = "Gauss") 
  kdeSurface <- focal(flightsRast, 
                      w = weightMatrix, 
                      fun = sum, 
                      na.rm = TRUE)
  
  return(kdeSurface)
})

names(species_rasters) <- species_list

# stack each sp raster and sum for the final mapping
stack_bird <- rast(species_rasters)
bird_cooccurrence <- app(stack_bird, 
                            fun = sum, 
                            na.rm = TRUE)

# adding log transformation for mapping
bird_cooccurrence_log <- app(bird_cooccurrence, 
                                fun = function(x) { log1p(x) })

writeRaster(stack_bird,
            "LKDEBirdRasterStack.tif",
            overwrite = TRUE)

mapview(bird_cooccurrence_log,
        col.regions = viridis::inferno(256),
        na.color = "transparent",
        layer.name = "Core migratory areas")

#==================================================================
# Multispecies Line Kernel Density Estimation Overlay Nocturnal Migrant Birds
#===================================================================
library(terra)
library(sf)
library(dplyr)
library(purrr)
library(mapview)
library(stringr)
library(spatstat.geom)
library(spatstat.explore)

#creating flight lines density
flight_noc_bird.pj <- st_transform(flight_lines, crs = 3978) %>% 
  filter (
    MigrateTime == "nocturnal",
    Animal == "Bird"
  )


# create a blank template of 2.5 km for the region
regionTemplate <- rast(ext(flight_noc_bird.pj), 
                       resolution = 2500, 
                       crs = st_crs(flight_noc_bird.pj)$wkt)

st_write(flight_noc_bird.pj, 
         "flight_noc_bird.pj.gpkg", 
         delete_layer = TRUE)

# get the list of species to loop through
species_list <- unique(flight_noc_bird.pj$species)

# loop through the species to create the rasters
species_rasters <- lapply(species_list, 
                          function(sp_name) {
                            
                            df_sp <- flight_noc_bird.pj %>% 
                              filter(species == sp_name)
                            
                            flightsRast <- rasterize(vect(df_sp), 
                                                     regionTemplate, 
                                                     field = 1, 
                                                     fun = "sum", 
                                                     background = 0)
                            
                            #this is for the smoothing window
                            weightMatrix <- focalMat(flightsRast, 
                                                     d = 5000, 
                                                     type = "Gauss") 
                            kdeSurface <- focal(flightsRast, 
                                                w = weightMatrix, 
                                                fun = sum, 
                                                na.rm = TRUE)
                            
                            return(kdeSurface)
                          })

names(species_rasters) <- species_list

# stack each sp raster and sum for the final mapping
stack_noc_bird <- rast(species_rasters)
bird_noc_cooccurrence <- app(stack_noc_bird, 
                         fun = sum, 
                         na.rm = TRUE)

# adding log transformation for mapping
bird_noc_cooccurrence_log <- app(bird_noc_cooccurrence, 
                             fun = function(x) { log1p(x) })

writeRaster(stack_noc_bird,
            "LKDEBirdNocRasterStack.tif",
            overwrite = TRUE)

mapview(bird_noc_cooccurrence_log,
        col.regions = viridis::inferno(256),
        na.color = "transparent",
        layer.name = "Core migratory areas")

#==================================================================
# Multispecies Line Kernel Density Estimation Overlay Nocturnal Takeoff
#===================================================================
library(terra)
library(sf)
library(dplyr)
library(purrr)
library(mapview)
library(stringr)
library(spatstat.geom)
library(spatstat.explore)

#creating flight lines density
flight_noc_takeoff.pj <- st_transform(noc.takeoff, 
                                      crs = 3978) 


# create a blank template of 2.5 km for the region
regionTemplate <- rast(ext(flight_noc_takeoff.pj), 
                       resolution = 2500, 
                       crs = st_crs(flight_noc_takeoff.pj)$wkt)

st_write(flight_noc_takeoff.pj, 
         "flight_noc_takeoff.pj.gpkg", 
         delete_layer = TRUE)

# get the list of species to loop through
species_list <- unique(flight_noc_takeoff.pj$species)

# loop through the species to create the rasters
species_rasters <- lapply(species_list, 
                          function(sp_name) {
                            
                            df_sp <- flight_noc_takeoff.pj %>% 
                              filter(species == sp_name)
                            
                            flightsRast <- rasterize(vect(df_sp), 
                                                     regionTemplate, 
                                                     field = 1, 
                                                     fun = "sum", 
                                                     background = 0)
                            
                            #this is for the smoothing window
                            weightMatrix <- focalMat(flightsRast, 
                                                     d = 5000, 
                                                     type = "Gauss") 
                            kdeSurface <- focal(flightsRast, 
                                                w = weightMatrix, 
                                                fun = sum, 
                                                na.rm = TRUE)
                            
                            return(kdeSurface)
                          })

names(species_rasters) <- species_list

# stack each sp raster and sum for the final mapping
stack_noc_bird_takeoff <- rast(species_rasters)
bird_noc_takeoff_cooccurrence <- app(stack_noc_bird_takeoff, 
                             fun = sum, 
                             na.rm = TRUE)

# adding log transformation for mapping
bird_noc_takeoff_cooccurrence_log <- app(bird_noc_takeoff_cooccurrence, 
                                 fun = function(x) { log1p(x) })

writeRaster(stack_noc_bird_takeoff,
            "LKDEBirdNocTakeoffRasterStack.tif",
            overwrite = TRUE)

mapview(bird_noc_takeoff_cooccurrence_log,
        col.regions = viridis::inferno(256),
        na.color = "transparent",
        layer.name = "Core migratory areas")



#==================================================================
# Multispecies Line Kernel Density Estimation Overlay Nocturnal Landing
#===================================================================
library(terra)
library(sf)
library(dplyr)
library(purrr)
library(mapview)
library(stringr)
library(spatstat.geom)
library(spatstat.explore)

#creating flight lines density
flight_noc_land.pj <- st_transform(noc.landing, 
                                      crs = 3978) 


# create a blank template of 2.5 km for the region
regionTemplate <- rast(ext(flight_noc_land.pj), 
                       resolution = 2500, 
                       crs = st_crs(flight_noc_land.pj)$wkt)

st_write(flight_noc_land.pj, 
         "flight_noc_land.pj.gpkg", 
         delete_layer = TRUE)

# get the list of species to loop through
species_list <- unique(flight_noc_land.pj$species)

# loop through the species to create the rasters
species_rasters <- lapply(species_list, 
                          function(sp_name) {
                            
                            df_sp <- flight_noc_land.pj %>% 
                              filter(species == sp_name)
                            
                            flightsRast <- rasterize(vect(df_sp), 
                                                     regionTemplate, 
                                                     field = 1, 
                                                     fun = "sum", 
                                                     background = 0)
                            
                            #this is for the smoothing window
                            weightMatrix <- focalMat(flightsRast, 
                                                     d = 5000, 
                                                     type = "Gauss") 
                            kdeSurface <- focal(flightsRast, 
                                                w = weightMatrix, 
                                                fun = sum, 
                                                na.rm = TRUE)
                            
                            return(kdeSurface)
                          })

names(species_rasters) <- species_list

# stack each sp raster and sum for the final mapping
stack_noc_bird_land <- rast(species_rasters)
bird_noc_land_cooccurrence <- app(stack_noc_bird_land, 
                                     fun = sum, 
                                     na.rm = TRUE)

# adding log transformation for mapping
bird_noc_land_cooccurrence_log <- app(bird_noc_land_cooccurrence, 
                                         fun = function(x) { log1p(x) })

writeRaster(stack_noc_bird_land,
            "LKDEBirdNocLandRasterStack.tif",
            overwrite = TRUE)

mapview(bird_noc_land_cooccurrence_log,
        col.regions = viridis::inferno(256),
        na.color = "transparent",
        layer.name = "Core migratory areas")

#==================================================================
# Multispecies Line Kernel Density Estimation Overlay diurnal birds
#===================================================================
library(terra)
library(sf)
library(dplyr)
library(purrr)
library(mapview)
library(stringr)
library(spatstat.geom)
library(spatstat.explore)

day_birds.pj <- flight_lines.pj %>%
  filter (
    Animal == "Bird",
    MigrateTime == "diurnal"
    
  )

# create a blank template of 2.5 km for the region
regionTemplate <- rast(ext(day_birds.pj), 
                       resolution = 2500, 
                       crs = st_crs(day_birds.pj)$wkt)

st_write(day_birds.pj, 
         "day_birds.pj.gpkg", 
         delete_layer = TRUE)

# get the list of species to loop through
species_list <- unique(day_birds.pj$species)

# loop through the species to create the rasters
species_rasters <- lapply(species_list, 
                          function(sp_name) {
                            
                            df_sp <- day_birds.pj %>% 
                              filter(species == sp_name)
                            
                            flightsRast <- rasterize(vect(df_sp), 
                                                     regionTemplate, 
                                                     field = 1, 
                                                     fun = "sum", 
                                                     background = 0)
                            
                            #this is for the smoothing window
                            weightMatrix <- focalMat(flightsRast, 
                                                     d = 5000, 
                                                     type = "Gauss") 
                            kdeSurface <- focal(flightsRast, 
                                                w = weightMatrix, 
                                                fun = sum, 
                                                na.rm = TRUE)
                            
                            return(kdeSurface)
                          })

names(species_rasters) <- species_list

# stack each sp raster and sum for the final mapping
stack_day_birds <- rast(species_rasters)
bird_day_cooccurrence <- app(stack_day_birds, 
                                  fun = sum, 
                                  na.rm = TRUE)

# adding log transformation for mapping
bird_day_cooccurrence_log <- app(bird_cooccurrence, 
                                      fun = function(x) { log1p(x) })

writeRaster(stack_day_birds,
            "LKDEBirdDayRasterStack.tif",
            overwrite = TRUE)

mapview(bird_day_cooccurrence_log,
        col.regions = viridis::inferno(256),
        na.color = "transparent",
        layer.name = "Core migratory areas")

#======================================================================
# Multispecies Line Kernel Density Estimation Overlay bats
#======================================================================
library(terra)
library(sf)
library(dplyr)
library(purrr)
library(mapview)
library(stringr)
library(spatstat.geom)
library(spatstat.explore)

all_bats.pj <- flight_lines.pj %>%
  filter (
    Animal == "Bat"
  )

# create a blank template of 2.5 km for the region
regionTemplate <- rast(ext(all_bats.pj), 
                       resolution = 2500, 
                       crs = st_crs(all_bats.pj)$wkt)

st_write(all_bats.pj, 
         "all_bats.pj.gpkg", 
         delete_layer = TRUE)

# get the list of species to loop through
species_list <- unique(all_bats.pj$species)

# loop through the species to create the rasters
species_rasters <- lapply(species_list, 
                          function(sp_name) {
                            
                            df_sp <- all_bats.pj %>% 
                              filter(species == sp_name)
                            
                            flightsRast <- rasterize(vect(df_sp), 
                                                     regionTemplate, 
                                                     field = 1, 
                                                     fun = "sum", 
                                                     background = 0)
                            
                            #this is for the smoothing window
                            weightMatrix <- focalMat(flightsRast, 
                                                     d = 5000, 
                                                     type = "Gauss") 
                            kdeSurface <- focal(flightsRast, 
                                                w = weightMatrix, 
                                                fun = sum, 
                                                na.rm = TRUE)
                            
                            return(kdeSurface)
                          })

names(species_rasters) <- species_list

# stack each sp raster and sum for the final mapping
stack_bats <- rast(species_rasters)
bats_cooccurrence <- app(stack_bats, 
                             fun = sum, 
                             na.rm = TRUE)

# adding log transformation for mapping
bats_cooccurrence_log <- app(bats_cooccurrence, 
                                 fun = function(x) { log1p(x) })

writeRaster(stack_bats,
            "LKDEBatsRasterStack.tif",
            overwrite = TRUE)

mapview(bats_cooccurrence_log,
        col.regions = viridis::inferno(256),
        na.color = "transparent",
        layer.name = "Core migratory areas")


