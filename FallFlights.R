# code for investigating fall migration phenology

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

fallData <- fullData %>% 
  filter(
    season == "Fall Migration",
  )
nrow(fallData) #10399 obs

fallTable <- fallData %>%  
  mutate(
    tsEnd_dt = as_datetime(tsEnd_dt, tz = "GMT"),
    # merging georgian bay with lake huron so they are not separate. 
    subbasin = if_else(subbasin == "geo_bay", "lk_huron", subbasin)
  )

#=========================================================
# creating flight paths
#=========================================================

library(sfheaders)

fall_flight_steps <- fallTable %>%
  filter(flight_type != "incidence") %>%
  group_by(tagDeployID, flight_ID) %>%
  arrange(tsEnd_dt, .by_group = TRUE) %>%
  ungroup()

geoms <- lapply(seq_len(nrow(fall_flight_steps)), function(i) {
  st_linestring(matrix(
    c(fall_flight_steps$lon_previous[i], fall_flight_steps$lon[i],
      fall_flight_steps$lat_previous[i], fall_flight_steps$lat[i]),
    ncol = 2
  ))
})
fall_flight_lines <- st_sf(fall_flight_steps, 
                      geometry = st_sfc(geoms, 
                                        crs = 4326))

#===========================================================
# creating a figure of search effort
# ======================================================

#plots of just effort

fall_flight_lines %>%
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
    title = "Longitude Distribution - Fall",
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

fall_flight_time_summary <- fall_flight_lines %>%
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
fall_plot_diurnal <- fall_flight_time_summary %>%
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
    title = "Fall Diurnal",
    x = "",
    y = "Percentage of Flights (%)",
    fill = "Diel Period"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# Mixed  Plot
fall_plot_mixed <- fall_flight_time_summary %>%
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
    title = "Fall Mixed",
    x = "",
    y = "Percentage of Flights (%)",
    fill = "Diel Period"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# nocturnal
fall_plot_nocturnal <- fall_flight_time_summary %>%
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
    title = "Fall Nocturnal",
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
fall_noc.takeoff <- fall_flight_lines %>%  
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

fall_noc.landing <- fall_flight_lines %>%  
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
# Multispecies Line Kernel Density Estimation Function
#===================================================================
library(terra)
library(sf)
library(dplyr)
library(purrr)
library(mapview)
library(stringr)
library(spatstat.geom)
library(spatstat.explore)

generate_lkde_surface <- function(data, 
                                  animal_group, 
                                  migrate_time, 
                                  output_prefix = "migration",
                                  target_crs = 3978) {
  
  # Filter data based on function arguments
  filtered_data <- data %>%
    filter(
      Animal == animal_group,
      MigrateTime == migrate_time
    )
  
  if (nrow(filtered_data) == 0) {
    stop("No rows match the specified Animal and MigrateTime combination.")
  }
  
  # Automatically reproject to a projected coordinate system in meters (EPSG 3978) 
  # if it's currently in geographic coordinates (lat/long)
  if (st_is_longlat(filtered_data) || is.na(st_crs(filtered_data)$epsg) || st_crs(filtered_data)$epsg != target_crs) {
    message(sprintf("Transforming input data to CRS EPSG:%d (meters)...", target_crs))
    filtered_data <- st_transform(filtered_data, crs = target_crs)
  }
  
  # Create a blank template of 2.5 km for the region
  regionTemplate <- rast(ext(filtered_data), 
                         resolution = 2500, 
                         crs = st_crs(filtered_data)$wkt)
  
  # write out filtered layer if needed
  gpkg_name <- paste0(output_prefix, 
                      "_", 
                      animal_group, 
                      "_", 
                      migrate_time,
                      ".gpkg")
  st_write(filtered_data, 
           gpkg_name, 
           delete_layer = TRUE,
           quiet = TRUE)
  
  # list of species to loop through
  species_list <- unique(filtered_data$species)
  
  # Loop through the species to create individual KDE rasters
  species_rasters <- lapply(species_list, function(sp_name) {
    
    df_sp <- filtered_data %>% 
      filter(species == sp_name)
    
    flightsRast <- rasterize(vect(df_sp), 
                             regionTemplate, 
                             field = 1, 
                             fun = "sum", 
                             background = 0)
    
    # SAFETY CHECK: Ensure the raster is large enough for the 5km focal window
    if (nrow(flightsRast) > 2 && ncol(flightsRast) > 2) {
      # Smoothing window (Gaussian Kernel)
      weightMatrix <- focalMat(flightsRast, 
                               d = 5000, 
                               type = "Gauss") 
      
      kdeSurface <- focal(flightsRast, 
                          w = weightMatrix, 
                          fun = sum, 
                          na.rm = TRUE)
    } else {
      # If the extent is too tiny, return the unsmoothed raster with a message
      message(sprintf("Note: Species '%s' has a spatial extent too small for 5km smoothing. Skipping focal smoothing.", sp_name))
      kdeSurface <- flightsRast
    }
    
    return(kdeSurface)
  })
  
  names(species_rasters) <- species_list
  
  # Stack each species raster and sum for co-occurrence
  stack_rasters <- rast(species_rasters)
  cooccurrence_surface <- app(stack_rasters, 
                              fun = sum, 
                              na.rm = TRUE)
  
  # Log transformation for mapping
  cooccurrence_log <- app(cooccurrence_surface, 
                          fun = function(x) { log1p(x) })
  
  # Save raster stack output
  tiff_name <- paste0("LKDE_RasterStack_Fall_", 
                      animal_group, 
                      "_", 
                      migrate_time, 
                      ".tif")
  writeRaster(stack_rasters, 
              tiff_name,
              overwrite = TRUE)
  
  # Return a list
  return(list(
    filtered_data = filtered_data,
    species_stack = stack_rasters,
    cooccurrence_surface = cooccurrence_surface,
    cooccurrence_log = cooccurrence_log
  ))
}

## running the function

# fall all nocturnal birds

fall_night_bird_results <- generate_lkde_surface(
  data = fall_flight_lines, 
  animal_group = "Bird", 
  migrate_time = "nocturnal",
  output_prefix= "fall_all"
)

mapview(fall_night_bird_results$cooccurrence_log,
        col.regions = viridis::inferno(256),
        na.color = "transparent",
        layer.name = "Core migratory areas")

# fall nocturnal takeoff

fall_takeoff_bird_results <- generate_lkde_surface(
  data = fall_noc.takeoff, 
  animal_group = "Bird", 
  migrate_time = "nocturnal",
  output_prefix= "fall_takeoff"
)

mapview(fall_takeoff_bird_results$cooccurrence_log,
        col.regions = viridis::inferno(256),
        na.color = "transparent",
        layer.name = "Core migratory areas")

#fall nocturnal landing

fall_land_bird_results <- generate_lkde_surface(
  data = fall_noc.landing, 
  animal_group = "Bird", 
  migrate_time = "nocturnal",
  output_prefix= "fall_land"
)

mapview(fall_land_bird_results$cooccurrence_log,
        col.regions = viridis::inferno(256),
        na.color = "transparent",
        layer.name = "Core migratory areas")

# fall diurnal migrants

fall_day_bird_results <- generate_lkde_surface(
  data = fall_flight_lines, 
  animal_group = "Bird", 
  migrate_time = "diurnal",
  output_prefix= "fall_all"
)

mapview(fall_day_bird_results$cooccurrence_log,
        col.regions = viridis::inferno(256),
        na.color = "transparent",
        layer.name = "Core migratory areas")

#bats

fall_bat_results <- generate_lkde_surface(
  data = fall_flight_lines, 
  animal_group = "Bat", 
  migrate_time = "unclassified",
  output_prefix= "fall_bat"
)

mapview(fall_bat_results$cooccurrence_log,
        col.regions = viridis::inferno(256),
        na.color = "transparent",
        layer.name = "Core migratory areas")
