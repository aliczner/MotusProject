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
library(patchwork)
library(sf)
library(terra)
library(tidyterra)

#defining the study are extent
bbox_wgs84 <- sf::st_bbox(c(xmin = -92.7, 
                            xmax = -72.9, 
                            ymin = 37, 
                            ymax = 49), 
                          crs = 4326)

# getting the ESRI basemap for plotting
tiles <- maptiles::get_tiles(x = bbox_wgs84, 
                             provider = "Esri.WorldGrayCanvas", 
                             zoom = 8, 
                             crop = TRUE)
b <- sf::st_bbox(tiles)

### making a bubble plot for detections

spring_station_coords <- bind_rows(
  flight_lines %>% 
    st_drop_geometry() %>% 
    select(lon, lat),
  
  flight_lines %>% 
    st_drop_geometry() %>% 
    select(lon = lon_previous, lat = lat_previous)
) %>% 
  filter(!is.na(lon) & !is.na(lat)) %>%
  group_by(lon, lat) %>%
  summarise(total_detections = n(), .groups = "drop") %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) 

# make the plot
spring_bubble_map <- ggplot() +
  tidyterra::geom_spatraster_rgb(data = tiles) +
  geom_sf(data = spring_station_coords, 
          aes(size = total_detections), 
          color = "#0F0E0E", 
          fill = "#5B1865", 
          alpha = 0.6, 
          shape = 21) +
  scale_size_continuous(
    range = c(1.2, 6), 
    name = "Station detections"
  ) +
  coord_sf(
    crs = 4326, 
    xlim = c(b["xmin"], b["xmax"]), 
    ylim = c(b["ymin"], b["ymax"]), 
    expand = FALSE
  ) +
  labs(
    title = "A",
    x = "", 
    y = ""
  )+
  theme_minimal() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.99,0.99),
    legend.justification = c(1,1),
    legend.direction = "horizontal",
    legend.title.position = "top",
    legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  )

### tagging site effort plot

#getting the coordinates
spring_tagging_plot_df <- flight_lines %>%
  st_drop_geometry() %>%
  select(lon = lon_tagSite, lat = lat_tagSite) %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  count(lon, lat, name = "num_tagged") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  mutate(x = st_coordinates(.)[, 1]) %>%
  st_drop_geometry()

# making the effort plot

spring_tagging_lon <- ggplot(spring_tagging_plot_df, 
                           aes(x = x, 
                               y = num_tagged)) +
  geom_segment(aes(xend = x, y = 0, yend = num_tagged),
               color = "black", linewidth = 2) +
  theme_minimal() +
  coord_cartesian(xlim = unname(c(b["xmin"], 
                                  b["xmax"])), 
                  expand = FALSE) +
  labs(title = "B",
       x = "longitude", 
       y = "Individuals Tagged") +
  theme(panel.grid.minor = element_blank())

#putting both figures together
xlims = c(b["xmin"], b["xmax"]) 
ylims = c(b["ymin"], b["ymax"])
asp_map <- diff(ylims) / (diff(xlims) * cos(mean(ylims) * pi / 180))
asp_b   <- asp_map * 0.35   # plot B height relative to the map; change 0.35 to taste

spring_tagging_lon <- spring_tagging_lon +
  theme(aspect.ratio = asp_b)

spring_composite_figure <- spring_bubble_map / spring_tagging_lon +
  plot_layout(heights = c(asp_map, asp_b)) &
  theme(plot.margin = margin(3, 5, 3, 5))

# Save at a size that fits the map's shape
w <- 9
h <- (w - 1) * (asp_map + asp_b) + 1.2   # ~1" for the y-axis text, ~1.2" for titles and the x-axis
ggsave("spring_composite.png", spring_composite_figure, 
       width = w, 
       height = h, 
       dpi = 300)
ggsave("spring_composite.pdf",
       spring_composite_figure, 
       width = w,
       height = h, 
       dpi = 300)

# just search effort
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


