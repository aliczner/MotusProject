# code for investigating spring migration phenology

library(dplyr)
library(lubridate)
library(ggplot2)
library(ggridges)
library(sf)

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

#======================================================
# line kernel density
# ====================================================

library(spatstat.geom)
library(spatstat.explore)
library(terra)

#creating flight lines density
flight_lines.pj <- st_transform(flight_lines, crs = 3978)

flight_lines.pj <- flight_lines.pj %>%
  group_by(previousStationName, stationName) %>%
  mutate(
    path_freq = n(),
    # Square Root Distance Weighting (Power = 0.5)
    weight_sqrt_dist = path_freq / (distance_km^0.5)
  ) %>% 
  ungroup()

flights_psp <- as.psp(st_geometry(flight_lines.pj))
marks(flights_psp) <- flight_lines.pj$weight_sqrt_dist

### testing different sigma values
process_flight_density <- function(psp_obj, 
                                   sigma_meters, 
                                   envelope_shape, 
                                   eps_meters = 1000, 
                                   fixed_min_val = 1e-7) {
  
  # Estimate kernel density weighted by step distance (spatstat)
  dens_im <- density.psp(
    psp_obj, 
    weights = marks(psp_obj),
    sigma = sigma_meters, 
    eps = eps_meters
  )
  
  # Convert spatstat image to SpatRaster
  r <- rast(dens_im)
  crs(r) <- "EPSG:3978"
  
  # cutoff to avoid log(0) -> -Inf
  r[r < 1e-10] <- NA
  
  # log transform
  r_log <- log(r)
  
  # Set visual baseline anchor & calculate max
  fixed_min <- log(fixed_min_val)
  l_max     <- max(values(r_log), na.rm = TRUE)
  
  # Rescale from 0 to 100
  r_index <- ((r_log - fixed_min) / (l_max - fixed_min)) * 100
  
  # make aint background values (< 1e-7) into 0
  r_index[r_index < 0] <- 0
  
  # for the envlope make internal values NA
  env_vect <- vect(envelope_shape)
  env_rast <- rasterize(env_vect, r_index, field = 0)
  
  r_index <- cover(r_index, env_rast)
  r_index <- mask(r_index, env_vect)
  
  return(r_index)
}

# Generate layers for 10 km, 20 km, and 30 km bandwidths
rast_10km <- process_flight_density(
  psp_obj = flights_psp, 
  weights_vec = flight_lines.pj$route_weight, 
  sigma_meters = 10000, 
  envelope_shape = flight_envelope
)

rast_20km <- process_flight_density(
  psp_obj = flights_psp, 
  weights_vec = flight_lines.pj$route_weight, 
  sigma_meters = 20000, 
  envelope_shape = flight_envelope
)

rast_30km <- process_flight_density(
  psp_obj = flights_psp, 
  weights_vec = flight_lines.pj$route_weight, 
  sigma_meters = 30000, 
  envelope_shape = flight_envelope
)

mapview(
  rast_10km, 
  col.regions = viridis::inferno(100), 
  layer.name = "Relative Movement Density (Sigma 10 km)",
  na.color = "transparent"
) +
  mapview(
    rast_20km, 
    col.regions = viridis::inferno(100), 
    layer.name = "Relative Movement Density (Sigma 20 km)",
    na.color = "transparent"
  ) +
  mapview(
    rast_30km, 
    col.regions = viridis::inferno(100), 
    layer.name = "Relative Movement Density (Sigma 30 km)",
    na.color = "transparent"
  )

##sigma 10 km is the best output

# Save it
terra::writeRaster(rast_10km, "flights_all_10km.tif")
