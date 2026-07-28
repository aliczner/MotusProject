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

#======================================================
# line kernel density all flights
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

##sigma 10 km is the best output, use that going forward

# Save it but what I made was for plotting
writeRaster(rast_10km, "Plotting_flights_all_10km.tif")

#need to save the raw data 

density_data <- density.psp(
  flights_psp, 
  weights = marks(flights_psp), #this keeps the weights from above
  sigma = 10000, 
  eps = 1000
)

density_data.rast <- rast(density_data)
crs(density_data.rast) <- "EPSG:3978"
terra::writeRaster(density_data.rast, "Raw_flights_all.tif", 
                   overwrite=T)

#=============================================
# line kernel density - tagging outside GLWS
#============================================
flight_steps.sf <- flight_steps %>%
  filter(!is.na(lon_tagSite) & !is.na(lat_tagSite)) %>%
  st_as_sf(
    coords = c("lon_tagSite", "lat_tagSite"),
    crs = 4326,
    remove = FALSE
  )

GLWatershed_matched <- st_transform(GLWatershed, 
                                    crs = 4326)

tagged_outside <- st_filter(
  flight_steps.sf, 
  GLWatershed_matched, 
  .predicate = st_disjoint)


geoms <- lapply(seq_len(nrow(tagged_outside)), function(i) {
  st_linestring(rbind(
    c(tagged_outside$lon_previous[i], tagged_outside$lat_previous[i]),
    c(tagged_outside$lon[i],          tagged_outside$lat[i])
  ))
})


tagout_flight_lines <- st_sf(
  st_drop_geometry(tagged_outside),
  geometry = st_sfc(geoms, crs = 4326))

#creating flight lines density
tagout_flight.pj <- st_transform(tagout_flight_lines, crs = 3978)

tagout_flight.pj <- tagout_flight.pj %>%
  group_by(previousStationName, stationName) %>%
  mutate(
    path_freq = n(),
    # Square Root Distance Weighting
    weight_sqrt_dist = path_freq / (distance_km^0.5)
  ) %>% 
  ungroup()

tagout_psp <- as.psp(st_geometry(tagout_flight.pj))
marks(tagout_psp) <- tagout_flight.pj$weight_sqrt_dist

#line kernel density 

tagout.kde <- density.psp(
  tagout_psp, 
  weights = marks(tagout_psp), #this keeps the weights from above
  sigma = 10000, 
  eps = 1000
)
tagout.rast<- terra::rast(tagout.kde)
crs(tagout.rast) <- "EPSG:3978"

writeRaster(tagout.rast, "LineKDE_TaggedOut.tif")

#longpoint is still too hot to see anything else

# cutoff to avoid log(0) -> -Inf
tagout.rast[tagout.rast < 1e-10] <- NA

# Log transform the density values
tagout.rast_log <- log(tagout.rast)

# Set the lower bound anchor (log of 1e-7)
fixed_min <- log(1e-7)

# Get the maximum value directly from the log-transformed raster
l_max <- max(terra::values(tagout.rast_log), na.rm = TRUE)

# Rescale values onto the 0-100 index
r_index <- ((tagout.rast_log - fixed_min) / (l_max - fixed_min)) * 100

# Clamp baseline noise to 0
r_index[r_index <= 0] <- NA


# Render map
mapview(
  r_index,
  col.regions = viridis::inferno(100)
)
