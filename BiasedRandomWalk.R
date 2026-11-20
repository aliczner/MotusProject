# script to create biased correlated random walks for sensitivity testing
# will create one for spring, and one for fall

#====================================================================
#  Calibrate directional bias (kappa) from empirical flight bearings
#====================================================================

library(dplyr)
library(sf)
library (lubridate)

#kappa = a measure of how closely the bird sticks to North/South
# k = 0, completely random direction
# k = 0.5 - 1, heads loosely north but pretty wide
# k = 1.5 - 3, mainly northward flight, some deviations
# k = >5 sticks to North

calculate_kappa <- function(flights_sf, season = "Spring") {
  
  # Set the main migration direction in radians:
  #    0 radians = Due North (Spring)
  #    pi radians = Due South (Fall)
  target_rad <- if (grepl("Spring", 
                          season, 
                          ignore.case = TRUE)) 0 else pi
  
  # Drop spatial geometry 
  df <- st_drop_geometry(flights_sf)
  
  # Convert degrees (coordinates) to radians 
  rad  <- pi / 180
  phi1 <- df$lat_previous * rad # previous tower latitude
  phi2 <- df$lat * rad # current tower latitude
  dlon <- (df$lon - df$lon_previous) * rad # Change in longitude to account
  #for longitude distortion at the poles
  
  # Calculate the compass bearing (in radians, -pi to +pi) for every row
  bearings <- atan2(
    sin(dlon) * cos(phi2), #spherical law of sines (E/W)
    cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(dlon) ## Spherical cosine identity (N/S)
  )
  
  # Measure how aligned the real flights are with due North (Spring) or due South (Fall).
  # R_proj ranges from 0 (no directional alignment) to 1 (straight due North/South).
  R_proj <- max(0.01, 
                min(0.95, 
                    mean(cos(bearings - target_rad), 
                         na.rm = TRUE)))
  
  # R_proj into the von Mises concentration parameter 'kappa' (Banerjee et al. 2005)
  kappa <- (R_proj * (2 - R_proj^2)) / (1 - R_proj^2)
  
  message(sprintf("[%s] Calculated kappa = %.2f (directional alignment = %.2f)",
                  season, kappa, R_proj))
  return(kappa)
}

##############################
# running the kappa calculation function

#spring nocturnal birds
spring_noc_kappa <- calculate_kappa(flights_sf = flight_noc_bird.pj, 
                              season = "Spring Migration")
#spring nocturnal takeoff
spring_takeoff_kappa <- calculate_kappa(flights_sf = flight_noc_takeoff.pj, 
                                        season = "Spring Migration")
#spring nocturnal landing
spring_land_kappa <- calculate_kappa(flights_sf = flight_noc_land.pj, 
                                        season = "Spring Migration")
#spring diurnal birds
spring_day_kappa <- calculate_kappa(flights_sf = day_birds.pj, 
                                     season = "Spring Migration")
#spring bats
spring_bats_kappa <- calculate_kappa(flights_sf = all_bats.pj, 
                                    season = "Spring Migration")
#fall nocturnal birds
fall_noc_kappa <- calculate_kappa(flights_sf = fall_noc_bird, 
                                    season = "Fall Migration")
#Fall nocturnal takeoff
fall_takeoff_kappa <- calculate_kappa(flights_sf = fall_takeoff_bird, 
                                        season = "Fall Migration")
#Fall nocturnal landing
fall_land_kappa <- calculate_kappa(flights_sf = fall_land_bird, 
                                     season = "Fall Migration")
#Fall diurnal birds
fall_day_kappa <- calculate_kappa(flights_sf = fall_day_bird, 
                                    season = "Fall Migration")
#Fall bats
fall_bats_kappa <- calculate_kappa(flights_sf = fall_bats, 
                                     season = "Fall Migration")

#====================================================================
# Simulate tower-to-tower biased random walk 
#====================================================================
#based on the first point in the GLWS not tagging location

simulate_tower_brw <- function(flights_sf,
                               effort_df,
                               glws_poly,
                               season = "Spring",
                               kappa = NULL,
                               max_step_km = NULL,
                               seed = 42) {
  set.seed(seed)
  target_bearing <- if (grepl("Spring", 
                              season, 
                              ignore.case = TRUE)) 0 else pi
  # kappa if missing
  if (is.null(kappa)) {
    kappa <- calculate_kappa(flights_sf, 
                             season = season)
  }
  
  # fitting step lengths of true flights to a gamma distribution
  # sort chronologically by individual, get their first Great Lakes entry point
  emp_df <- st_drop_geometry(flights_sf) %>%
    arrange(tagDeployID, as.POSIXct(tsEnd_dt))
  emp_km <- emp_df$distance_km
  #remove any NAs and any zeros
  emp_km <- emp_km[!is.na(emp_km) & emp_km > 1] 
  
  # Use 99th percentile of flights as the max for the distribution
  step_cap <- if (is.null(max_step_km)) quantile(emp_km, 
                                                 0.99, 
                                                 na.rm = TRUE) else max_step_km
  emp_km <- emp_km[emp_km <= step_cap]
  
  #this determines where the curve's peak is
  g_shape <- (mean(emp_km, na.rm = TRUE) / sd(emp_km, na.rm = TRUE))^2
  #this determines the scale of the axis 
  g_scale <- (sd(emp_km, na.rm = TRUE)^2) / mean(emp_km, na.rm = TRUE)
  
  # Summarize each animal's Great Lakes entry spot, species, year, and number of steps
  birds <- emp_df %>%
    filter(!is.na(lon_previous), 
           !is.na(lat_previous), !is.na(year_end)) %>%
    group_by(tagDeployID, 
             species, 
             year_end, 
             lon_tagSite, 
             lat_tagSite) %>%
    summarise(
      start_lon = first(lon_previous),
      start_lat = first(lat_previous),
      n_steps   = n(),
      first_ts  = min(as.POSIXct(tsEnd_dt), na.rm = TRUE),
      .groups   = "drop"
    )
  
  # Pre-filter active towers by year once so it doesnt re-filter inside the loop
  valid_towers <- effort_df %>%
    filter(active_days >= 1, !is.na(longitude), !is.na(latitude)) %>%
    mutate(log_effort = log(active_days + 1))
  towers_by_year <- split(valid_towers, valid_towers$year)
  
  # make emplty vectors to store steps fast
  max_segs <- sum(birds$n_steps)
  out_tagDeployID <- vector("character", max_segs)
  out_species <- vector("character", max_segs)
  out_stationID <- vector("integer", max_segs)
  out_lon_previous <- vector("numeric", max_segs)
  out_lat_previous <- vector("numeric", max_segs)
  out_lon <- vector("numeric", max_segs)
  out_lat <- vector("numeric", max_segs)
  out_lon_tagSite <- vector("numeric", max_segs)
  out_lat_tagSite <- vector("numeric", max_segs)
  out_year_end <- vector("integer", max_segs)
  out_distance_km <- vector("numeric", max_segs)
  out_tsEnd_dt <- rep(birds$first_ts[1], max_segs)
  
  rad <- pi / 180
  k <- 0L
  
  for (b in seq_len(nrow(birds))) {
    # Active towers in this individual's migration year
    #to make sure they didn't fly to towers that didn't exist yet
    towers_yr <- towers_by_year[[as.character(birds$year_end[b])]]
    if (is.null(towers_yr) || nrow(towers_yr) < 2) next
    
    t_id      <- towers_yr$stationID
    t_lon     <- towers_yr$longitude
    t_lat     <- towers_yr$latitude
    t_log_eff <- towers_yr$log_effort
    
    # Start the walk at the bird's first observed Great Lakes entry point
    curr_lon <- birds$start_lon[b]
    curr_lat <- birds$start_lat[b]
    curr_id  <- NA_integer_
    
    for (s in seq_len(birds$n_steps[b])) {
      # Exclude the station the animal is currently standing at so they
      #have to move
      not_curr <- if (is.na(curr_id)) rep(TRUE, length(t_id)) else (t_id != curr_id)
      
      # Great-circle Haversine distance (km) from current location to 
      #all candidate stations
      dlat <- (t_lat - curr_lat) * rad
      dlon <- (t_lon - curr_lon) * rad
      #accounts for the curvature of the earth
      a <- sin(dlat / 2)^2 + cos(curr_lat * rad) * cos(t_lat * rad) * sin(dlon / 2)^2
      #converts a back to rad and then km
      dist_km <- 2 * 6371 * asin(pmin(1, sqrt(a)))
      
      # Keep candidate towers within min/max
      idx <- which(not_curr & dist_km > 1 & dist_km <= step_cap)
      if (length(idx) == 0) break
      #keep only matching rows within the min/max steps
      d_sub <- dist_km[idx]
      c_lon <- t_lon[idx]
      c_lat <- t_lat[idx]
      
      # Compass bearing (radians) from current location to each candidate tower
      phi1 <- curr_lat * rad
      phi2 <- c_lat * rad
      dlon_sub <- (c_lon - curr_lon) * rad
      bearing <- atan2(
        sin(dlon_sub) * cos(phi2),
        cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(dlon_sub)
      )
      
      # Assign each station a probability score
      # Selection weight: directional bias 
      #and step-length likelihood 
      # also considering station effort
      w <- exp(kappa * cos(bearing - target_bearing)) *
        dgamma(d_sub, shape = g_shape, scale = g_scale) *
        t_log_eff[idx]
      
      #this looks at the weights and selects one station
      if (sum(w, na.rm = TRUE) <= 0) break
      chosen_pos <- sample.int(length(idx), 1, prob = w)
      chosen_idx <- idx[chosen_pos]
      
      # Store plain values for this chosen step (avoids slow sf creation inside loop)
      k <- k + 1L
      out_tagDeployID[k]  <- as.character(birds$tagDeployID[b])
      out_species[k]      <- as.character(birds$species[b])
      out_stationID[k]    <- t_id[chosen_idx]
      out_lon_previous[k] <- curr_lon
      out_lat_previous[k] <- curr_lat
      out_lon[k]          <- t_lon[chosen_idx]
      out_lat[k]          <- t_lat[chosen_idx]
      out_lon_tagSite[k]  <- birds$lon_tagSite[b]
      out_lat_tagSite[k]  <- birds$lat_tagSite[b]
      out_year_end[k]     <- birds$year_end[b]
      out_distance_km[k]  <- d_sub[chosen_pos]
      out_tsEnd_dt[k]     <- birds$first_ts[b]
      
      # Update current position and station ID for the next step
      curr_lon <- t_lon[chosen_idx]
      curr_lat <- t_lat[chosen_idx]
      curr_id  <- t_id[chosen_idx]
    }
  }
  
  # Stop if nothing simulated
  if (k == 0L) stop("No simulated tracks generated—check year matching and tower coordinates.")
  
  # Trim previous vectors to actual simulated segments
  seq_k <- seq_len(k)
  sim_df <- data.frame(
    tagDeployID = out_tagDeployID[seq_k],
    species  = out_species[seq_k],
    stationID = out_stationID[seq_k],
    lon_previous = out_lon_previous[seq_k],
    lat_previous = out_lat_previous[seq_k],
    lon = out_lon[seq_k],
    lat = out_lat[seq_k],
    lon_tagSite = out_lon_tagSite[seq_k],
    lat_tagSite = out_lat_tagSite[seq_k],
    year_end = out_year_end[seq_k],
    distance_km = out_distance_km[seq_k],
    tsEnd_dt = out_tsEnd_dt[seq_k],
    stringsAsFactors = FALSE
  )
  
  # Build all spatial linestrings at once instead of one segment at a time
  geoms <- lapply(seq_k, function(j) {
    st_linestring(matrix(c(sim_df$lon_previous[j], sim_df$lon[j],
                           sim_df$lat_previous[j], sim_df$lat[j]), ncol = 2))
  })
  sim_all <- st_sf(sim_df, geometry = st_sfc(geoms, crs = 4326))
  
  # Keep only simulated flights with at least one station inside the GLWS
  glws_proj <- if (st_crs(glws_poly) == st_crs(3348)) glws_poly else st_transform(st_union(glws_poly), crs = 3348)
  end_pts   <- st_as_sf(sim_df, 
                        coords = c("lon", "lat"), 
                        crs = 4326) %>% st_transform(3348)
  start_pts <- st_as_sf(sim_df, 
                        coords = c("lon_previous", 
                                   "lat_previous"), 
                        crs = 4326) %>% st_transform(3348)
  
  keep <- (lengths(st_intersects(end_pts, glws_proj)) > 0) |
    (lengths(st_intersects(start_pts, glws_proj)) > 0)
  
  # Reproject filtered tracks to match flights_sf CRS
  sim_all[keep, ] %>% st_transform(crs = st_crs(flights_sf))
}

####################################
#running the simulate tower function

# Load and prepare the Great Lakes watershed polygon
GLWatershed <- st_read("./greatlakes_subbasins/greatlakes_subbasins.shp")
GLWS_union  <- st_transform(st_union(GLWatershed), crs = 3348)

# Load effort files
spring_effort_by_year <- readRDS("SpringDepsDuration.rds")
fall_effort_by_year <- readRDS("FallDepsDuration.rds")

# Single test run for Spring nocturnal birds
sim_spring_noc <- simulate_tower_brw(
  flights_sf = flight_noc_bird.pj,
  effort_df  = spring_effort_by_year,
  glws_poly = GLWS_union,
  season  = "Spring Migration",
  kappa = spring_noc_kappa,
  seed = 42
) # yay it worked! so run 200 reps

# Run 200 null replicates for Spring nocturnal birds
n_reps <- 200

#adding a progress bar because it takes awhile
sim_reps_spring_noc <- lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0 || i == 1) cat(sprintf("Running replicate %d of %d...\n", i, n_reps))
  
  simulate_tower_brw(
    flights_sf = flight_noc_bird.pj,
    effort_df = spring_effort_by_year,
    glws_poly = GLWS_union,
    season = "Spring Migration",
    kappa = spring_noc_kappa,
    seed = 1000 + i
  ) %>%
    mutate(sim_rep = i)
})

sim_all_spring_noc <- bind_rows(sim_reps_spring_noc)
saveRDS(sim_all_spring_noc, "sim_all_spring_noc.rds")

### spring nocturnal takeoff

# Run 200 null replicates
n_reps <- 200

#adding a progress bar because it takes awhile
sim_reps_spring_takeoff <- lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0 || i == 1) cat(sprintf("Running replicate %d of %d...\n", 
                                          i, 
                                          n_reps))
  simulate_tower_brw(
    flights_sf = flight_noc_takeoff.pj,
    effort_df = spring_effort_by_year,
    glws_poly = GLWS_union,
    season = "Spring Migration",
    kappa = spring_takeoff_kappa,
    seed = 1000 + i
  ) %>%
    mutate(sim_rep = i)
})

sim_spring_takeoff <- bind_rows(sim_reps_spring_takeoff)
saveRDS(sim_spring_takeoff, "sim_spring_takeoff.rds")

### spring nocturnal landing

# Run 200 null replicates
n_reps <- 200

sim_reps_spring_land <- lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0 || i == 1) cat(sprintf("Running replicate %d of %d...\n", 
                                          i, 
                                          n_reps))
  simulate_tower_brw(
    flights_sf = flight_noc_land.pj,
    effort_df = spring_effort_by_year,
    glws_poly = GLWS_union,
    season = "Spring Migration",
    kappa = spring_land_kappa,
    seed = 1000 + i
  ) %>%
    mutate(sim_rep = i)
})

sim_spring_land <- bind_rows(sim_reps_spring_land)
saveRDS(sim_spring_land, "sim_spring_land.rds")


### spring diurnal birds

# Run 200 null replicates
n_reps <- 200

sim_reps_spring_day <- lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0 || i == 1) cat(sprintf("Running replicate %d of %d...\n", 
                                          i, 
                                          n_reps))
  simulate_tower_brw(
    flights_sf = day_birds.pj,
    effort_df = spring_effort_by_year,
    glws_poly = GLWS_union,
    season = "Spring Migration",
    kappa = spring_day_kappa,
    seed = 1000 + i
  ) %>%
    mutate(sim_rep = i)
})

sim_spring_day <- bind_rows(sim_reps_spring_day)
saveRDS(sim_spring_day, "sim_spring_day.rds")


###  spring bats

# Run 200 null replicates
n_reps <- 200

sim_reps_spring_bats <- lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0 || i == 1) cat(sprintf("Running replicate %d of %d...\n", 
                                          i, 
                                          n_reps))
  simulate_tower_brw(
    flights_sf = all_bats.pj,
    effort_df = spring_effort_by_year,
    glws_poly = GLWS_union,
    season = "Spring Migration",
    kappa = spring_bats_kappa,
    seed = 1000 + i
  ) %>%
    mutate(sim_rep = i)
})

sim_spring_bats <- bind_rows(sim_reps_spring_bats)
saveRDS(sim_spring_bats, "sim_spring_bats.rds")

####  Fall season

# Run 200 null replicates for Fall nocturnal birds
n_reps <- 200

#adding a progress bar because it takes awhile
sim_reps_fall_noc <- lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0 || i == 1) cat(sprintf("Running replicate %d of %d...\n", i, n_reps))
  
  simulate_tower_brw(
    flights_sf = fall_noc_bird,
    effort_df = fall_effort_by_year,
    glws_poly = GLWS_union,
    season = "Fall Migration",
    kappa = fall_noc_kappa,
    seed = 1000 + i
  ) %>%
    mutate(sim_rep = i)
})

sim_all_fall_noc <- bind_rows(sim_reps_fall_noc)
saveRDS(sim_all_fall_noc, "sim_all_fall_noc.rds")

### Fall nocturnal takeoff

# Run 200 null replicates
n_reps <- 200

#adding a progress bar because it takes awhile
sim_reps_fall_takeoff <- lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0 || i == 1) cat(sprintf("Running replicate %d of %d...\n", 
                                          i, 
                                          n_reps))
  simulate_tower_brw(
    flights_sf = fall_takeoff_bird,
    effort_df = fall_effort_by_year,
    glws_poly = GLWS_union,
    season = "Fall Migration",
    kappa = fall_takeoff_kappa,
    seed = 1000 + i
  ) %>%
    mutate(sim_rep = i)
})

sim_fall_takeoff <- bind_rows(sim_reps_fall_takeoff)
saveRDS(sim_fall_takeoff, "sim_fall_takeoff.rds")

### Fall nocturnal landing

# Run 200 null replicates
n_reps <- 200

sim_reps_fall_land <- lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0 || i == 1) cat(sprintf("Running replicate %d of %d...\n", 
                                          i, 
                                          n_reps))
  simulate_tower_brw(
    flights_sf = fall_land_bird,
    effort_df = fall_effort_by_year,
    glws_poly = GLWS_union,
    season = "Fall Migration",
    kappa = fall_land_kappa,
    seed = 1000 + i
  ) %>%
    mutate(sim_rep = i)
})

sim_fall_land <- bind_rows(sim_reps_fall_land)
saveRDS(sim_fall_land, "sim_fall_land.rds")


### Fall diurnal birds

# Run 200 null replicates
n_reps <- 200

sim_reps_fall_day <- lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0 || i == 1) cat(sprintf("Running replicate %d of %d...\n", 
                                          i, 
                                          n_reps))
  simulate_tower_brw(
    flights_sf = fall_day_bird,
    effort_df = fall_effort_by_year,
    glws_poly = GLWS_union,
    season = "Fall Migration",
    kappa = fall_day_kappa,
    seed = 1000 + i
  ) %>%
    mutate(sim_rep = i)
})

sim_fall_day <- bind_rows(sim_reps_fall_day)
saveRDS(sim_fall_day, "sim_fall_day.rds")


###  Fall bats

# Run 200 null replicates
n_reps <- 200

sim_reps_fall_bats <- lapply(seq_len(n_reps), function(i) {
  if (i %% 10 == 0 || i == 1) cat(sprintf("Running replicate %d of %d...\n", 
                                          i, 
                                          n_reps))
  simulate_tower_brw(
    flights_sf = fall_bats,
    effort_df = fall_effort_by_year,
    glws_poly = GLWS_union,
    season = "Fall Migration",
    kappa = fall_bats_kappa,
    seed = 1000 + i
  ) %>%
    mutate(sim_rep = i)
})

sim_fall_bats <- bind_rows(sim_reps_fall_bats)
saveRDS(sim_fall_bats, "sim_fall_bats.rds")
