## Rerunning the prioritizr analysis with sensitivity testing
## by removing any individuals tagged at longpoint

#=================================================
# identifying the longpoint birds
#=================================================
library(dplyr)
library (sf)
library(lubridate)

#creating a dataframe containing the lat/lon for the towers at longpoint
lp_towers <- data.frame(
    lon = c(-80.3966, -80.3966, -80.3981, -80.2844, -80.3987, -80.3983,
            -80.3852, -80.0559),
    lat = c(42.5832, 42.5833, 42.5828, 42.5611, 42.5828, 42.5829, 42.5798, 
            42.5467)
  )


# Convert towers to an sf object (EPSG:3978)
lp_towers_sf <- lp_towers %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = st_crs(flight_lines.pj))

# Create a 2 km (2000 m) buffer around the 8 towers
lp_buffer <- st_buffer(lp_towers_sf, dist = 2000) %>%
  st_union()


# ID which tagDeployIDs are within the longpoint tagging location
tagged_at_lp <- flight_lines.pj %>%
  st_drop_geometry() %>%
  distinct(tagDeployID, 
           lon_tagSite, 
           lat_tagSite) %>%
  filter(!is.na(lon_tagSite) & !is.na(lat_tagSite)) %>%
  st_as_sf(coords = c("lon_tagSite", "lat_tagSite"), crs = 4326) %>%
  st_transform(crs = st_crs(flight_lines.pj)) %>%
  filter(lengths(st_intersects(., lp_buffer)) > 0) %>%
  pull(tagDeployID)

# remove those individuals from spring
flight_lines_no_lp <- flight_lines.pj %>%
  filter(!tagDeployID %in% tagged_at_lp)
#flight_lines.pj = 4281
#flight_lines_no_lp = 3051

#remove those individual from fall
fall_flight_lines_no_lp <- fall_flight_lines.pj%>%
  filter(!tagDeployID %in% tagged_at_lp)
#fall flight lines = 4404
#fall flight lines no lp = 4386

##need to remake the flight groups:

#spring nocturnal birds
flight_noLP_noc_bird.pj <- flight_lines_no_lp %>%  
  filter(season == "Spring Migration",
    MigrateTime == "nocturnal",
         Animal == "Bird")

#spring nocturnal takeoff
flight_noLP_takeoff.pj <- flight_lines_no_lp %>%  
  filter(MigrateTime == "nocturnal",
         Animal == "Bird",
         season == "Spring Migration") %>%  
  
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

# spring nocturnal landing

flight_noLP_land.pj <- flight_lines_no_lp %>%  
  filter(
    MigrateTime == "nocturnal",
    Animal == "Bird",
    season == "Spring Migration") %>%  
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

# diurnal birds
flight_noLP_day_bird.pj <- flight_lines_no_lp %>% 
  filter(season == "Spring Migration",
         MigrateTime == "diurnal",
         Animal == "Bird")

#spring bats
flight_noLP_bats.pj <- flight_lines_no_lp %>% 
  filter(season == "Spring Migration",
         Animal == "Bat")

### fall 

# fall nocturnal birds
flight_fall_noLP_noc_bird.pj <- fall_flight_lines_no_lp %>%  
  filter(season == "Fall Migration",
         MigrateTime == "nocturnal",
         Animal == "Bird")

#fall nocturnal takeoff
flight_fall_noLP_takeoff.pj <- fall_flight_lines_no_lp %>%  
  filter(MigrateTime == "nocturnal",
         Animal == "Bird",
         season == "Fall Migration") %>%  
  
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

# Fall nocturnal landing

flight_fall_noLP_land.pj <- fall_flight_lines_no_lp %>%  
  filter(
    MigrateTime == "nocturnal",
    Animal == "Bird",
    season == "Fall Migration") %>%  
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

# fall diurnal birds
flight_fall_noLP_day_bird.pj <- fall_flight_lines_no_lp %>% 
  filter(season == "Fall Migration",
         MigrateTime == "diurnal",
         Animal == "Bird")

#fall bats
flight_fall_noLP_bats.pj <- fall_flight_lines_no_lp %>% 
  filter(season == "Fall Migration",
         Animal == "Bat")



#=======================================================
# creating the group layers pre-prioritizr
# ======================================================

library(dplyr)
library(terra)
library(prioritizr)
library(gurobi)
library(sf)

prioritizrPrep_pipeline <- function(group_sf, 
                                    group_name, 
                                    season = "Spring", 
                                    percentiles = c(0.75, 
                                                    0.50, 
                                                    0.30)) {
  
  message(sprintf("Processing per-species upstream-weighted layers for %s (%s)...", 
                  group_name, 
                  season))
  
  output_dir <- file.path("prioritizrOutput", 
                          "noLongpoint", 
                          season, 
                          group_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load template
  master_template_path <- "masterRegionTemplate.tif"
  if (!file.exists(master_template_path)) {
    stop("Could not find 'masterRegionTemplate.tif'. Run your template creation script first.")
  }
  regionTemplate <- rast(master_template_path)
  
  # Load year-specific station effort file
  effort_rds <- ifelse(grepl("Spring", season, ignore.case = TRUE), 
                       "SpringDepsDuration.rds", 
                       "FallDepsDuration.rds")
  if (!file.exists(effort_rds)) {
    stop(sprintf("Could not find station effort file: '%s'.", effort_rds))
  }
  station_effort_df <- readRDS(effort_rds)
  
  # Extract detection arrival year and join on stationID + year
  group_sf_weighted <- group_sf %>%
    mutate(
      arrival_dt = as.POSIXct(tsEnd_dt),
      year = lubridate::year(arrival_dt)
    ) %>%
    left_join(
      station_effort_df %>% dplyr::select(stationID, year, active_days), 
      by = c("stationID", "year")
    ) %>%
    mutate(
      active_days = ifelse(is.na(active_days) | active_days < 1, 1, active_days),
      station_weight = 1 / log(active_days + 1)
    )
  
  species_list <- unique(group_sf_weighted$species)
  
  # Cumulative Threshold Function (applied per species surface)
  apply_cumulative_threshold <- function(r, top_pct) {
    r_thresh <- r
    vals <- values(r_thresh, mat = FALSE)
    valid_idx <- which(!is.na(vals) & vals > 0)
    
    if (length(valid_idx) > 0) {
      pos_vals <- vals[valid_idx]
      sorted_order <- order(pos_vals, decreasing = TRUE)
      sorted_vals <- pos_vals[sorted_order]
      cum_sum <- cumsum(sorted_vals)
      total_sum <- tail(cum_sum, 1)
      
      cutoff_idx <- which(cum_sum >= (total_sum * top_pct))[1]
      cutoff_val <- sorted_vals[cutoff_idx]
      
      r_thresh[r_thresh < cutoff_val] <- 0
    }
    return(r_thresh)
  }
  
  # Process each percentile tier independently
  results_list <- list()
  
  for (pct in percentiles) {
    pct_label <- paste0("top", pct * 100)
    
    # Loop through each species, process, and threshold it individually
    species_rasters_thresh <- lapply(species_list, 
                                     function(sp_name) {
                                       df_sp <- group_sf_weighted %>% filter(species == sp_name)
                                       
                                       # Rasterize individual species
                                       flightsRast <- rasterize(vect(df_sp), 
                                                                regionTemplate, 
                                                                field = "station_weight", 
                                                                fun = "sum", 
                                                                background = 0)
                                       
                                       # Gaussian smoothing (5 km bandwidth)
                                       weightMatrix <- focalMat(flightsRast,
                                                                d = 5000, type = "Gauss") 
                                       kdeSurface <- focal(flightsRast, w = weightMatrix, fun = sum, na.rm = TRUE)
                                       kdeSurface[is.na(kdeSurface) | is.infinite(kdeSurface)] <- 0
                                       
                                       # Apply cumulative threshold 
                                       thresh_surface <- apply_cumulative_threshold(kdeSurface, pct)
                                       return(thresh_surface)
                                     })
    
    names(species_rasters_thresh) <- species_list
    
    # Stack thresholded species together for this percentile
    thresh_stack <- rast(species_rasters_thresh)
    
    file_path <- file.path(output_dir, 
                           paste0(group_name, "_", 
                                  pct_label, ".tif"))
    writeRaster(thresh_stack, 
                file_path, overwrite = TRUE)
    
    results_list[[pct_label]] <- rast(file_path)
  }
  
  # Create uniform cost layer from region template
  cost_layer <- regionTemplate
  cost_layer[!is.na(cost_layer)] <- 1
  names(cost_layer) <- "cost"
  
  results_list[["cost_layer"]] <- cost_layer
  
  message(sprintf("Finished processing per-species weighted stacks for: %s (%s)!", group_name, season))
  return(results_list)
}

########################################################
### running the above function

### spring

# nocturnal spring birds
spring_noc_noLP <- prioritizrPrep_pipeline(
  group_sf = flight_noLP_noc_bird.pj, 
  group_name = "nocturnal_birds", 
  season = "Spring Migration",
  percentiles = c(0.75, 0.50, 0.30)
)

#spring nocturnal takeoff birds
spring_noc_noLP <- prioritizrPrep_pipeline(
  group_sf = flight_noLP_takeoff.pj, 
  group_name = "nocturnal_takeoff", 
  season = "Spring Migration",
  percentiles = c(0.75, 0.50, 0.30)
)

#spring nocturnal landing birds
spring_land_noLP <- prioritizrPrep_pipeline(
  group_sf = flight_noLP_land.pj, 
  group_name = "nocturnal_land", 
  season = "Spring Migration",
  percentiles = c(0.75, 0.50, 0.30)
)

#spring diurnal birds
spring_day_noLP <- prioritizrPrep_pipeline(
  group_sf = flight_noLP_day_bird.pj, 
  group_name = "day_birds", 
  season = "Spring Migration",
  percentiles = c(0.75, 0.50, 0.30)
)

#spring bats
spring_day_noLP <- prioritizrPrep_pipeline(
  group_sf = flight_noLP_bats.pj, 
  group_name = "bats", 
  season = "Spring Migration",
  percentiles = c(0.75, 0.50, 0.30)
)

### fall

#fall nocturnal birds
fall_noc_noLP <- prioritizrPrep_pipeline(
  group_sf = flight_fall_noLP_noc_bird.pj, 
  group_name = "nocturnal_birds", 
  season = "Fall Migration",
  percentiles = c(0.75, 0.50, 0.30)
)

#fall nocturnal takeoff birds
fall_noc_takeoff <- prioritizrPrep_pipeline(
  group_sf = flight_fall_noLP_takeoff.pj, 
  group_name = "nocturnal_takeoff", 
  season = "Fall Migration",
  percentiles = c(0.75, 0.50, 0.30)
)

#fall nocturnal landing birds
fall_land_noLP <- prioritizrPrep_pipeline(
  group_sf = flight_fall_noLP_land.pj, 
  group_name = "nocturnal_land", 
  season = "Fall Migration",
  percentiles = c(0.75, 0.50, 0.30)
)

#fall diurnal birds
fall_day_noLP <- prioritizrPrep_pipeline(
  group_sf = flight_fall_noLP_day_bird.pj, 
  group_name = "day_birds", 
  season = "Fall Migration",
  percentiles = c(0.75, 0.50, 0.30)
)

#fall bats
fall_bats_noLP <- prioritizrPrep_pipeline(
  group_sf = flight_fall_noLP_bats.pj, 
  group_name = "bats", 
  season = "Fall Migration",
  percentiles = c(0.75, 0.50, 0.30)
)

#========================================================
# prioritizr 
# =======================================================

prioritizr_pipeline <- function(group_name, 
                                season = "Spring", 
                                thresholds = c(0.75, 0.50, 0.30),
                                targets = c(0.17, 0.30, 0.50)) {
  
  # Updated path to include 'effortCorrected'
  output_dir <- file.path("prioritizrOutput", 
                          "noLongpoint",
                          season, 
                          group_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load the pre-saved clean thresholded stacks from the nested folder
  stack_top75 <- rast(file.path(output_dir, 
                                paste0(group_name, 
                                       "_top75.tif")))
  stack_top50 <- rast(file.path(output_dir, 
                                paste0(group_name, 
                                       "_top50.tif")))
  stack_top30 <- rast(file.path(output_dir, 
                                paste0(group_name, 
                                       "_top30.tif")))
  
  # Create uniform cost layer from the loaded stack
  cost_layer <- stack_top75[[1]]
  cost_layer[!is.na(cost_layer)] <- 1
  names(cost_layer) <- "cost"
  
  results_list <- list()
  
  scenarios <- expand.grid(
    threshold = thresholds,
    target = targets
  )
  
  # Run the scenarios
  for (i in 1:nrow(scenarios)) {
    
    curr_thresh <- scenarios$threshold[i]
    curr_target <- scenarios$target[i]
    
    # Select the correct raster stack 
    current_features <- switch(as.character(curr_thresh),
                               "0.75" = stack_top75,
                               "0.5" = stack_top50,  
                               "0.50" = stack_top50,
                               "0.3" = stack_top30,
                               "0.30" = stack_top30,
                               stop("Unknown threshold value encountered!")
    )
    
    # Build the prioritizr problem
    p <- problem(cost_layer, current_features) %>%
      add_min_set_objective() %>%
      add_relative_targets(curr_target) %>%
      add_binary_decisions()
    
    # Set the Gurobi solver 
    p <- p %>% add_gurobi_solver(gap = 0.10, 
                                 time_limit = 600, 
                                 verbose = TRUE)
    
    # Solve the problem
    solution <- solve(p, force = TRUE)
    
    # Gather the outputs
    scenario_output <- list(
      scenario_id = i,
      group = group_name,
      season = season,
      threshold = curr_thresh,
      target = curr_target,
      solution_raster = solution
    )
    
    results_list[[i]] <- scenario_output
    
    # Save wrapped output (clean filename without _corrected)
    scenario_output_wrapped <- scenario_output
    scenario_output_wrapped$solution_raster <- terra::wrap(scenario_output$solution_raster)
    
    save_path <- file.path(output_dir, paste0("scenario_", 
                                              group_name, 
                                              "_id_", 
                                              i, 
                                              ".rds"))
    saveRDS(scenario_output_wrapped, 
            file = save_path)
    
    message(sprintf("Completed & Saved Scenario %d of %d for [%s - %s] -> Thresh: %.2f | Target: %.2f", 
                    i, 
                    nrow(scenarios), 
                    group_name, 
                    season,
                    curr_thresh, 
                    curr_target))
  }
  
  return(results_list)
}

###########################
## running the above function

### spring

#spring nocturnal birds
spring_noc_noLP <- prioritizr_pipeline(
  group_name = "nocturnal_birds", 
  season = "Spring Migration"
)

#spring nocturnal takeoff birds
spring_takeoff_noLP <- prioritizr_pipeline(
  group_name = "nocturnal_takeoff", 
  season = "Spring Migration"
)

#spring nocturnal landing birds
springland_noLP <- prioritizr_pipeline(
  group_name = "nocturnal_land", 
  season = "Spring Migration"
)

#spring dirunal birds
springland_noLP <- prioritizr_pipeline(
  group_name = "day_birds", 
  season = "Spring Migration"
)

#spring bats
springland_noLP <- prioritizr_pipeline(
  group_name = "bats", 
  season = "Spring Migration"
)

### Fall

#Fall nocturnal birds
Fall_noc_noLP <- prioritizr_pipeline(
  group_name = "nocturnal_birds", 
  season = "Fall Migration"
)

#Fall nocturnal takeoff birds
Fall_takeoff_noLP <- prioritizr_pipeline(
  group_name = "nocturnal_takeoff", 
  season = "Fall Migration"
)

#Fall nocturnal landing birds
Fallland_noLP <- prioritizr_pipeline(
  group_name = "nocturnal_land", 
  season = "Fall Migration"
)

#Fall dirunal birds
Fallland_noLP <- prioritizr_pipeline(
  group_name = "day_birds", 
  season = "Fall Migration"
)

#Fall bats
Fallland_noLP <- prioritizr_pipeline(
  group_name = "bats", 
  season = "Fall Migration"
)

#================================================
# evaluating prioritizr outputs (Effort-Corrected)
#================================================

evaluate_scenarios <- function(group_name, 
                               season = "Spring",
                               total_scenarios = 9) {
  
  # Updated file path to include 'effortCorrected'
  output_dir <- file.path("prioritizrOutput",
                          "noLongpoint",
                          season,
                          group_name)
  
  # Load thresholded raster stacks from the effortCorrected folder
  stack_top75 <- rast(file.path(output_dir, paste0(group_name, "_top75.tif")))
  stack_top50 <- rast(file.path(output_dir, paste0(group_name, "_top50.tif")))
  stack_top30 <- rast(file.path(output_dir, paste0(group_name, "_top30.tif")))
  
  # Creating cost layer
  cost_layer <- stack_top75[[1]]
  cost_layer[!is.na(cost_layer)] <- 1
  names(cost_layer) <- "cost"
  
  summary_table <- data.frame() 
  
  # Loop through the scenarios
  for (i in 1:total_scenarios) {
    
    file_path <- file.path(output_dir, 
                           paste0("scenario_", 
                                  group_name, 
                                  "_id_", 
                                  i, 
                                  ".rds"))
    
    if (!file.exists(file_path)) {
      next # Skip if the file doesn't exist
    }
    
    message(sprintf("Evaluating Scenario ID %d for [%s - %s]...", 
                    i, 
                    group_name,
                    season))
    
    # Load and unwrap the saved scenario
    wrapped_obj <- readRDS(file_path)
    sol_raster <- terra::unwrap(wrapped_obj$solution_raster)
    
    # Select the correct feature stack based on the saved threshold
    curr_thresh <- wrapped_obj$threshold
    current_features <- switch(as.character(curr_thresh),
                               "0.75" = stack_top75,
                               "0.5" = stack_top50,  
                               "0.50" = stack_top50,
                               "0.3" = stack_top30,
                               "0.30" = stack_top30,
                               stop(sprintf("Unknown threshold value in file: %s", file_path))
    )
    
    # Re-create the problem object
    p <- problem(cost_layer, 
                 current_features) %>%
      add_min_set_objective() %>%
      add_relative_targets(wrapped_obj$target) %>%
      add_binary_decisions()
    
    # Run evaluation summary functions
    cost_sum <- eval_cost_summary(p, sol_raster)
    n_sum <- eval_n_summary(p, sol_raster)
    bound_sum <- eval_boundary_summary(p, sol_raster)
    target_sum <- eval_target_coverage_summary(p, sol_raster)
    
    # Compile metrics
    row_df <- data.frame(
      scenario_id = wrapped_obj$scenario_id,
      threshold = curr_thresh,
      target = wrapped_obj$target,
      total_cost = cost_sum$cost,
      num_units = n_sum$n,
      boundary_len = bound_sum$boundary,
      all_targets_met = all(target_sum$met),
      min_held = min(target_sum$relative_held)
    )
    
    summary_table <- rbind(summary_table, row_df)
  }
  
  # Save the evaluation summary as a csv file inside the group folder
  output_csv_name <- file.path(output_dir, paste0(group_name, "_evaluation.csv"))
  write.csv(summary_table, output_csv_name, row.names = FALSE)
  message(sprintf("Successfully generated evaluation table for [%s - %s]. Saved to '%s'", group_name, season, output_csv_name))
  
  return(summary_table)
}

#############################
##running the function

## spring

#spring nocturnal birds
spring_eval <- evaluate_scenarios(group_name = "nocturnal_birds", 
                                  season = "Spring Migration")

#spring nocturnal takeoff birds
takeoff_eval <- evaluate_scenarios(group_name = "nocturnal_takeoff", 
                                  season = "Spring Migration")

#spring nocturnal landing birds
land_eval <- evaluate_scenarios(group_name = "nocturnal_land", 
                                   season = "Spring Migration")

#spring diurnal birds
day_eval <- evaluate_scenarios(group_name = "day_birds", 
                                season = "Spring Migration")

#spring bats
bat_eval <- evaluate_scenarios(group_name = "bats", 
                               season = "Spring Migration")

### fall

#Fall nocturnal birds
Fall_eval <- evaluate_scenarios(group_name = "nocturnal_birds", 
                                  season = "Fall Migration")

#Fall nocturnal takeoff birds
takeoff_eval <- evaluate_scenarios(group_name = "nocturnal_takeoff", 
                                   season = "Fall Migration")

#Fall nocturnal landing birds
land_eval <- evaluate_scenarios(group_name = "nocturnal_land", 
                                season = "Fall Migration")

#Fall diurnal birds
day_eval <- evaluate_scenarios(group_name = "day_birds", 
                               season = "Fall Migration")

#Fall bats
bat_eval <- evaluate_scenarios(group_name = "bats", 
                               season = "Fall Migration")
