#script for running prioritizr analysis for migratory birds

library(dplyr)
library(terra)
library(prioritizr)
library(gurobi)

#=========================================
# preparing the all birds data layers
#=========================================

# features layer (species line KDE for birds)
allBirds <- rast("LKDEBirdRasterStack.tif")

#applying a threshold to remove very small values from being selected
apply_threshold <- function(rast_stack, top_pct) {
  thresh_stack <- rast_stack
  for (i in 1:nlyr(thresh_stack)) {
    r <- thresh_stack[[i]]
    
    # Extract values 
    vals <- values(r, mat = FALSE)
    vals_positive <- vals[vals > 0]
    
    if (length(vals_positive) > 0) {
      cutoff <- quantile(vals_positive, 
                         probs = (1 - top_pct), 
                         names = FALSE, 
                         na.rm = TRUE)
      
      # Modify raster values to be above the threshold
      r[r < cutoff] <- 0
      thresh_stack[[i]] <- r
    }
    
    # Explicitly clear temporary objects to free RAM for the next loop
    rm(vals, vals_positive)
    gc(verbose = FALSE)
  }
  return(thresh_stack)
}
# Generate the thresholded stacks
stack_top75 <- apply_threshold(allBirds, 0.75)
stack_top50 <- apply_threshold(allBirds, 0.50)
stack_top30 <- apply_threshold(allBirds, 0.30)

writeRaster(stack_top75, "stack_top75.tif", overwrite = TRUE)
writeRaster(stack_top50, "stack_top50.tif", overwrite = TRUE)
writeRaster(stack_top30, "stack_top30.tif", overwrite = TRUE)

stack_top75 <- rast("stack_top75.tif")
stack_top50 <- rast("stack_top50.tif")
stack_top30 <- rast("stack_top30.tif")

# create cost layer with even cost
cost_layer <- allBirds[[1]]
cost_layer[!is.na(cost_layer)] <- 1
names(cost_layer) <- "cost"

#=======================================================
# creating the group layers pre-prioritizr
# ======================================================

regionTemplate <- rast("LKDEBirdRasterStack.tif")[[1]]
regionTemplate[!is.na(regionTemplate)] <- 0

### turning the data prep from above into a function for rerunning

prioritizrPrep_pipeline <- function(group_sf, 
                                    regionTemplate, 
                                    group_name, 
                                    percentiles = c(0.75, 
                                                    0.50, 
                                                    0.30)) {
  
  message(sprintf("Processing spatial layers for: %s...", 
                  group_name))
  
  # create file path
  output_dir <- file.path("prioritizrOutput", 
                          group_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get the species list 
  species_list <- unique(group_sf$species)
  
  # Loop through species to create rasters and apply Gaussian smoothing (KDE)
  species_rasters <- lapply(species_list, 
                            function(sp_name) {
    
    df_sp <- group_sf %>% filter(species == sp_name)
    
    # Rasterize counts
    flightsRast <- rasterize(vect(df_sp), 
                             regionTemplate, 
                             field = 1, 
                             fun = "sum", 
                             background = 0)
    
    # Gaussian smoothing kernel (5 km bandwidth)
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
  
  # stack species rasters
  all_stack <- rast(species_rasters)
  
  # thresholding function 
  apply_threshold <- function(rast_stack, top_pct) {
    thresh_stack <- rast_stack
    for (i in 1:nlyr(thresh_stack)) {
      r <- thresh_stack[[i]]
      vals <- values(r, mat = FALSE)
      vals_positive <- vals[vals > 0]
      
      if (length(vals_positive) > 0) {
        cutoff <- quantile(vals_positive, 
                           probs = (1 - top_pct), 
                           names = FALSE, 
                           na.rm = TRUE)
        r[r < cutoff] <- 0
        thresh_stack[[i]] <- r
      }
      rm(vals, vals_positive)
      gc(verbose = FALSE)
    }
    return(thresh_stack)
  }
  
  # Make and save thresholded stacks 
  results_list <- list()
  for (pct in percentiles) {
    pct_label <- paste0("top", pct * 100)
    thresh_stack <- apply_threshold(all_stack,
                                    pct)
    
    file_path <- file.path(output_dir, paste0(group_name, 
                                              "_", 
                                              pct_label, 
                                              ".tif"))
    writeRaster(thresh_stack, 
                file_path, 
                overwrite = TRUE)
    
    results_list[[pct_label]] <- rast(file_path)
  }
  
  # Create uniform cost layer from the first layer of stack
  cost_layer <- all_stack[[1]]
  cost_layer[!is.na(cost_layer)] <- 1
  names(cost_layer) <- "cost"
  
  results_list[["cost_layer"]] <- cost_layer
  results_list[["raw_stack"]] <- all_stack
  
  message(sprintf("Finished processing and saved files for: %s!", group_name))
  return(results_list)
}
### applying function to each group

#all nocturnal birds
nocturnal_results <- prioritizrPrep_pipeline(flight_noc_bird.pj, 
                                            regionTemplate, 
                                            group_name = "nocturnal_birds")

#nocturnal birds takeoff
flight_noc_takeoff.pj<- st_read("flight_noc_takeoff.pj.gpkg")
nocTakeOff_results <- prioritizrPrep_pipeline(flight_noc_takeoff.pj, 
                                             regionTemplate, 
                                             group_name = "nocturnal_takeoff")

#nocturnal birds landing
flight_noc_land.pj<- st_read("flight_noc_land.pj.gpkg")
nocLand_results <- prioritizrPrep_pipeline(flight_noc_land.pj, 
                                              regionTemplate, 
                                              group_name = "nocturnal_land")
#diurnal birds all
flight_day_birds.pj <- st_read("day_birds.pj.gpkg")
daybird_results <-  prioritizrPrep_pipeline(flight_day_birds.pj,
                                            regionTemplate,
                                            group_name = "day_birds")

# all bats
all_bats.pj <- st_read("all_bats.pj.gpkg")
bats_results <-  prioritizrPrep_pipeline(all_bats.pj,
                                         regionTemplate,
                                         group_name = "all_bats")

#========================================================
# prioritizr
# =======================================================

prioritizr_pipeline <- function(group_name, 
                                thresholds = c(0.75, 0.50, 0.30),
                                targets = c(0.17, 0.30, 0.50)) {
  
  # get output path
  output_dir <- file.path("prioritizrOutput", group_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  #  Load the pre-saved thresholded stacks from the group's nested directory
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
      threshold = curr_thresh,
      target = curr_target,
      solution_raster = solution
    )
    
    results_list[[i]] <- scenario_output
    
    # save wrapped output
    scenario_output_wrapped <- scenario_output
    scenario_output_wrapped$solution_raster <- terra::wrap(scenario_output$solution_raster)
    
    save_path <- file.path(output_dir, paste0("scenario_", 
                                              group_name, 
                                              "_id_", 
                                              i, 
                                              ".rds"))
    saveRDS(scenario_output_wrapped, 
            file = save_path)
    
    message(sprintf("Completed & Saved Scenario %d of %d for [%s] -> Thresh: %.2f | Target: %.2f", 
                    i, 
                    nrow(scenarios), 
                    group_name, 
                    curr_thresh, 
                    curr_target))
  }
  
  return(results_list)
}

### running the function above ###

#nocturnal birds

nocturnal_bird_results <- prioritizr_pipeline(group_name = "nocturnal_birds")

#nocturnal birds takeoff

nocturnal_bird_takeoff <- prioritizr_pipeline(group_name = "nocturnal_takeoff")

#nocturnal birds landing

nocturnal_bird_land <- prioritizr_pipeline(group_name = "nocturnal_land")

#diurnal birds all

diurnal_bird_results <- prioritizr_pipeline(group_name = "day_birds")

#all bats

all_bats_results <- prioritizr_pipeline(group_name = "all_bats")

#================================================
# evaluating prioritizr outputs
#================================================

evaluate_scenarios <- function(group_name, 
                               total_scenarios = 9) {
  
  # file path
  output_dir <- file.path("prioritizrOutput", group_name)
  
  # load thresholded rasterstacks
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
    
    message(sprintf("Evaluating Scenario ID %d for [%s]...", 
                    i, 
                    group_name))
    
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
  message(sprintf("Successfully generated evaluation table for [%s]. Saved to '%s'", group_name, output_csv_name))
  
  return(summary_table)
}

# evaluating results
noc_birds_evaluate <- evaluate_scenarios("nocturnal_birds")
noc_takeoff_evaluate <- evaluate_scenarios("nocturnal_takeoff")
noc_land_evaluate <- evaluate_scenarios("nocturnal_land")
day_birds_evaluate <- evaluate_scenarios("day_birds")
all_bats <- evaluate_scenarios("all_bats")

#=============================================================
#plotting the prioritizr results
#=============================================================

library(terra)
library(sf)
library(maptiles)
library(ggplot2)
library(tidyterra)

plot_group_scenario_maps <- function(group_name, 
                                     total_scenarios = 9) {
  
  # creating the output paths
  output_dir <- file.path("prioritizrOutput", group_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(sprintf("Created nested directory: '%s'", output_dir))
  }
  
  # Set up the mapping extent and reproject
  smaller_extent <- st_as_sfc(st_bbox(
    c(xmin = -94.46, 
      xmax = -74.24, 
      ymin = 40, 
      ymax = 48), 
    crs = 4326 # mapview default CRS
  ))
  
  smaller_extent.pj <- st_transform(smaller_extent, 
                                    crs = 3978)
  
  # get basemap tiles once for the region
  tiles <- get_tiles(
    x = smaller_extent.pj, 
    provider = "CartoDB.PositronNoLabels", 
    zoom = 8, 
    crop = TRUE
  )
  
  b <- st_bbox(tiles)
  
  # Loop through scenarios
  for (i in 1:total_scenarios) {
    
    # Look for the .rds file 
    file_path <- file.path(output_dir, 
                           paste0("scenario_", 
                                  group_name, 
                                  "_id_", 
                                  i, 
                                  ".rds"))
    
    if (!file.exists(file_path)) {
      next # Skip if the file doesn't exist in that folder
    }
    
    message(sprintf("Processing map for Scenario ID %d [%s]...", i, group_name))
    
    # Load and unwrap data
    data <- readRDS(file_path)
    rast <- terra::unwrap(data$solution_raster)
    
    # Extract metadata for labeling
    scenario_name <- if (!is.null(data$name)) data$name else paste0(group_name, 
                                                                    "_id_",
                                                                    i)
    threshold <- data$threshold
    target <- data$target
    penalty <- ifelse(is.null(data$penalty), 
                      0, 
                      data$penalty)
    
    # Crop raster to match map tiles extent
    rast_cropped <- terra::crop(rast, 
                                ext(b["xmin"], 
                                    b["xmax"], 
                                    b["ymin"], 
                                    b["ymax"]))
    
    # Generate ggplot
    p <- ggplot() +
      geom_spatraster_rgb(data = tiles, maxcell = 500000) +
      geom_spatraster(data = as.factor(rast_cropped), maxcell = 500000) + 
      scale_fill_manual(values = c("0" = "#ffffff00", "1" = "#2E933C")) +
      labs(
        title = paste("Scenario:", scenario_name),
        subtitle = paste("Threshold:", threshold, 
                         "| Target:", target, 
                         "| Penalty:", penalty)
      ) +
      theme_minimal() +
      coord_sf(
        xlim = c(b["xmin"], b["xmax"]), 
        ylim = c(b["ymin"], b["ymax"]), 
        expand = FALSE
      )
    
    # Save PDF directly into the nested group folder
    pdf_filename <- file.path(output_dir, paste0("Map_scenario_id_", i, ".pdf"))
    ggsave(pdf_filename, plot = p, width = 8, height = 6)
  }
  
  message(sprintf("All maps for [%s] successfully saved to '%s/'", group_name, output_dir))
}

nocturnal_bird_map <- plot_group_scenario_maps("nocturnal_birds")
nocturnal_takeoff_map <- plot_group_scenario_maps("nocturnal_takeoff")
nocturnal_land_map <- plot_group_scenario_maps("nocturnal_land")
day_birds_map <- plot_group_scenario_maps("day_birds")
all_bats_map <- plot_group_scenario_maps("all_bats")

#=========================================================
# Selection frequency plot
#=========================================================
library(terra)
library(sf)
library(maptiles)
library(ggplot2)
library(tidyterra)

#create the station file for plotting later
stations <- read.csv("StationDownloads/Motus-data-region-undefined_stations_downloaded-2026-09-09.csv")
stations.sf <- st_as_sf(stations,
                        coords = c("longitude", 
                                   "latitude"),
                        crs = 4326)
stations_pj <- st_transform(stations.sf, crs = 3978)

GLWatershed <- st_read("./greatlakes_subbasins/greatlakes_subbasins.shp", 
                       quiet = TRUE)
GLWS_proj2 <- st_transform(GLWatershed, 
                           crs = 3978)

# Filter points outside of the watershed
stations_pj2 <- stations_pj %>%
  st_filter(GLWS_proj2)

## selection frequency plot function 

plot_selection_frequency <- function(group_name) {
  
  # define file path
  output_dir <- file.path("prioritizrOutput", 
                          group_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, 
               recursive = TRUE)
  }
  
  # mapping extent and basemap
  smaller_extent <- st_as_sfc(st_bbox(
    c(xmin = -94.46, 
      xmax = -74.24, 
      ymin = 40, 
      ymax = 48), 
    crs = 4326
  ))
  smaller_extent.pj <- st_transform(smaller_extent, 
                                    crs = 3978)
  
  tiles <- get_tiles(
    x = smaller_extent.pj, 
    provider = "Esri.WorldGrayCanvas", 
    zoom = 8, 
    crop = TRUE
  )
  b <- st_bbox(tiles)
  
  # load the files
  scenario_files <- list.files(path = output_dir, 
                               pattern = "^scenario_.*\\.rds$", 
                               full.names = TRUE)
  
  if (length(scenario_files) == 0) {
    stop(sprintf("No scenario files found in directory: '%s'", 
                 output_dir))
  }
  
  message(sprintf("Processing %d scenarios for [%s]...", 
                  length(scenario_files), 
                  group_name))
  
  # use the first file as a template for the summed freq file
  first_data <- readRDS(scenario_files[1])
  sum_raster <- terra::unwrap(first_data$solution_raster)
  
  # Loop through remaining files and add them together
  if (length(scenario_files) > 1) {
    for (i in 2:length(scenario_files)) {
      data <- readRDS(scenario_files[i])
      r <- terra::unwrap(data$solution_raster)
      sum_raster <- sum_raster + r
    }
  }
  
  # Convert sum to proportion(0 to 1)
  freq_raster <- sum_raster / length(scenario_files)
  
  # Crop to smaller extent
  freq_cropped <- terra::crop(freq_raster, 
                              ext(b["xmin"], 
                                  b["xmax"], 
                                  b["ymin"], 
                                  b["ymax"]))

  
  # plot
  p_freq <- ggplot() +
    geom_spatraster_rgb(data = tiles, 
                        maxcell = 500000) +
    geom_spatraster(data = freq_cropped, 
                    #alpha = 0.75,
                    na.rm = TRUE,
                    maxcell = 500000) + 
    scale_fill_viridis_c(
      option = "plasma", 
      name = "Selection\nFrequency",
      na.value = "transparent",
      limits = c(0.01, 1)
    ) +
    geom_sf(
      data = stations_pj2, 
      color = "black",
      fill = NA,
      size = 2,
      stroke = 0.3,
      shape = 21
    ) +
    labs(
      title = paste("Selection Frequency Across Scenarios:", 
                    group_name)
    ) +
    theme_minimal() +
    coord_sf(
      xlim = c(b["xmin"], b["xmax"]), 
      ylim = c(b["ymin"], b["ymax"]), 
      expand = FALSE
    )
  
  # Save PDF into the group's folder
  pdf_filename <- file.path(output_dir,
                            paste0("Selection_Frequency_Map_", 
                                   group_name, ".pdf"))
  ggsave(pdf_filename, 
         plot = p_freq, 
         width = 8, 
         height = 6)
  
  message(sprintf("Selection frequency map successfully saved to '%s'", 
                  pdf_filename))
  return(p_freq)
}

noc_birds_selection <- plot_selection_frequency("nocturnal_birds")
noc_take_selection <- plot_selection_frequency("nocturnal_takeoff")
noc_land_selection <- plot_selection_frequency("nocturnal_land")
day_bird_selection <- plot_selection_frequency("day_birds")
all_bats_selection <- plot_selection_frequency("all_bats")

#================================================================
#irreplaceability map
#================================================================

library(terra)
library(sf)
library(maptiles)
library(ggplot2)
library(tidyterra)
library(prioritizr)

#I want to plot scenario with 50% target, 75% threshold
# to find all scenarios

all_scenarios <- dplyr::bind_rows(lapply(
  list.files(path = file.path("prioritizrOutput", "nocturnal_birds"), 
             pattern = "^scenario_.*\\.rds$", 
             full.names = TRUE),
  function(f) {
    d <- readRDS(f)
  data.frame(
    file = f,
    target = if (!is.null(d$target)) d$target else NA,
    threshold = if (!is.null(d$threshold)) d$threshold else NA,
    penalty = if (!is.null(d$penalty)) d$penalty else NA
  )
}))

print(all_scenarios) #scenario 7 matches

### function to make irreplaceabiltiy maps for each group

plot_group_irreplaceability_7 <- function(group_name, 
                                          raster_stack_path = NULL) {
  
  message(sprintf("Processing irreplaceability map for: %s...", group_name))
  
  # define the path
  output_dir <- file.path("prioritizrOutput", group_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load Scenario 7 for the group
  scenario_path <- file.path(output_dir, 
                             paste0("scenario_", 
                                    group_name, 
                                    "_id_7.rds"))
  if (!file.exists(scenario_path)) {
    stop(sprintf("Could not find file: '%s'", 
                 scenario_path))
  }
  irreplace7 <- readRDS(scenario_path)
  
  # find raster stack path
  if (is.null(raster_stack_path)) {
    pct_label <- paste0("top", 
                        irreplace7$threshold * 100)
    raster_stack_path <- file.path(output_dir, 
                                   paste0(group_name, 
                                          "_", 
                                          pct_label, 
                                          ".tif"))
  }
  
  if (!file.exists(raster_stack_path)) {
    stop(sprintf("Could not find matching raster stack file at: '%s'.", 
                 raster_stack_path))
  }
  
  # Set mapping extent and basemap tiles
  smaller_extent <- st_as_sfc(st_bbox(
    c(xmin = -94.46, 
      xmax = -74.24, 
      ymin = 40, 
      ymax = 48), 
    crs = 4326
  ))
  smaller_extent.pj <- st_transform(smaller_extent, 
                                    crs = 3978)
  
  tiles <- get_tiles(
    x = smaller_extent.pj, 
    provider = "Esri.WorldGrayCanvas", 
    zoom = 8, 
    crop = TRUE
  )
  b <- st_bbox(tiles)
  
  #redo the prioritizr problem
  group_rast <- rast(raster_stack_path)
  
  cost_layer <- group_rast[[1]]
  cost_layer[!is.na(cost_layer)] <- 1
  names(cost_layer) <- "cost"
  
  p_problem <- problem(cost_layer, 
                       group_rast) %>%
    add_min_set_objective() %>%
    add_relative_targets(irreplace7$target) %>%
    add_binary_decisions()
  
  # Extract solution
  p_solution <- terra::unwrap(irreplace7$solution_raster)
  
  # Calculate rank importance scores
  message(sprintf("Calculating rank importance scores (n = 10)..."))
  importance_scores <- eval_rank_importance(
    p_problem, 
    p_solution, 
    n = 10,
    force = TRUE
  )
  
  # make the plot
  p_rank <- ggplot() +
    geom_spatraster_rgb(data = tiles, 
                        maxcell = 500000) +
    geom_spatraster(data = importance_scores, 
                    maxcell = 500000) + 
    scale_fill_viridis_c(
      option = "plasma", 
      name = "Importance Rank",
      na.value = "transparent",
      n.breaks = 10
    ) +
    geom_sf(
      data = stations_pj2,
      color = "black",
      fill = NA,
      size = 2,
      shape = 21
    ) +
    labs(
      title = paste("Irreplaceability Rank Importance Map:", group_name),
      subtitle = paste("Scenario 7 | Target:", irreplace7$target, "| Threshold:", irreplace7$threshold)
    ) +
    theme_minimal() +
    coord_sf(
      xlim = c(b["xmin"], b["xmax"]), 
      ylim = c(b["ymin"], b["ymax"]), 
      expand = FALSE
    )
  
  # Save PDF 
  pdf_filename <- file.path(output_dir, "Scenario_7_Rank_Importance_Map.pdf")
  ggsave(pdf_filename, plot = p_rank, width = 8, height = 6)
  
  message(sprintf("Irreplaceability map successfully saved to '%s'", pdf_filename))
  return(p_rank)
}

noc_bird_irr <- plot_group_irreplaceability_7("nocturnal_birds")
noc_take_irr <- plot_group_irreplaceability_7("nocturnal_takeoff")
noc_land_irr <- plot_group_irreplaceability_7("nocturnal_land")
day_bird_irr <- plot_group_irreplaceability_7("day_birds")
all_bats_irr <- plot_group_irreplaceability_7("all_bats")
