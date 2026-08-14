#script for running prioritizr analysis for migratory birds

library(dplyr)
library(terra)
library(prioritizr)
library(highs)

#=========================================
# preparing the data layers
#=========================================

# features layer (species line KDE)
allSpecies <- rast("LKDERasterStack.tif")

apply_threshold <- function(rast_stack, top_pct) {
  thresh_stack <- rast_stack
  for (i in 1:nlyr(thresh_stack)) {
    r <- thresh_stack[[i]]
    
    # Extract values directly as a vector without extra baggage
    vals <- values(r, mat = FALSE)
    vals_positive <- vals[vals > 0]
    
    if (length(vals_positive) > 0) {
      cutoff <- quantile(vals_positive, probs = (1 - top_pct), names = FALSE, na.rm = TRUE)
      
      # Modify raster values efficiently using terra's internal indexing
      r[r < cutoff] <- 0
      thresh_stack[[i]] <- r
    }
    
    # Explicitly clear temporary objects to free RAM for the next loop
    rm(vals, vals_positive)
    gc(verbose = FALSE)
  }
  return(thresh_stack)
}
# Generate your thresholded stacks
stack_top75 <- apply_threshold(allSpecies, 0.75)
stack_top50 <- apply_threshold(allSpecies, 0.50)
stack_top30 <- apply_threshold(allSpecies, 0.30)

writeRaster(stack_top75, "stack_top75.tif", overwrite = TRUE)
writeRaster(stack_top50, "stack_top50.tif", overwrite = TRUE)
writeRaster(stack_top30, "stack_top30.tif", overwrite = TRUE)

stack_top75 <- rast("stack_top75.tif")
stack_top50 <- rast("stack_top50.tif")
stack_top30 <- rast("stack_top30.tif")

# create cost layer with even cost
cost_layer <- allSpecies[[1]]
cost_layer[!is.na(cost_layer)] <- 1
names(cost_layer) <- "cost"

#========================================================
# prioritizr
# =======================================================

results_list <- list()

thresholds <- c(0.75, 0.50, 0.30)
targets <- c(0.17, 0.30, 0.50)

scenarios <- expand.grid(
  threshold = thresholds,
  target = targets
)

# Run the scenarios (without boundary constraints for now)
for (i in 1:nrow(scenarios)) {
  
  curr_thresh <- scenarios$threshold[i]
  curr_target <- scenarios$target[i]
  
  # Select the correct raster stack 
  current_features <- switch(as.character(curr_thresh),
                             "0.75" = stack_top75,
                             "0.5"  = stack_top50,  
                             "0.50" = stack_top50,
                             "0.3"  = stack_top30,
                             "0.30" = stack_top30,
                             stop("Unknown threshold value encountered!")
  )
  
  # Build the prioritizr problem
  p <- problem(cost_layer, current_features) %>%
    add_min_set_objective() %>%
    add_relative_targets(curr_target) %>%
    add_binary_decisions()
  
  # Set the solver (HiGHS)
  p <- p %>% add_highs_solver(gap = 0.10, 
                              time_limit = 600, 
                              verbose = TRUE)
  
  # Solve the problem
  solution <- solve(p, force = TRUE)
  
  # Bundle output
  scenario_output <- list(
    scenario_id = i,
    threshold = curr_thresh,
    target = curr_target,
    solution_raster = solution
  )
  
  results_list[[i]] <- scenario_output
  assign(paste0("scenario_", i), scenario_output)
  
  # Save file 
  scenario_output_wrapped <- scenario_output
  scenario_output_wrapped$solution_raster <- terra::wrap(scenario_output$solution_raster)
  saveRDS(scenario_output_wrapped, file = paste0("scenario_", i, ".rds"))
  
  message(sprintf("Completed & Saved Scenario %d of %d -> Thresh: %.2f | Target: %.2f", 
                  i, nrow(scenarios), curr_thresh, curr_target))
}


## trying to add in the boundary constraint with highs solver

constrained_results_list <- list()

thresholds <- c(0.75, 0.50, 0.30)
targets <- c(0.17, 0.30, 0.50)
penalties <- c(0.001, 0.05) # Only the new boundary penalties

# Create the grid for just the constrained scenarios
scenarios_constrained <- expand.grid(
  threshold = thresholds,
  target = targets,
  penalty = penalties
)

# Run the loop for the constrained scenarios
for (i in 1:nrow(scenarios_constrained)) {
  
  curr_thresh <- scenarios_constrained$threshold[i]
  curr_target <- scenarios_constrained$target[i]
  curr_penalty <- scenarios_constrained$penalty[i]
  
  # Select the correct raster stack safely using character conversion
  current_features <- switch(as.character(curr_thresh),
                             "0.75" = stack_top75,
                             "0.5"  = stack_top50,
                             "0.50" = stack_top50,
                             "0.3"  = stack_top30,
                             "0.30" = stack_top30,
                             stop("Unknown threshold value encountered!")
  )
  
  # Build the prioritizr problem with boundary penalties using knapsack
  p <- problem(cost_layer, current_features) %>%
    add_min_set_objective() %>%
    add_relative_targets(curr_target) %>%
    add_binary_decisions() %>%
    add_boundary_penalties(penalty = curr_penalty, 
                           edge_factor = 0.5, 
                           formulation = "knapsack") %>%
    add_highs_solver(gap = 0.10, 
                     time_limit = 600, 
                     verbose = TRUE)
  
  # Solve the problem
  solution <- solve(p, force = TRUE)
  
  # Bundle output (offsetting index name so it doesn't overwrite 1-9)
  scenario_id_offset <- i + 9 
  
  scenario_output <- list(
    scenario_id = scenario_id_offset,
    threshold = curr_thresh,
    target = curr_target,
    penalty = curr_penalty,
    solution_raster = solution
  )
  
  constrained_results_list[[i]] <- scenario_output
  assign(paste0("scenario_", scenario_id_offset), scenario_output)
  
  # Save file 
  scenario_output_wrapped <- scenario_output
  scenario_output_wrapped$solution_raster <- terra::wrap(scenario_output$solution_raster)
  saveRDS(scenario_output_wrapped, file = paste0("scenario_", scenario_id_offset, ".rds"))
  
  message(sprintf("Completed & Saved Constrained Scenario %d (ID %d) -> Thresh: %.2f | Target: %.2f | Penalty: %.3f", 
                  i, scenario_id_offset, curr_thresh, curr_target, curr_penalty))
}


#================================================
# evaluating prioritizr outputs
#================================================

# Initialize an empty data frame to store results
summary_table <- data.frame()
total_scenarios <- 27 

for (i in 1:total_scenarios) {
  
  file_path <- paste0("scenario_", i, ".rds")
  if (!file.exists(file_path)) next
  
  #load the saved scenarios
  wrapped_obj <- readRDS(file_path)
  sol_raster <- terra::unwrap(wrapped_obj$solution_raster)
  
  # Select the correct feature stack based on the saved threshold
  curr_thresh <- wrapped_obj$threshold
  current_features <- switch(as.character(curr_thresh),
                             "0.75" = stack_top75,
                             "0.5"  = stack_top50,
                             "0.50" = stack_top50,
                             "0.3"  = stack_top30,
                             "0.30" = stack_top30,
                             stop("Unknown threshold value in saved file!")
  )
  
  # Re-create the base problem object
  p <- problem(cost_layer, current_features) %>%
    add_min_set_objective() %>%
    add_relative_targets(wrapped_obj$target) %>%
    add_binary_decisions()
  
  #add boundary penalties only if they were used in this scenario
  if (!is.null(wrapped_obj$penalty) && !is.na(wrapped_obj$penalty) && wrapped_obj$penalty > 0) {
    p <- p %>% add_boundary_penalties(penalty = wrapped_obj$penalty, 
                                      edge_factor = 0.5)
  }
  
  # Run evaluation functions
  cost_sum  <- eval_cost_summary(p, sol_raster)
  n_sum <- eval_n_summary(p, sol_raster)
  bound_sum <- eval_boundary_summary(p, sol_raster)
  target_sum <- eval_target_coverage_summary(p, sol_raster)
  
  # combines the evaluation metrics
  summary_table <- rbind(summary_table, data.frame(
    scenario_id = wrapped_obj$scenario_id,
    threshold = curr_thresh,
    target = wrapped_obj$target,
    penalty = ifelse(is.null(wrapped_obj$penalty), 0, wrapped_obj$penalty),
    total_cost = cost_sum$cost,
    num_units = n_sum$n,
    boundary_len = bound_sum$boundary,
    all_targets_met = all(target_sum$met),
    min_held = min(target_sum$relative_held)
  ))
}

# View the final comparison table
print(summary_table)
write.csv(summary_table, "allFlightsPrioritizrEvaluation.csv")

#=============================================================
#plotting the prioritizr results
#=============================================================

library(terra)
library(sf)
library(maptiles)
library(ggplot2)
library(tidyterra)

#make a list of the scenarios
scenario_files <- paste0("scenario_", 1:27, ".rds")

#make a smaller extent for mapping (got the values from mapview)
smaller_extent <- st_as_sfc(st_bbox(
  c(xmin = -94.46, 
    xmax = -74.24, 
    ymin = 40, 
    ymax = 48), 
  crs = 4326 #mapview default CRS
))
#reproject it 
smaller_extent.pj <- st_transform(smaller_extent, 
                                  crs = 3978)
#for the basemap
tiles <- get_tiles(
  x = smaller_extent.pj, 
  provider = "CartoDB.PositronNoLabels", 
  zoom = 8, 
  crop = TRUE
)

b <- st_bbox(tiles)

#loop through the scenarios

for (file in scenario_files) {
  # Load data
  data <- readRDS(file)
  rast <- terra::unwrap(data$solution_raster)
  
  # get the threshold, targets, penalties for naming
  scenario_name <- if (!is.null(data$name)) data$name else gsub("\\.rds$", 
                                                                "", 
                                                                file)
  threshold <- data$threshold
  target <- data$target
  penalty <- data$penalty
  
  # Crop
  rast_cropped <- terra::crop(rast, 
                              ext(b["xmin"], 
                                  b["xmax"], 
                                  b["ymin"], 
                                  b["ymax"]))
  
  # Plot
  p <- ggplot() +
    geom_spatraster_rgb(data = tiles, maxcell = 500000) +
    geom_spatraster(data = as.factor(rast_cropped), maxcell = 500000) + 
    scale_fill_manual(values = c("0" = "#ffffff00", "1" = "#2E933C")) +
    labs(
      title = paste("Scenario:", 
                    scenario_name),
      subtitle = paste("Threshold:", threshold, 
                       "| Target:", target, 
                       "| Penalty:", penalty)
    ) +
    theme_minimal() +
    coord_sf(
      xlim = c(b["xmin"],
               b["xmax"]), 
      ylim = c(b["ymin"], 
               b["ymax"]), 
      expand = FALSE
    )
  
  # Save
  filename <- paste0("Map_", scenario_name, ".pdf")
  ggsave(filename, plot = p, width = 8, height = 6)
}


#=========================================================
# Selection frequency plot
#=========================================================
library(terra)
library(sf)
library(maptiles)
library(ggplot2)
library(tidyterra)

# get all the scenarios
scenario_files <- list.files(pattern = "^scenario_.*\\.rds$")

#need to remake tiles
tiles <- get_tiles(
  x = smaller_extent.pj, #same as prev section
  provider = "CartoDB.PositronNoLabels", 
  zoom = 8, 
  crop = TRUE
)
b <- st_bbox(tiles)

# going to loop through but start with 1 to get a template
first_data <- readRDS(scenario_files[1])
sum_raster <- terra::unwrap(first_data$solution_raster)

# Loop through the remaining files and add them together
for (i in 2:length(scenario_files)) {
  data <- readRDS(scenario_files[i])
  r <- terra::unwrap(data$solution_raster)
  sum_raster <- sum_raster + r
}

# Convert the sum to selection frequency (proportion from 0 to 1)
freq_raster <- sum_raster / length(scenario_files)

# crop it to the smaller extent made in previous section
freq_cropped <- terra::crop(freq_raster, 
                            ext(b["xmin"], 
                                b["xmax"], 
                                b["ymin"], 
                                b["ymax"]))

#  Plot the selection frequency map
p_freq <- ggplot() +
  geom_spatraster_rgb(data = tiles, maxcell = 500000) +
  geom_spatraster(data = freq_cropped, maxcell = 500000) + 
  scale_fill_viridis_c(
    option = "plasma", 
    name = "Selection\nFrequency",
    na.value = "transparent",
    limits = c(0.01, 1)
  ) +
  labs(
    title = "Selection Frequency Across Scenarios"
  ) +
  theme_minimal() + # Uses default dark text so it's clearly visible
  coord_sf(
    xlim = c(b["xmin"], b["xmax"]), 
    ylim = c(b["ymin"], b["ymax"]), 
    expand = FALSE
  )

print(p_freq)
ggsave("Selection_Frequency_Map.pdf", plot = p_freq, width = 8, height = 6)

#================================================================
#irreplaceability map
#================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(prioritizr)

#I want to plot scenario with 50% target, 75% threshold, no constraint
# to find all scenarios

all_scenarios <- bind_rows(lapply(list.files(pattern = "^scenario_.*\\.rds$"), function(f) {
  d <- readRDS(f)
  data.frame(
    file = f,
    target = if (!is.null(d$target)) d$target else NA,
    threshold = if (!is.null(d$threshold)) d$threshold else NA,
    penalty = if (!is.null(d$penalty)) d$penalty else NA
  )
}))

print(all_scenarios) #scenario 7 matches

irreplace7 <- readRDS("scenario_7.rds")

#need to rebuild the problem
allSpecies <- rast("LKDERasterStack.tif")

cost_layer <- allSpecies[[1]]
cost_layer[!is.na(cost_layer)] <- 1
names(cost_layer) <- "cost"

p_problem <- problem(cost_layer, allSpecies) %>%
  add_min_set_objective() %>%
  add_relative_targets(irreplace7$target) %>%
  add_binary_decisions()

#get the solution
p_solution <- terra::unwrap(irreplace7$solution_raster)

# Calculates rank importance
importance_scores <- eval_rank_importance(p_problem, 
                                          p_solution, 
                                          n = 10,
                                          force = TRUE)

# plotting the results

importance_scores[importance_scores == 0] <- NA

p_rank <- ggplot() +
  geom_spatraster_rgb(data = tiles)+
  geom_spatraster(data = importance_scores) + 
  scale_fill_viridis_c(
    option = "plasma", 
    name = "Importance Rank",
    na.value = "transparent",
    n.breaks = 10
  ) +
  labs(
    title = "Irreplaceability Rank Importance Map",
    subtitle = "Scenario 7: Target 0.5 | Threshold 0.75 | No Boundary Penalty"
  ) +
  theme_minimal() +
  coord_sf(
    xlim = c(b["xmin"], b["xmax"]), 
    ylim = c(b["ymin"], b["ymax"]), 
    expand = FALSE
  )

print(p_rank)

ggsave("Scenario_7_Rank_Importance_Map.pdf", plot = p_rank, width = 8, height = 6)

#=====================================================
# prioritizr with flight paths as planning units
#=====================================================

library(prioritizr)
library(sf)
library(terra)
library(highs)
library(ggplot2)
library(ggspatial)
library(prettymapr)

flight_lines.pj <- st_read("flight_lines.gpkg")

allSpecies <- rast("LKDERasterStack.tif")
flight_lines.pj$cost <- 1
cost_column = "cost"

p_problem <- problem(flight_lines.pj, 
                     allSpecies, 
                     cost_column) %>%
  add_min_set_objective() %>%
  add_relative_targets(0.17) %>% # 
  add_binary_decisions() %>% 
  add_highs_solver(gap = 0.1, 
                   verbose = TRUE)

# Solve the network optimization
p_solution <- solve(p_problem)

ggplot(p_solution) +
  annotation_map_tile(type = "cartolight", zoom = 5)+
  geom_sf(aes(color = as.factor(solution_1)), 
          linewidth = 0.8) +
  scale_color_manual(values = c("0" = alpha("grey80", 
                                            0.3), 
                                "1" = "darkblue"), 
                     name = "Selected") +
  theme_minimal() +
  labs(title = "Prioritizr Flight Path Solution")

mapview(p_solution, zcol = "solution_1", layer.name = "Prioritized Paths")
