# comparing the solutions for spring and fall

#================================================
# within group comparisons
#================================================
library(terra)     
library(sf)        
library(dplyr)     
library(ggplot2)   
library(tidyterra) 
library(maptiles)

compare_selection_frequency <- function(group_name) {
  
  spring_group <- group_name
  fall_group <- group_name
  
  if (group_name == "nocturnal_takeoff") {
    fall_group <- "takeoff_birds"
  } else if (group_name == "takeoff_birds") {
    spring_group <- "nocturnal_takeoff"
  } else if (group_name == "nocturnal_land") {
    fall_group <- "land_birds"
  } else if (group_name == "land_birds") {
    spring_group <- "nocturnal_land"
  } else if (group_name == "day_birds") {
    fall_group <- "diurnal_birds"
  } else if (group_name == "diurnal_birds") {
    spring_group <- "day_birds"
  } else if (group_name == "all_bats") {
    fall_group <- "bats"
  } else if (group_name == "bats") {
    spring_group <- "all_bats"
  }
  spring_dir <- file.path("prioritizrOutput", 
                          "Spring", 
                          spring_group)
  fall_dir   <- file.path("prioritizrOutput", 
                          "Fall Migration", 
                          fall_group)
  
  if (!dir.exists(spring_dir) || !dir.exists(fall_dir)) {
    stop(sprintf("Directory missing. Checked:\n- Spring: %s\n- Fall: %s", 
                 spring_dir, 
                 fall_dir))
  }
  
  # Load Spring files and calc frequency
  spring_files <- list.files(path = spring_dir, 
                             pattern = "^scenario_.*\\.rds$", 
                             full.names = TRUE)
  if (length(spring_files) == 0) stop(sprintf("No Spring scenarios found for %s", 
                                              spring_group))
  
  first_spring <- readRDS(spring_files[1])
  spring_sum <- terra::unwrap(first_spring$solution_raster)
  if (length(spring_files) > 1) {
    for (i in 2:length(spring_files)) {
      spring_sum <- spring_sum + terra::unwrap(readRDS(spring_files[i])$solution_raster)
    }
  }
  spring_freq <- spring_sum / length(spring_files)
  
  # Load Fall files and calc frequency
  fall_files <- list.files(path = fall_dir, 
                           pattern = "^scenario_.*\\.rds$", 
                           full.names = TRUE)
  if (length(fall_files) == 0) stop(sprintf("No Fall scenarios found for %s", 
                                            fall_group))
  
  first_fall <- readRDS(fall_files[1])
  fall_sum <- terra::unwrap(first_fall$solution_raster)
  if (length(fall_files) > 1) {
    for (i in 2:length(fall_files)) {
      fall_sum <- fall_sum + terra::unwrap(readRDS(fall_files[i])$solution_raster)
    }
  }
  fall_freq <- fall_sum / length(fall_files)
  
  # Define extent bounding box (-94.46 to -74.24 Lon, 40 to 48 Lat in EPSG 3978)
  bbox_wgs84 <- sf::st_bbox(c(xmin = -94.46,
                              xmax = -74.24, 
                              ymin = 40, 
                              ymax = 48), 
                            crs = 4326)
  bbox_proj  <- sf::st_transform(sf::st_as_sfc(bbox_wgs84), 
                                 crs = 3978)
  
  # Fetch the Esri background tiles
  tiles <- maptiles::get_tiles(
    x = bbox_proj, 
    provider = "Esri.WorldGrayCanvas", 
    zoom = 8, 
    crop = TRUE
  )
  
  # Get bounding box directly from the tiles
  b <- sf::st_bbox(tiles)
  
  # Crop and align rasters to the tile bounding box
  spring_cropped <- terra::crop(spring_freq, 
                                terra::ext(b["xmin"], 
                                           b["xmax"], 
                                           b["ymin"], 
                                           b["ymax"]))
  fall_cropped <- terra::crop(fall_freq, terra::ext(b["xmin"], 
                                                    b["xmax"], 
                                                    b["ymin"], 
                                                    b["ymax"]))
  fall_aligned <- terra::resample(fall_cropped, 
                                  spring_cropped, 
                                  method = "near")
  
#classification of pixes
  cat_raster <- spring_cropped * 0
  
  s_valid <- !is.na(spring_cropped) & (spring_cropped > 0)
  f_valid <- !is.na(fall_aligned) & (fall_aligned > 0)
  
  is_both  <- s_valid & f_valid
  is_spring <- s_valid & !f_valid
  is_fall <- !s_valid & f_valid
  
  cat_raster[is_spring] <- 1
  cat_raster[is_fall] <- 2
  cat_raster[is_both] <- 3
  cat_raster[!(is_spring | is_fall | is_both)] <- NA
  
  # Convert to factor 
  cat_raster <- as.factor(cat_raster)
  levels(cat_raster) <- data.frame(
    id = c(1, 2, 3), 
    category = c("Spring Only", "Fall Only", "Both Seasons")
  )
  names(cat_raster) <- "Seasonality"
  
  # Plot with na.translate = FALSE to make background fully transparent
  p <- ggplot() +
    tidyterra::geom_spatraster_rgb(data = tiles) +
    tidyterra::geom_spatraster(data = cat_raster, na.rm = TRUE, alpha = 0.9) + 
    scale_fill_manual(
      values = c(
        "Spring Only" = "#2E5266", 
        "Fall Only" = "#E9C46A", 
        "Both Seasons" = "#519872"
      ),
      name = "Migration Seasonality",
      na.value = "transparent",
      na.translate = FALSE,
      drop = FALSE
    ) +
    theme_minimal() +
    coord_sf(
      crs = 3978,
      xlim = c(b["xmin"], b["xmax"]), 
      ylim = c(b["ymin"], b["ymax"]), 
      expand = FALSE
    ) 
  
  # Save output
  output_dir <- file.path("prioritizrOutput", 
                          "Seasonal_Comparisons")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pdf_filename <- file.path(output_dir, paste0("Selection_Cat_Comparison_", 
                                               group_name, 
                                               ".pdf"))
  ggsave(pdf_filename, plot = p, width = 8, height = 6)
  
  message(sprintf("Categorical selection comparison map successfully saved to '%s'", 
                  pdf_filename))
  return(p)
}

# All nocturnal birds
compare_selection_frequency(group_name = "nocturnal_birds")

# Nocturnal birds takeoff
compare_selection_frequency(group_name = "takeoff_birds")

# Nocturnal birds landing
compare_selection_frequency(group_name = "land_birds")

# Diurnal birds
compare_selection_frequency(group_name = "day_birds")

# Bats
compare_selection_frequency(group_name = "all_bats")

#===========================================================
# identifying seasonal hotspots
#===========================================================
library(terra)    
library(sf)       
library(maptiles)
library(tidyterra)  
library(ggplot2)

synthesize_migration_hotspots <- function() {
  
  # Define group lists and their matching directory 
  guild_map <- list(
    list(spring = "nocturnal_birds", 
         fall = "nocturnal_birds"),
    list(spring = "nocturnal_takeoff", 
         fall = "takeoff_birds"),
    list(spring = "nocturnal_land",    
         fall = "land_birds"),
    list(spring = "day_birds",         
         fall = "diurnal_birds"),
    list(spring = "all_bats",          
         fall = "bats")
  )
  
  # mean frequency across scenario for a group/season
  get_group_freq <- function(dir_path) {
    files <- list.files(path = dir_path, 
                        pattern = "^scenario_.*\\.rds$", 
                        full.names = TRUE)
    if (length(files) == 0) stop(sprintf("No scenarios found in %s", 
                                         dir_path))
    
    sum_rast <- terra::unwrap(readRDS(files[1])$solution_raster)
    if (length(files) > 1) {
      for (i in 2:length(files)) {
        sum_rast <- sum_rast + terra::unwrap(readRDS(files[i])$solution_raster)
      }
    }
    return(sum_rast / length(files))
  }
  
  # Define extent bounding box (-94.46 to -74.24 Lon, 40 to 48 Lat in EPSG 3978)
  bbox_wgs84 <- sf::st_bbox(c(xmin = -94.46, 
                              xmax = -74.24, 
                              ymin = 40, ymax = 48), 
                            crs = 4326)
  bbox_proj  <- sf::st_transform(sf::st_as_sfc(bbox_wgs84), 
                                 crs = 3978)
  
  # Esri background tiles
  tiles <- maptiles::get_tiles(x = bbox_proj, 
                               provider = "Esri.WorldGrayCanvas", 
                               zoom = 8, 
                               crop = TRUE)
  b <- sf::st_bbox(tiles)
  ext_template <- terra::ext(b["xmin"], 
                             b["xmax"], 
                             b["ymin"], 
                             b["ymax"])
  
  # Template raster
  sample_raw <- get_group_freq(file.path("prioritizrOutput", 
                                         "Spring", 
                                         "nocturnal_birds"))
  template_rast <- terra::rast(ext_template, 
                               resolution = terra::res(sample_raw), 
                               crs = "EPSG:3978")
  
  spring_rasters <- list()
  fall_rasters <- list()
  
  # Loop through, crop, and resample to template grid
  for (g in guild_map) {
    s_dir <- file.path("prioritizrOutput", 
                       "Spring", 
                       g$spring)
    f_dir <- file.path("prioritizrOutput", 
                       "Fall Migration", 
                       g$fall)
    
    # Load and crop Spring
    s_raw <- get_group_freq(s_dir)
    s_raw[s_raw < 0.5] <- 0
    s_raw[s_raw >= 0.5] <- 1
    s_crop <- terra::crop(s_raw, ext_template)
    s_resamp <- terra::resample(s_crop, 
                                template_rast, 
                                method = "near")
    
    # Load and crop Fall
    f_raw <- get_group_freq(f_dir)
    f_raw[f_raw < 0.5] <- 0
    f_raw[f_raw >= 0.5] <- 1
    f_crop <- terra::crop(f_raw, 
                          ext_template)
    f_resamp <- terra::resample(f_crop, 
                                template_rast, 
                                method = "near")
    
    spring_rasters[[g$spring]] <- s_resamp
    fall_rasters[[g$fall]] <- f_resamp
  }
  
  # Stack and sum across all groups within each season
  spring_stack <- terra::rast(spring_rasters)
  fall_stack   <- terra::rast(fall_rasters)
  
  spring_hotspot <- terra::app(spring_stack, 
                               fun = "sum", 
                               na.rm = TRUE)
  fall_hotspot   <- terra::app(fall_stack, 
                               fun = "sum", 
                               na.rm = TRUE)
  
  # Fully combined master stack
  master_stack <- c(spring_stack, 
                    fall_stack)
  master_hotspot <- terra::app(master_stack, 
                               fun = "sum", 
                               na.rm = TRUE)
  
  # Mask out zero-values so they become NAs
  spring_hotspot[spring_hotspot <= 0] <- NA
  fall_hotspot[fall_hotspot <= 0] <- NA
  master_hotspot[master_hotspot <= 0] <- NA
  
  spring_hotspot <- terra::as.factor(spring_hotspot)
  fall_hotspot <- terra::as.factor(fall_hotspot)
  master_hotspot <- terra::as.factor(master_hotspot)
  
  # Save rasters 
  out_dir <- file.path("prioritizrOutput", 
                       "Synthesis")
  if (!dir.exists(out_dir)) dir.create(out_dir, 
                                       recursive = TRUE)
  
  terra::writeRaster(spring_hotspot, 
                     file.path(out_dir, 
                               "Spring_MultiGuild_Hotspots.tif"), 
                     overwrite = TRUE)
  terra::writeRaster(fall_hotspot, 
                     file.path(out_dir, 
                               "Fall_MultiGuild_Hotspots.tif"), 
                     overwrite = TRUE)
  terra::writeRaster(master_hotspot, 
                     file.path(out_dir, 
                               "Master_Invariant_Core_Refugia.tif"), 
                     overwrite = TRUE)
  
#plotting spring, fall, and both combined
  plot_hotspot <- function(rast_obj, 
                           title_text, 
                           filename_suffix) {
    
    p <- ggplot() +
      tidyterra::geom_spatraster_rgb(data = tiles) +
      tidyterra::geom_spatraster(data = rast_obj, 
                                 na.rm = TRUE, 
                                 alpha = 0.95) + 
      scale_fill_viridis_d(
        option = "inferno",
        direction = -1,
        name = "Overlapping\nGuilds",
        na.value = "transparent",
        na.translate = FALSE
      ) +
      theme_minimal() +
      coord_sf(
        crs = 3978,
        xlim = c(b["xmin"], 
                 b["xmax"]), 
        ylim = c(b["ymin"], 
                 b["ymax"]), 
        expand = FALSE
      ) +
      labs(title = title_text)
    
    pdf_path <- file.path(out_dir, 
                          paste0("Synthesis_", 
                                 filename_suffix, 
                                 ".pdf"))
    ggsave(pdf_path, 
           plot = p, 
           width = 8, 
           height = 6)
    return(p)
  }
  
  p_spring <- plot_hotspot(spring_hotspot, 
                           "Spring", 
                           "Spring_Hotspots")
  p_fall   <- plot_hotspot(fall_hotspot, 
                           "Fall", 
                           "Fall_Hotspots")
  p_master <- plot_hotspot(master_hotspot, 
                           "Both Seasons", 
                           "Master_Refugia")
  
  message("Discrete synthesis and plotting complete! Files saved to 'prioritizrOutput/Synthesis/'")
  
  return(list(
    spring_raster = spring_hotspot, 
    fall_raster = fall_hotspot, 
    master_raster = master_hotspot,
    spring_plot = p_spring,
    fall_plot = p_fall,
    master_plot = p_master
  ))
}

hotspot_results <- synthesize_migration_hotspots()

