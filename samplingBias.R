# creating the sampling bias file to correct the flight data by station
# spatial sampling effort


#================================================================
# get the station deployments through motus
#====================================================================

#the duration of the station deployments will be used to create a bias file

library(motus)
library(dplyr)
library(lubridate)

Sys.setenv(TZ='GMT')
proj.num <- 437
srvTimeout(300)
sql.motus <- tagme(projRecv = proj.num, new = TRUE, 
                   update = TRUE) 
sql.motus <- DBI::dbConnect(RSQLite::SQLite(), "project-437.motus")

tbl.recvDeps <- tbl(sql.motus, "recvDeps")
# View information associated with a deployID
tbl_vars(tbl.recvDeps)

df.recvDeps <- tbl.recvDeps %>%
  collect() %>%
  as.data.frame() %>%
  mutate(timeStart = as_datetime(tsStart),
         timeEnd = as_datetime(tsEnd)) %>%
  dplyr::select(deployID, 
                receiverType,
                stationID,
                deviceID, 
                name, 
                latitude, 
                longitude, 
                isMobile, 
                timeStart, 
                timeEnd, 
                projectID, 
                elevation) 

df.recvDeps$yearEnd = year(df.recvDeps$timeEnd)
df.recvDeps %>% count(yearEnd)

DepsDuration <- df.recvDeps %>%
  mutate(
    start_date = as.POSIXct(timeStart),
    end_date = as.POSIXct(timeEnd),
    active_days = as.numeric(difftime(end_date,
                                      start_date, 
                                      units = "days"))
  )

### function to separate effort by migration season
get_seasonal_station_effort <- function(df.recvDeps, season = "Spring") {
  
  df.recvDeps %>%
    # Filter out rows with missing timestamps to prevent NA errors in year()
    filter(!is.na(timeStart) & !is.na(timeEnd)) %>%
    rowwise() %>%
    do({
      row <- .
      start_dt <- as.POSIXct(row$timeStart)
      end_dt   <- as.POSIXct(row$timeEnd)
      
      start_yr <- year(start_dt)
      end_yr   <- year(end_dt)
      
      season_records <- list()
      for (yr in start_yr:end_yr) {
        
        # Define the seasonal window for this year
        if (season %in% c("Spring", "Spring_Migration")) {
          win_start <- ymd_hms(paste0(yr, "-03-01 00:00:00"))
          win_end   <- ymd_hms(paste0(yr, "-06-16 23:59:59"))
        } else if (season %in% c("Fall", "Fall_Migration")) {
          win_start <- ymd_hms(paste0(yr, "-08-15 00:00:00"))
          win_end   <- ymd_hms(paste0(yr, "-11-15 23:59:59"))
        } else {
          stop("Specify 'Spring' or 'Fall'")
        }
        
        # Calculate overlap between deployment and the seasonal window
        overlap_start <- max(start_dt, win_start)
        overlap_end   <- min(end_dt, win_end)
        
        if (overlap_start < overlap_end) {
          days_active <- as.numeric(difftime(overlap_end, overlap_start, units = "days"))
          
          season_records[[length(season_records) + 1]] <- data.frame(
            stationID = row$stationID,
            name = row$name,
            latitude = row$latitude,
            longitude = row$longitude,
            year = yr,
            active_days_season = days_active
          )
        }
      }
      
      if (length(season_records) > 0) {
        bind_rows(season_records)
      } else {
        data.frame()
      }
    }) %>%
    ungroup() %>%
    # Aggregate total active days across all years into a single row per station
    group_by(stationID, name, latitude, longitude) %>%
    summarise(active_days = sum(active_days_season), .groups = "drop")
}


# Generate Spring effort 
spring_effort_stations <- get_seasonal_station_effort(DepsDuration, season = "Spring")

# Generate Fall effort 
fall_effort_stations <- get_seasonal_station_effort(DepsDuration, season = "Fall")

saveRDS(spring_effort_stations,"SpringDepsDuration.rds")
saveRDS(fall_effort_stations,"FallDepsDuration.rds")

#==================================================
# creating spatial bias file with kernel denisty
#=================================================

library(sf)
library(dplyr)
library(terra)

spring_effort_sf <- spring_effort_stations %>%
  filter(!is.na(longitude) & !is.na(latitude)) %>%
  st_as_sf(
    coords = c("longitude", "latitude"), 
    crs = 4326
  ) %>%
  st_transform(crs = 3978) # reproject to Statistics Canada Lambert

fall_effort_sf <- fall_effort_stations %>%
  filter(!is.na(longitude) & !is.na(latitude)) %>%
  st_as_sf(
    coords = c("longitude", "latitude"), 
    crs = 4326
  ) %>%
  st_transform(crs = 3978)

#use the bird LDKE to make sure its the same resolution and extent
template_raster <- rast("LKDEBirdRasterStack.tif") 

# convert station data to terra vector object
spring_stations_vect <- vect(spring_effort_sf)
fall_stations_vect <- vect(fall_effort_sf)

# Extract coordinates and weights
spring_coords <- crds(spring_stations_vect)
spring_weights <- spring_stations_vect$total_active_days

fall_coords <- crds(fall_stations_vect)
fall_weights <- fall_stations_vect$total_active_days

# Generate the weighted Kernel Density Estimation
# Using terra's spatial smoothing 

springbias_raster <- rasterize(spring_stations_vect, 
                         template_raster, 
                         field = "active_days", 
                         fun = sum, 
                         background = 0)

# Apply a Gaussian smoothing
springEffortWeightMatrix <- focalMat(
  springbias_raster, 
  d = 20000, 
  type = "Gauss"
) 

springbiasSurface <- focal(
  springbias_raster, 
  w = springEffortWeightMatrix, 
  fun = sum, 
  na.rm = TRUE
)

springbiasSurface[springbiasSurface < 0.1] <- 1.0

writeRaster(springbiasSurface, "springbiasRaster.tif",
            overwrite = TRUE)

# fall migration
fallbias_raster <- rasterize(fall_stations_vect, 
                         template_raster, 
                         field = "active_days", 
                         fun = sum, 
                         background = 0)

# Apply a Gaussian smoothing
fallEffortWeightMatrix <- focalMat(
  fallbias_raster, 
  d = 20000, 
  type = "Gauss"
) 

fallbiasSurface <- focal(
  fallbias_raster, 
  w = fallEffortWeightMatrix, 
  fun = sum, 
  na.rm = TRUE
)


writeRaster(fallbiasSurface, "fallbiasRaster.tif",
            overwrite = TRUE)

#######################################################3

#SEE IF THIS WORKS INSTEAD!

library(motus)
library(dplyr)
library(lubridate)

Sys.setenv(TZ='GMT')
proj.num <- 437
srvTimeout(300)
sql.motus <- tagme(projRecv = proj.num, new = TRUE, update = TRUE) 
sql.motus <- DBI::dbConnect(RSQLite::SQLite(), "project-437.motus")

tbl.recvDeps <- tbl(sql.motus, "recvDeps")

df.recvDeps <- tbl.recvDeps %>%
  collect() %>%
  as.data.frame() %>%
  mutate(timeStart = as_datetime(tsStart),
         timeEnd = as_datetime(tsEnd)) %>%
  dplyr::select(deployID, 
                receiverType,
                stationID,
                deviceID, 
                name, 
                latitude, 
                longitude, 
                isMobile, 
                timeStart, 
                timeEnd, 
                projectID, 
                elevation) 

DepsDuration <- df.recvDeps %>%
  mutate(
    start_date = as.POSIXct(timeStart),
    end_date = as.POSIXct(timeEnd),
    active_days = as.numeric(difftime(end_date,
                                      start_date, 
                                      units = "days"))
  )

### function to separate effort by migration season
get_seasonal_station_effort <- function(df.recvDeps, season = "Spring") {
  
  df.recvDeps %>%
    filter(!is.na(timeStart) & !is.na(timeEnd)) %>%
    rowwise() %>%
    do({
      row <- .
      start_dt <- as.POSIXct(row$timeStart)
      end_dt   <- as.POSIXct(row$timeEnd)
      
      start_yr <- year(start_dt)
      end_yr   <- year(end_dt)
      
      season_records <- list()
      for (yr in start_yr:end_yr) {
        
        if (season %in% c("Spring", "Spring_Migration")) {
          win_start <- ymd_hms(paste0(yr, "-03-01 00:00:00"))
          win_end   <- ymd_hms(paste0(yr, "-06-16 23:59:59"))
        } else if (season %in% c("Fall", "Fall_Migration")) {
          win_start <- ymd_hms(paste0(yr, "-08-15 00:00:00"))
          win_end   <- ymd_hms(paste0(yr, "-11-15 23:59:59"))
        } else {
          stop("Specify 'Spring' or 'Fall'")
        }
        
        overlap_start <- max(start_dt, win_start)
        overlap_end   <- min(end_dt, win_end)
        
        if (overlap_start < overlap_end) {
          days_active <- as.numeric(difftime(overlap_end, overlap_start, units = "days"))
          
          season_records[[length(season_records) + 1]] <- data.frame(
            stationID = row$stationID,
            name = row$name,
            latitude = row$latitude,
            longitude = row$longitude,
            year = yr,
            active_days_season = days_active
          )
        }
      }
      
      if (length(season_records) > 0) {
        bind_rows(season_records)
      } else {
        data.frame()
      }
    }) %>%
    ungroup() %>%
    group_by(stationID, name, latitude, longitude) %>%
    summarise(active_days = sum(active_days_season), .groups = "drop")
}

# Generate and save Spring effort 
spring_effort_stations <- get_seasonal_station_effort(DepsDuration, season = "Spring")
saveRDS(spring_effort_stations, "SpringDepsDuration.rds")

# Generate and save Fall effort 
fall_effort_stations <- get_seasonal_station_effort(DepsDuration, season = "Fall")
saveRDS(fall_effort_stations, "FallDepsDuration.rds")
