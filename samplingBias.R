# creating the sampling bias file to correct the flight data by station
# spatial sampling effort

library(dplyr)
library(lubridate)
library(stringr)

#station metadata files
bias_ON <- read.csv("./StationDownloads/stationMetadata/Motus-data-region-CA-ON_stations_downloaded-2026-11-04.csv")
bias_NY <- read.csv("./StationDownloads/stationMetadata/Motus-data-region-US-NY_stations_downloaded-2026-11-04.csv")
bias_PA <- read.csv("./StationDownloads/stationMetadata/Motus-data-region-US-PA_stations_downloaded-2026-11-04.csv")
bias_OH <- read.csv("./StationDownloads/stationMetadata/Motus-data-region-US-OH_stations_downloaded-2026-11-04.csv")
bias_IN <- read.csv("./StationDownloads/stationMetadata/Motus-data-region-US-IN_stations_downloaded-2026-11-04.csv")
bias_IL <- read.csv("./StationDownloads/stationMetadata/Motus-data-region-US-IL_stations_downloaded-2026-11-04.csv")
bias_WI <- read.csv("./StationDownloads/stationMetadata/Motus-data-region-US-WI_stations_downloaded-2026-11-04.csv")
bias_MI <- read.csv("./StationDownloads/stationMetadata/Motus-data-region-US-MI_stations_downloaded-2026-11-04.csv")
bias_MN <- read.csv("./StationDownloads/stationMetadata/Motus-data-region-US-MN_stations_downloaded-2026-11-04.csv")

bias_raw <- bind_rows(bias_ON,
                      bias_NY,
                      bias_PA,
                      bias_OH,
                      bias_IN,
                      bias_IL,
                      bias_WI,
                      bias_MI,
                      bias_MN) %>%  
  distinct() 

# Format the date columns using lubridate functions

bias_cleaned <- bias_raw %>%
  mutate(
    # 1. Clean and parse start date
    dtStart_raw = na_if(dtStart, "null"),
    dtStart_clean = str_remove(dtStart_raw, " \\(.*?\\)$"),
    dtStart_dt = parse_date_time(dtStart_clean, orders = "a b d Y H:M:S z", tz = "UTC"),
    
    # 2. Clean and parse lastData (endpoint)
    dtEnd_raw = case_when(
      lastData %in% c("null", "false", "") ~ NA_character_,
      TRUE ~ lastData
    ),
    dtEnd_clean = str_remove(dtEnd_raw, " \\(.*?\\)$"),
    dtEnd_dt = parse_date_time(dtEnd_clean, orders = "a b d Y H:M:S z", tz = "UTC"),
    
    # \Calculate duration:
    # Option 1: Strict duration\
    duration_days = as.numeric(difftime(dtEnd_dt, dtStart_dt, units = "days")),
    
    # Option 2: Assume active "null" stations run up to today's date
    duration_days_up_to_today = as.numeric(difftime(coalesce(dtEnd_dt, Sys.time()), dtStart_dt, units = "days"))
  )


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


saveRDS(DepsDuration,"DepsDuration.rds")
DepsDuration <- readRDS("DepsDuration.rds")

#summarize the effort by stationID rather than deployID
station_effort <- DepsDuration %>%
  group_by(stationID) %>%
  summarise(
    name = first(name),
    latitude = mean(latitude, na.rm = TRUE),
    longitude = mean(longitude, na.rm = TRUE),
    total_active_days = sum(active_days, na.rm = TRUE),
    n_deployments = n(),
    .groups = "drop"
  )

#==================================================
# creating spatial bias file with kernel denisty
#=================================================

library(sf)
library(dplyr)
library(terra)

station_effort_sf <- station_effort %>%
  filter(!is.na(longitude) & !is.na(latitude)) %>%
  st_as_sf(
    coords = c("longitude", "latitude"), 
    crs = 4326
  ) %>%
  st_transform(crs = 3978) # reproject to Statistics Canada Lambert

#use the bird LDKE to make sure its the same resolution and extent
template_raster <- rast("LKDEBirdRasterStack.tif") 

# convert station data to terra vector object
stations_vect <- vect(station_effort_sf)

# Extract coordinates and weights
coords <- crds(stations_vect)
weights <- stations_vect$total_active_days

# Generate the weighted Kernel Density Estimation
# Using terra's spatial smoothing (or kde functions depending on package versions)

bias_raster <- rasterize(stations_vect, 
                         template_raster, 
                         field = "total_active_days", 
                         fun = sum, 
                         background = 0)

# Apply a Gaussian smoothing
effortWeightMatrix <- focalMat(
  bias_raster, 
  d = 30000, 
  type = "Gauss"
) 

biasSurface <- focal(
  bias_raster, 
  w = effortWeightMatrix, 
  fun = sum, 
  na.rm = TRUE
)
# Normalize the bias surface 
bias_normalized <- biasSurface / global(biasSurface, 
                                        "mean", 
                                        na.rm = TRUE)[[1]]

writeRaster(bias_normalized, "biasRaster.tif",
            overwrite = TRUE)
