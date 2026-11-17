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
library(dplyr)
library(lubridate)

get_seasonal_station_effort_by_year <- function(df.recvDeps, season = "Spring") {
  
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
        
        # Define exact seasonal operating windows
        if (season %in% c("Spring", "Spring_Migration")) {
          win_start <- ymd_hms(paste0(yr, "-03-01 00:00:00"))
          win_end   <- ymd_hms(paste0(yr, "-06-16 23:59:59"))
        } else if (season %in% c("Fall", "Fall_Migration")) {
          win_start <- ymd_hms(paste0(yr, "-08-15 00:00:00"))
          win_end   <- ymd_hms(paste0(yr, "-11-15 23:59:59"))
        } else {
          stop("Specify 'Spring' or 'Fall'")
        }
        
        # Calculate overlap between deployment and the seasonal window for THIS year
        overlap_start <- max(start_dt, win_start)
        overlap_end   <- min(end_dt, win_end)
        
        if (overlap_start < overlap_end) {
          days_active <- as.numeric(difftime(overlap_end, 
                                             overlap_start, 
                                             units = "days"))
          
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
    # GROUP BY STATION AND YEAR (Preserves year-specific seasonal effort!)
    group_by(stationID, name, latitude, longitude, year) %>%
    summarise(active_days = sum(active_days_season), .groups = "drop")
}


# Generate Spring effort table (by station and year)
spring_effort_by_year <- get_seasonal_station_effort_by_year(DepsDuration, 
                                                             season = "Spring")

# Generate Fall effort table (by station and year)
fall_effort_by_year <- get_seasonal_station_effort_by_year(DepsDuration, 
                                                           season = "Fall")

# Save
saveRDS(spring_effort_by_year, "SpringDepsDuration.rds")
saveRDS(fall_effort_by_year, "FallDepsDuration.rds")

#add active days to the detections
# Join using the current station and the arrival year
spring_effort_lookup <- readRDS("SpringDepsDuration.rds") %>% 
  rename(active_days_spring = active_days) %>%
  select(stationID, year, active_days_spring)

fall_effort_lookup <- readRDS("FallDepsDuration.rds") %>% 
  rename(active_days_fall = active_days) %>%
  select(stationID, year, active_days_fall)

# Join t based on the record's season
flight_lines_weighted <- flight_lines.pj %>%
  mutate(
    arrival_dt = as.POSIXct(tsEnd_dt),
    year = lubridate::year(arrival_dt)
  ) %>%
  # Join Spring effort
  left_join(spring_effort_lookup, by = c("stationID", "year")) %>%
  # Join Fall effort
  left_join(fall_effort_lookup, by = c("stationID", "year")) %>%
  mutate(
    # Select the correct active days pool depending on the row's season
    active_days = case_when(
      grepl("Spring", season, ignore.case = TRUE) ~ active_days_spring,
      grepl("Fall", season, ignore.case = TRUE) ~ active_days_fall,
      TRUE ~ NA_real_
    ),
    # Fallback safety and final weight calculation
    active_days = tidyr::replace_na(active_days, 1),
    station_weight = 1 / log(active_days + 1)
  ) %>%
  # Clean up intermediate columns
  select(-active_days_spring, -active_days_fall)



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
