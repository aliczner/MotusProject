library(motus)
library(dplyr)
library(lubridate)

Sys.setenv(TZ = "GMT")
proj.num <- 437
srvTimeout(300)

# Connect to SQLite database and pull global receiver metadata
sql.motus <- tagme(projRecv = proj.num, new = FALSE, update = FALSE)
metadata(sql.motus)

df.recvDeps <- tbl(sql.motus, "recvDeps") %>%
  collect() %>%
  as.data.frame() %>%
  mutate(
    timeStart = as_datetime(tsStart),
    # Replace NA end dates (currently active towers) with data download cutoff date
    timeEnd   = if_else(is.na(tsEnd), ymd_hms("2026-09-09 23:59:59"), as_datetime(tsEnd))
  ) %>%
  dplyr::select(
    deployID, receiverType, stationID, deviceID, name,
    latitude, longitude, isMobile, timeStart, timeEnd, projectID, elevation
  )

#====================================================================
# Calculate year-specific seasonal station effort
#====================================================================

get_seasonal_station_effort_by_year <- function(df.recvDeps, season = "Spring") {
  
  max_window_days <- if (season %in% c("Spring", "Spring_Migration")) 108 else 93
  
  df.recvDeps %>%
    filter(!is.na(timeStart) & !is.na(timeEnd) & !is.na(stationID)) %>%
    rowwise() %>%
    do({
      row <- .
      start_dt <- as.POSIXct(row$timeStart, tz = "GMT")
      end_dt   <- as.POSIXct(row$timeEnd, tz = "GMT")
      
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
      
      if (length(season_records) > 0) bind_rows(season_records) else data.frame()
    }) %>%
    ungroup() %>%
    # Group by stationID and year (omitting name/lat/lon so minor coordinate tweaks don't split a station)
    group_by(stationID, year) %>%
    summarise(
      name        = first(na.omit(name)),
      latitude    = first(na.omit(latitude)),
      longitude   = first(na.omit(longitude)),
      active_days = pmin(sum(active_days_season), max_window_days),
      .groups     = "drop"
    )
}

# Generate and save Spring and Fall effort tables (by station and year)
spring_effort_by_year <- get_seasonal_station_effort_by_year(df.recvDeps, season = "Spring")
fall_effort_by_year   <- get_seasonal_station_effort_by_year(df.recvDeps, season = "Fall")

saveRDS(spring_effort_by_year, "SpringDepsDuration.rds")
saveRDS(fall_effort_by_year,   "FallDepsDuration.rds")

DBI::dbDisconnect(sql.motus)

#====================================================================
# Verification check: 
#====================================================================

flight_lines_no_lp %>%
  sf::st_drop_geometry() %>%
  mutate(year = lubridate::year(as.POSIXct(tsEnd_dt))) %>%
  left_join(spring_effort_by_year %>% select(stationID, year, active_days),
            by = c("stationID", "year")) %>%
  summarise(total_spring = n(), unmatched_spring = sum(is.na(active_days)))

fall_flight_lines_no_lp %>%
  sf::st_drop_geometry() %>%
  mutate(year = lubridate::year(as.POSIXct(tsEnd_dt))) %>%
  left_join(fall_effort_by_year %>% select(stationID, year, active_days),
            by = c("stationID", "year")) %>%
  summarise(total_fall = n(), unmatched_fall = sum(is.na(active_days)))
