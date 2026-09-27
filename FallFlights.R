# code for investigating fall migration phenology

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

fallData <- fullData %>% 
  filter(
    season == "Fall Migration",
  )
nrow(fallData) #10399 obs

fallTable <- fallData %>%  
  mutate(
    tsEnd_dt = as_datetime(tsEnd_dt, tz = "GMT"),
    # merging georgian bay with lake huron so they are not separate. 
    subbasin = if_else(subbasin == "geo_bay", "lk_huron", subbasin)
  )

#=========================================================
# creating flight paths
#=========================================================

library(sfheaders)

flight_steps <- fallTable %>%
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

#===========================================================
# creating a figure of search effort
# ======================================================
library(dplyr)
library(tidyr)
library(ggplot2)

flight_lines %>%
  select(lon, lon_previous, lon_tagSite) %>%
  pivot_longer(
    cols = c(lon, lon_previous, lon_tagSite),
    names_to = "longitude_type",
    values_to = "longitude"
  ) %>%
  mutate(
    longitude_type = case_match(
      longitude_type,
      "lon" ~ "Current Station",
      "lon_previous" ~ "Previous Station",
      "lon_tagSite" ~ "Tag Site"
    )
  ) %>%
  filter(!is.na(longitude) & longitude < 0) %>%
  ggplot(aes(x = longitude, 
             color = longitude_type, 
             fill = longitude_type)) +
  geom_density(alpha = 0.3, 
               linewidth = 0.9) +
  facet_wrap(~ longitude_type, 
             ncol = 1, 
             scales = "free_y") +
  scale_color_manual(values = c(
    "Current Station" = "#2b5c8f", 
    "Previous Station" = "#2d6a4f", 
    "Tag Site" = "#d95f02"
  )) +
  scale_fill_manual(values = c(
    "Current Station" = "#2b5c8f", 
    "Previous Station" = "#2d6a4f", 
    "Tag Site" = "#d95f02"
  )) +
  labs(
    title = "Longitude Distribution - Fall",
    x = "Longitude (°W)",
    y = "Density of Flights"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    strip.text = element_text(size = 11)
  )

#================================================
# checking flight timing
#=================================================

flight_time_summary <- flight_lines %>%
  group_by(MigrateTime, species, diel_period) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(MigrateTime, species) %>%
  mutate(
    total_flights = sum(n),
    percentage = (n / total_flights) * 100,
    species_label = paste0(species, " (", total_flights, ")")
  ) %>%
  ungroup()

library(ggplot2)

# diurnal barplot
plot_diurnal <- flight_time_summary %>%
  filter(MigrateTime == "diurnal") %>%
  ggplot(aes(x = reorder(species_label, percentage), 
             y = percentage, 
             fill = diel_period)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_manual (
    values = c (
      "daylight" = "#ECA72C",
      "night" = "#31263E"
    )) +
  labs(
    title = "Fall Diurnal",
    x = "",
    y = "Percentage of Flights (%)",
    fill = "Diel Period"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# Mixed  Plot
plot_mixed <- flight_time_summary %>%
  filter(MigrateTime == "mixed") %>%
  ggplot(aes(x = reorder(species_label, percentage), 
             y = percentage, 
             fill = diel_period)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_manual (
    values = c (
      "daylight" = "#ECA72C",
      "night" = "#31263E"
    )) +
  labs(
    title = "Fall Mixed",
    x = "",
    y = "Percentage of Flights (%)",
    fill = "Diel Period"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# nocturnal
plot_nocturnal <- flight_time_summary %>%
  filter(MigrateTime == "nocturnal") %>%
  ggplot(aes(x = reorder(species_label, percentage), 
             y = percentage, 
             fill = diel_period)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_manual (
    values = c (
      "daylight" = "#ECA72C",
      "night" = "#31263E"
    )) +
  labs(
    title = "Fall Nocturnal",
    x = "",
    y = "Percentage of Flights (%)",
    fill = "Diel Period"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

#==============================================================
# partitioning night flights for nocturnal migrants
#==============================================================
# If the flight starts in the morning  its true starting sunset 
# was actually on the previous calendar day (- 1 day).
# Otherwise, it belongs to the sunset of the current calendar day.

### nocturnal take off
noc.takeoff <- flight_lines %>%  
  filter(MigrateTime == "nocturnal",
         Animal == "Bird") %>%  
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

### nocturnal landing

noc.landing <- flight_lines %>%  
  filter(
    MigrateTime == "nocturnal",
    Animal == "Bird") %>%  
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


