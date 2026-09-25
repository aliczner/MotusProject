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
