# code for investigating spring migration phenology

library(dplyr)
library(lubridate)
library(ggplot2)
library(ggridges)

fullData <- read.csv ("StationPairsFiltered.csv")

springData <- fullData %>% 
  filter(
    season == "Spring Migration",
  )
nrow(springData)

springTable <- springData %>%  
  mutate(
    tsEnd_dt = as_datetime(tsEnd_dt, tz = "GMT"),
    # merging georgian bay with lake huron so they are not separate. 
    subbasin = if_else(subbasin == "geo_bay", "lk_huron", subbasin)
  )

# 1. Prepare and filter the data
erie_data <- springTable %>%
  # Remove year dependency to look at day-of-year
  mutate(month_day = format(tsEnd_dt, "%m-%d")) %>%
  filter(subbasin == "lk_erie") %>%
  group_by(species) %>%
  ungroup()

erie_plot_data <- erie_data %>%
  # Filter for your >= 3 threshold requirement
  group_by(species) %>%
  ungroup() %>%
  # Summarise to get the count of detections for every unique month-day
  group_by(species, month_day) %>%
  summarise(count = n(), .groups = "drop") %>%
  # Ensure the date is a proper format for the X-axis
  mutate(date = as.Date(paste0("2026-", month_day)))

# 2. Plotting
erieplot <- ggplot(erie_plot_data, aes(x = date, y = species, height = count, fill = species)) +
  geom_ridgeline(scale = 0.7, alpha = 0.6, color = "white") +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 weeks") +
  theme_ridges() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7)
  ) +
  labs(title = "Spring Phenology: Lake Erie", x = "Date", y = NULL)

library(plotly)
ggplotly(erieplot)

#================================================================
# map of spring first arrival points
#=================================================================

library(sf)
library(sfheaders)
library(mapview)

# springTable doesn't have coords so need to add that in 
first_arrival_paths <- springData %>%
  mutate(
    tsStart_dt = as_datetime(tsStart_dt, 
                             tz = "GMT"),
    subbasin = if_else(subbasin == "geo_bay", 
                       "lk_huron", subbasin)
  ) %>%
  semi_join(springTable, 
            by = c("species", 
                   "subbasin", 
                   "year_start", 
                   "tsStart_dt" = "first_arrival")) 

# separate current and previous station coords
pts_previous <- st_as_sf(first_arrival_paths, 
                         coords = c("lon_previous", 
                                    "lat_previous"), 
                         crs = 4326)
pts_current  <- st_as_sf(first_arrival_paths, 
                         coords = c("lon", 
                                    "lat"), 
                         crs = 4326)

# make the paths, setting to TRUE goes by row not distance
paths_first_arrival.sf <- pts_previous %>% 
  mutate(geometry = st_nearest_points(st_geometry(pts_previous), 
                                      st_geometry(pts_current), 
                                      pairwise = TRUE))
#add label for mapping
pts_prev_labeled <- pts_previous %>% mutate(location_type = "Previous")
pts_curr_labeled <- pts_current  %>% mutate(location_type = "Current")

combined_stations.sf <- bind_rows(pts_prev_labeled, 
                                  pts_curr_labeled) %>% 
  mutate(location_type = factor(location_type, 
                                levels = c("Previous", 
                                           "Current"))) %>% 
  st_as_sf()

subbasin_list <- c("lk_erie", 
                   "lk_ont", 
                   "lk_mich", "
                   lk_huron", 
                   "lk_sup")

final_map <- mapview()

# need to loop through the subbasins so they can be toggled in map
for (target_subbasin in subbasin_list) {
  
  # filter lines and stations
  subbasin_lines <- paths_first_arrival.sf %>% 
    filter(subbasin == target_subbasin)
  subbasin_dots  <- combined_stations.sf %>%
    filter(subbasin == target_subbasin)
  
  if (nrow(subbasin_dots) == 0) next
  
  # make unique names so they don't overwrite each other
  lines_name <- paste0(target_subbasin, 
                       " (Lines)")
  dots_name  <- paste0(target_subbasin,
                       " (Stations)")
  
  # Add lines and dots to the final map layout (hidden by default)
  final_map <- final_map + 
    mapview(
      subbasin_lines, 
      color = "black", 
      lwd = 1.5, 
      layer.name = lines_name,
      hide = TRUE
    ) + 
    mapview(
      subbasin_dots,
      zcol = "location_type", 
      #brown for previous, teal for current
      col.regions = c("Start" = "#846C5B", "End" = "#6DD3CE"), 
      size = 4,
      legend = FALSE,         
      layer.name = dots_name,
      hide = TRUE
    )
}

# Render the interactive map
final_map
