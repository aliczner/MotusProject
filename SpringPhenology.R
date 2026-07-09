# code for investigating spring migration phenology

library(dplyr)
library(lubridate)
library(ggplot2)

fullData <- read.csv ("StationPairsFiltered.csv")

#========================================================================
# filtering data based on tagging site
#=======================================================================

# will remove based on tagging location any flags 2, 3 and 4
  # flag 4 = missing tagging location info
  # flag 3 = tagged during summer or winter
  # flag 2 = tagged during migration but > 5 km away from shoreline


springData <- fullData %>% 
  filter(
    season == "Spring Migration",
    tagSite_Flags %in% c("None", "Flag_1"),
    current_in_GLWS == "TRUE"
  )
nrow(springData)

springTable <- springData %>%  
  mutate(
    tsStart_dt = as_datetime(tsStart_dt, tz = "GMT"),
    # merging georgian bay with lake huron so they are not separate. 
    subbasin = if_else(subbasin == "geo_bay", "lk_huron", subbasin)
  ) %>%
  group_by(species, subbasin, year_start) %>%  
  summarise(
    first_arrival = min(tsStart_dt),
    .groups = "drop"
  )

#============================================================
# horizontal bar chart of spring arrival timing by subbasin
#============================================================
library(ggrepel)

# rename the subbasins
subbasin_names <- c(
  "lk_ont"   = "Lake Ontario",
  "lk_erie"  = "Lake Erie",
  "lk_huron" = "Lake Huron",
  "lk_sup"   = "Lake Superior",
  "lk_mich"  = "Lake Michigan"
)

# get each subbasin to loop through to make a plot per subbasin
subbasin_list <- unique(springTable$subbasin[!is.na(springTable$subbasin)])

for (target_subbasin in subbasin_list) {
  
  clean_title <- if_else(
    target_subbasin %in% names(subbasin_names), 
    subbasin_names[target_subbasin], 
    target_subbasin
  )
  
  # normalizing dates by julian day
  spring_normalised <- springTable %>%
    filter(subbasin == target_subbasin) %>%
    filter(!is.na(first_arrival)) %>% 
    mutate(
      dummy_date = as.Date(format(first_arrival, "2026-%m-%d")),
      julian_day = yday(first_arrival),
      arrival_year = format(first_arrival, "%Y")
    )
  
  if (nrow(spring_normalised) == 0) next
  
  # summarizing the data to make the date range bars
  subbasin_data <- spring_normalised %>%
    group_by(species) %>%
    summarise(
      earliest = min(dummy_date),
      latest   = max(dummy_date),
      record_count = n(),
      .groups = "drop"
    ) %>%
    mutate(
      latest_visual = if_else(record_count == 1, 
                              earliest + days(1),#to ensure it is shown
                              latest)
    ) %>%
    arrange(earliest, species) %>%
    mutate(
      row_num = row_number(),
      #to colour each row differently
      stripe_group = if_else(row_num %% 2 == 0, "Even", "Odd"), 
      species = factor(species, 
                       levels = unique(species))
    )
  
  # put together the new variables for plotting
  spring_normalised <- spring_normalised %>%
    left_join(select(subbasin_data, species, stripe_group), by = "species") %>%
    mutate(species = factor(species, levels = levels(subbasin_data$species)))
  
  # now plotting
  p <- ggplot() + 
    geom_segment(
      data = subbasin_data,
      aes(x = earliest, 
          xend = latest_visual, 
          y = species, 
          yend = species, 
          colour = stripe_group), 
      linewidth = 4, 
      alpha = 0.25,
      show.legend = FALSE
    ) +
    geom_point(
      data = spring_normalised,
      aes(x = dummy_date, 
          y = species, 
          colour = stripe_group),
      size = 2.5,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    geom_text_repel(
      data = spring_normalised,
      aes(x = dummy_date, 
          y = species, 
          label = arrival_year),
      size = 2.2,                  
      colour = "black",             
      segment.size = 0.2,          
      segment.color = "grey50",
      direction = "both",          
      max.overlaps = 50            
    ) +
    scale_colour_manual(values = c("Odd" = "#631D76", 
                                  "Even" = "#62A87C")) +
    scale_x_date(
      date_labels = "%b %d", 
      date_breaks = "1 week", #which dates are shown
      date_minor_breaks = "1 day"
    ) +
    labs(
      title = paste("Spring Arrival Phenology -", clean_title), 
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_line(colour = "grey92", 
                                        linewidth = 0.3), 
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.x = element_line(colour = "grey90", 
                                        linewidth = 0.4), 
      axis.line.x = element_line(colour = "grey50", 
                                 linewidth = 0.5),      
      axis.line.y = element_line(colour = "grey50", 
                                 linewidth = 0.5),
      axis.ticks.x = element_line(colour = "grey50", 
                                  linewidth = 0.5),
      axis.ticks.length.x = unit(0.15, 
                                 "cm"),
      axis.text.y = element_text(size = 8, 
                                 face = "italic", 
                                 colour = "black"), 
      axis.text.x = element_text(angle = 45, 
                                 hjust = 1, 
                                 colour = "black")
    )
  #export pdfs
  file_name <- paste0("springArrival_", target_subbasin, ".pdf")
  plot_height <- max(6, nrow(subbasin_data) * 0.25)
  
  ggsave(
    filename = file_name, 
    plot = p, 
    width = 11, 
    height = plot_height,
    device = "pdf"
  )

}

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
