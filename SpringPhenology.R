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


