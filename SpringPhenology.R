# code for investigating spring migration phenology

library (dplyr)
library (lubridate)

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
  mutate(tsStart_dt = as_datetime(tsStart_dt, tz = "GMT")) %>%
    group_by(species, subbasin, year_start) %>% 
  summarise(first_arrival = min (tsStart_dt),
            .groups = "drop"
  )

