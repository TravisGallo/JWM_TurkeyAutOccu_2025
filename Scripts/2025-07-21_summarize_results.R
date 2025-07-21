source("./Scripts/TurkeyMultiOccu_package_load.R")

turkey_dat <- read.csv("./Data/2020-2023_turkey_obs.csv")

# number of observations
nrow(turkey_dat)

# number of sites
n_distinct(turkey_dat$locationAbbr)

# what sites?
unique(turkey_dat$fullName)

# average number of detections per season
mean_season <- turkey_dat |>
  group_by()

# load in detection history to calculate frequency of detection by season

# load in the detection histories
det <- read.csv("./Data/TurkeyOccupancyReport.csv")

sites_data <- det[1:77,3:5]

# remove weird last column
det2 <- det[ ,-ncol(det)]

# get total days detected by site and season
det2$total_det <- rowSums(det2[,6:ncol(det2)], na.rm = TRUE)

sum(det2$total_det) / 14

# group by season and get the mean
result <- det2 %>%
  group_by(Season) %>%
  summarise(across(Day_1:Day_45, ~ sum(.x, na.rm = TRUE), .names = "sum_{col}")) |>
  mutate(across(-Season, ~ ifelse(. > 1, 1, .))) |>
  mutate(row_sum = rowSums(across(-Season), na.rm = TRUE))

sum(result$row_sum, na.rm = TRUE) / 14

result[,1:ncol(result)]

seasons <- c("SU","FA","WI","SP","SU","FA","WI","SP","SU","FA","WI","SP","SU","FA")

result$season <- seasons

season_avg <- result |>
  group_by(season) |>
  summarize(mean = mean(row_sum))

# odds ratio

# roads
exp(0.2608580)

# canopy height
exp(-0.5416810)

# distance to water
exp(-0.4282965)

# habitat
exp(0.4013105)






# Table 1
covs <- read.csv("./Data/2024-10-17_covariates.csv")

ranges <- apply(covs, 2, range)

# for map

sites_data <- det[1:77,3:5]

all_sites <- sites_data$Site

obs_by_site <- turkey_dat |>
  group_by(locationAbbr) |>
  summarize(n = sum(numIndividuals), .groups = "drop") |>
  complete(locationAbbr = all_sites, fill = list(n = 0)) |>
  arrange(locationAbbr) |>
  mutate(Latitude = sites_data$Latitude,
         Longitude = sites_data$Longitude)

# Convert to an sf object
sf_data <- st_as_sf(obs_by_site, coords = c("Longitude", "Latitude"), crs = 4326)

st_write(sf_data, dsn = "./Results",
         layer = "turkey_observations",
         driver = "ESRI Shapefile")

esri_light_gray <- "https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}"


tm <- tm_shape(sf_data) +
  tm_symbols(
    size = "n",               # Map point size to the 'value' column
    breaks = c(0, 10, 50, 100, 200, 500),    # Set size range (e.g., min size = 1, max size = 5)
    col = "blue",                 # Set point color
    alpha = 0.7                   # Adjust transparency
  ) +
  tm_tiles(esri_light_gray) +     # Add ESRI Light Gray basemap
  tm_layout(title = "Map with ESRI Light Gray Basemap")

# Display the map
tmap_mode("view") # Use interactive mode
tm
