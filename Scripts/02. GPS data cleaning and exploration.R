## ########################################################################## ##
## 02. GPS data cleaning                                                      ##
##                                                                            ##
## Project: MoveSIN GPS data cleaning and exploration                         ##
##                                                                            ##
## Created by: Luke Emerson                                                   ##
## Created: 25th April 2026                                                   ##
##                                                                            ##
## Edited by: Luke Emerson                                                    ##
## Edited: 23rd June 2026                                                     ##
## ########################################################################## ##

# Clear environment ------------------------------------------------------------
rm(list = ls())
gc()

# Load packages ----------------------------------------------------------------
library(dplyr)
library(sp)
library(move2)
library(ggplot2)
library(sf)
library(lubridate)
library(geosphere)
library(adehabitatLT)
library(circular)
library(car)
library(stringr)
library(leaflet)
library(htmltools)
library(ctmm)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(cowplot)
library(glue)
library(httr)
library(utils)
library(terra)
library(tmap)
library(gstat)
library(suncalc)
library(readr)

## ########################################################################## ##
# MULTIPLE FERAL CATS ----------------------------------------------------------
## ########################################################################## ##

cat("\n================ Multiple feral cats ================\n")

# Check working directory ------------------------------------------------------
getwd()

# Load simulated feral cat data ------------------------------------------------
# NOTE THIS IS NOT REAL DATA AND SHOULD NOT BE USED FOR FORMAL ANALYSES!!!
cat("\nNOTE: THIS IS NOT REAL DATA AND SHOULD NOT BE USED FOR FORMAL ANALYSES!!!\n")

# Load the list
files <- list.files("Telemetry_data", full.names = TRUE)
files

# Choose the file from the list
cat_data <- read.csv(files[1]) 

# View data structure ----------------------------------------------------------
str(cat_data)
cat_data %>%
  dplyr::distinct(Device.Name, Device.ID)
plot(cat_data$Longitude, cat_data$Latitude)
nrow(cat_data)

# Initial cleaning -------------------------------------------------------------

# Renaming columns -------------------------------------------------------------
# Makes valid R column names but already done by read.csv in this case, so no change
colnames(cat_data)
names(cat_data) <- make.names(names(cat_data), unique = TRUE) 
colnames(cat_data)

# Make lower case column names
cat_data <- cat_data %>%
  rename_with(tolower)

# Reformat column names
colnames(cat_data)
colnames(cat_data) <- str_replace_all(colnames(cat_data), "\\.+", ".")  # Collapse multiple dots
colnames(cat_data) <- str_replace_all(colnames(cat_data), "\\.$", "")   # Trim trailing dot
colnames(cat_data)

# Count number of rows
nrow(cat_data)

# Keep only complete cases to remove NA lat and long
cat_data <- cat_data %>%
  filter(complete.cases(latitude, longitude))

# Count number of rows
nrow(cat_data)

# Look at the headers and the data structures
head(cat_data)
summary(cat_data)

cat("\nDOP (Dilution of Precision) is a unitless ratio that describes how satellite
geometry affects the precision of a GPS location estimate.\n")


# Define timezone --------------------------------------------------------------
str(cat_data)
cat_data <- cat_data %>%
  # Convert the GMT column to POSIXct with correct timezone
  mutate(date.time.gmt = as.POSIXct(date.time.gmt, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")) %>%
  # Create a local time column in Perth time (GMT+8)
  mutate(date.time.local = with_tz(date.time.gmt, tzone = "Australia/Perth"))
str(cat_data)

# Remove any duplicates --------------------------------------------------------
# Order data by individual and time
cat_data <- cat_data[order(cat_data$device.id, cat_data$date.time.local),]

# Check for duplicates
which(cat_data$device.id[-nrow(cat_data)] == cat_data$device.id[-1] &
        cat_data$date.time.local[-nrow(cat_data)] == cat_data$date.time.local[-1])

# Remove duplicates if required - keeps the first occurrence
cat_data <- cat_data %>%
  distinct(device.id, date.time.local, .keep_all = TRUE)

# Instead you may wish to keep the location with lowest DOP
cat_data <- cat_data %>%
  arrange(device.id, date.time.local, dop) %>%
  distinct(device.id, date.time.local, .keep_all = TRUE)


# Remove data for individuals based on when collars were deployed --------------
head(cat_data)
cats <- unique(cat_data$device.name)
cats

# Vector of individual IDs
individuals <- c("C25F14", "C25M32", "C25M33", "C25M34", "C25M35", "C25M36", "C25M37", "C25M38")

# Corresponding deployment dates 
deployment_start <- as.POSIXct(c(
  "1984-05-28 18:00:00",
  "1984-05-25 18:00:00",
  "1984-05-27 18:00:00",
  "1984-05-27 18:00:00",
  "1984-05-28 18:00:00",
  "1984-05-29 18:00:00",
  "1984-06-01 18:00:00",
  "1984-06-02 18:00:00"
))

# Corresponding retrieval (end) dates
deployment_end <- as.POSIXct(c(
  "1984-08-31 18:00:00",
  "1984-09-15 18:00:00",
  "1984-09-15 18:00:00",
  "1984-09-15 18:00:00",
  "1984-09-15 18:00:00",
  "1984-08-31 18:00:00",
  "1984-09-15 18:00:00",
  "1984-09-15 18:00:00"
))

# Remove pre-deployment and post-deployment fixes
# Remove 3 days at start of deployment to reduce potential effects of trapping on movement

# Create lookup table
deployment_df <- data.frame(
  device.name = individuals,
  start_date = deployment_start,
  end_date = deployment_end
)


# Diagnostics for deployment date filtering ------------------------------------
# Keep original data before filtering
cat_data_raw <- cat_data

# Check deployment dates and observed fix ranges
print(
  cat_data_raw %>%
    left_join(deployment_df, by = "device.name") %>%
    group_by(device.name) %>%
    summarise(
      n = n(),
      first_fix = min(date.time.local, na.rm = TRUE),
      last_fix = max(date.time.local, na.rm = TRUE),
      deploy_start = first(start_date),
      deploy_start_plus_3days = first(start_date) + lubridate::days(3),
      deploy_end = first(end_date),
      .groups = "drop"
    ),
  width = Inf
)

# Count retained observations per individual
cat_data_raw %>%
  left_join(deployment_df, by = "device.name") %>%
  mutate(
    retained = date.time.local >= start_date + lubridate::days(3) &
      date.time.local <= end_date
  ) %>%
  group_by(device.name) %>%
  summarise(
    total_fixes = n(),
    retained_fixes = sum(retained, na.rm = TRUE),
    removed_fixes = sum(!retained, na.rm = TRUE),
    .groups = "drop"
  )


# Apply final filtering ------------------------------------------------------
cat_data <- cat_data_raw %>%
  left_join(deployment_df, by = "device.name") %>%
  filter(
    date.time.local >= start_date + lubridate::days(3),
    date.time.local <= end_date
  )

nrow(cat_data)

cats2 <- unique(cat_data$device.name)
cats2
setdiff(cats, cats2)


# Check DOP and remove large potentially inaccurate values ---------------------
hist(cat_data$dop, breaks = 30, main = "Distribution of DOP", xlab = "DOP Value")
abline(v = 5, col = "red", lty = 2)
summary(cat_data$dop)

nrow(cat_data)
cat_data <- cat_data %>%
  filter(dop <= 10) # Arbitrary filter >10 typically considered unreliable
# Some methods/packages like ctmm account for error so no need to filter
nrow(cat_data)


# Visualise the data again -----------------------------------------------------
# Plot the data
plot(cat_data$longitude,cat_data$latitude,xlab="Longitude",ylab="Latitude")


# Interactive map -------------------------------------------------------------
# Force standard WGS84 lon/lat
cat.sf <- st_as_sf(
  cat_data,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

st_crs(cat.sf)

# Make device.name a factor
cat.sf$device.name <- as.factor(cat.sf$device.name)

# Create a color palette
pal <- colorFactor(palette = "Set1", domain = cat.sf$device.name)

# Plot in leaflet with hollow circles
leaflet(cat.sf) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(
    radius = 3,
    color = ~pal(device.name),
    fill = FALSE,
    weight = 2,
    label = ~paste(
      "Device:", device.name,
      "\nDate:", date.time.local
    ),
    popup = ~paste(
      "Device:", device.name,
      "<br>Date:", date.time.local
    )
  ) %>%
  addLegend(
    "topright",
    pal = pal,
    values = ~device.name,
    title = "Device"
  )

# Remove fixes after collar recovery
cat_data <- cat_data %>%
  filter(
    !(device.name %in% c("C25M36", "C25F14")) |
      date.time.local <= as.POSIXct("1984-08-23 23:59:59")
  )

# Check
cat_data %>%
  group_by(device.name) %>%
  summarise(
    start = min(date.time.local),
    end = max(date.time.local),
    n = n()
  )

# Force standard WGS84 long/lat
cat.sf <- st_as_sf(
  cat_data,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

st_crs(cat.sf)

# Make device.name a factor
cat.sf$device.name <- as.factor(cat.sf$device.name)

# Create a color palette
pal <- colorFactor(palette = "Set1", domain = cat.sf$device.name)

# Plot in leaflet with hollow circles
leaflet(cat.sf) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(
    radius = 3,
    color = ~pal(device.name),
    fill = FALSE,
    weight = 2,
    label = ~paste(
      "Device:", device.name,
      "\nDate:", date.time.local
    ),
    popup = ~paste(
      "Device:", device.name,
      "<br>Date:", date.time.local
    )
  ) %>%
  addLegend(
    "topright",
    pal = pal,
    values = ~device.name,
    title = "Device"
  )


# Investigate movement ---------------------------------------------------------
# Order data per individual and time
cat_data <- cat_data %>%
  arrange(device.name, date.time.local)

# Convert to dataframe
cat_df <- as.data.frame(cat_data)

# Ensure row names are unique
rownames(cat_df) <- NULL

# Create move2 object in WGS84 lon/lat
cat.move <- mt_as_move2(
  cat_df,
  coords = c("longitude", "latitude"),
  time = "date.time.local",
  track_id = "device.name"
)

# Assign WGS84 longitude/latitude CRS
st_crs(cat.move) <- 4326

# Check
st_crs(cat.move)

# Transform to projected CRS (metres) for distances/speeds
cat.move.utm <- st_transform(
  cat.move,
  7850
)

# Check projected CRS
st_crs(cat.move.utm)

# Plot geographic tracks (degrees)
plot(
  cat.move$geometry,
  xlab = "Longitude",
  ylab = "Latitude",
  main = "GPS tracks"
)

# Plot projected tracks (metres)
plot(
  st_geometry(cat.move.utm),
  xlab = "Easting (m)",
  ylab = "Northing (m)",
  main = "Projected tracks"
)

# Split data into a track per individual
tracks <- split(cat.move.utm, cat.move.utm$device.name)
names(tracks)
class(tracks[[1]])

# Calculate step lengths and time lags
get_valid_steps <- function(track) {
  
  df <- as.data.frame(track)
  
  coords <- st_coordinates(track)
  
  df <- df %>%
    mutate(
      x = coords[,1],
      y = coords[,2]
    ) %>%
    arrange(date.time.local) %>%
    mutate(
      previous_x = lag(x),
      previous_y = lag(y),
      previous_time = lag(date.time.local),
      
      time_diff = as.numeric(
        difftime(
          date.time.local,
          previous_time,
          units = "hours"
        )
      ),
      
      step_dist = sqrt(
        (x - previous_x)^2 +
          (y - previous_y)^2
      ),
      
      step_km = step_dist / 1000
    ) %>%
    filter(
      is.na(time_diff) | 
        (time_diff >= 0.25 & time_diff <= 1.25)
    )
  
  return(df)
}


# Apply to all individuals
valid_steps_all <- bind_rows(
  lapply(tracks, get_valid_steps)
)


# Step length distribution -  have not filtered out particular individuals or step lengths
ggplot(valid_steps_all, aes(x = step_km)) +
  geom_histogram(
    binwidth = 0.1,
    fill = "#009999",
    colour = "black"
  ) +
  labs(
    title = "Distribution of hourly step distances",
    x = "Step distance (km)",
    y = "Count"
  ) +
  theme_minimal()

# OR a density plot -  have not filtered out particular individuals or step lengths
ggplot(valid_steps_all, aes(x = step_km)) +
  geom_density(fill = "#009999") +
  labs(
    title = "Density of Hourly Step Distances",
    x = "Step distance (km)",
    y = "Density"
  ) +
  theme_minimal()

# Step summary
summary(valid_steps_all$step_km)


# Filter implausible steps based on biological knowledge -----------------------
# valid_steps_all <- valid_steps_all %>%
#   filter(step_km <= 3)

valid_steps_all <- valid_steps_all %>%
  mutate(step_hour = hour(date.time.local)) # Use the hour after the displacement has occurred rather than before - more accurate visual representation

head(valid_steps_all)
unique(valid_steps_all$device.name)

# Displacement for all individuals pooled using raw values
ggplot(valid_steps_all, aes(x = factor(step_hour), y = step_km)) +
  geom_boxplot(fill = "#009999", color = "black", outlier.alpha = 0.3) +
  labs(
    title = "Hourly Displacement",
    x = "Hour of Day",
    y = "Displacement (km)"
  ) +
  theme_minimal()


# Calculate mean displacement per individual per hour
mean_disp <- valid_steps_all %>%
  group_by(device.name, step_hour) %>%
  summarize(mean_step_km = mean(step_km, na.rm = TRUE), .groups = "drop")


# Plot line graph
ggplot(mean_disp, aes(x = step_hour, y = mean_step_km, color = device.name, group = device.name)) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    title = "Mean Hourly Displacement per Individual",
    x = "Hour of Day",
    y = "Mean Displacement (km)",
    color = "Individual"
  ) +
  scale_x_continuous(breaks = 0:23) +
  theme_minimal()

# Question
cat("\nQUESTION: Any red flags???\n")


# Assess displacement per week per individual to determine individual variability ----
valid_steps_all <- valid_steps_all %>%
  mutate(
    week = isoweek(date.time.local),       # Week of year
    year = year(date.time.local),          # Add year so weeks don't reset
    week_id = paste(year, week, sep = "-"), 
    day = as.Date(date.time.local)         # Exact date
  )

# Summarise by week individual and hour displacement
mean_disp_week <- valid_steps_all %>%
  group_by(device.name, week_id, step_hour) %>%
  summarise(
    mean_step_km = mean(step_km, na.rm = TRUE),
    .groups = "drop"
  )

# Plot hourly displacement per week per individual
ggplot(mean_disp_week, aes(x = step_hour, y = mean_step_km, color = device.name, group = device.name)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~ week_id, ncol = 3) +
  labs(
    title = "Mean Hourly Displacement per Week",
    x = "Hour of Day",
    y = "Mean Step Distance (km)",
    color = "Individual"
  ) +
  scale_x_continuous(breaks = 0:23) +
  theme_minimal()  +
  theme(
    axis.text.x = element_text(size = 6),   # change tick label size
    axis.title.x = element_text(size = 12), # change x-axis title size
    axis.text.y = element_text(size = 8),   # optional: y-axis tick size
    axis.title.y = element_text(size = 12)  # optional: y-axis title size
  )

# QUESYION
cat("\nQUESTION:Some individuals show unusually small displacement... possible reasons?\n")


# Isolate and investigate C25F14 C25M33 and C25M36 data
# Filter to just the individuals of interest
valid_steps_subset <- valid_steps_all %>%
  filter(device.name %in% c("C25F14", "C25M36","C25M33"))

# Summarise by week, individual, and hour displacement
mean_disp_week <- valid_steps_subset %>%
  group_by(device.name, week_id, step_hour) %>%
  summarize(mean_step_km = mean(step_km, na.rm = TRUE), .groups = "drop")

# Plot hourly displacement per week per individual
ggplot(mean_disp_week, aes(x = step_hour, y = mean_step_km, color = device.name, group = device.name)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~ week_id, ncol = 3) +
  labs(
    title = "Mean Hourly Displacement per Week",
    x = "Hour of Day",
    y = "Mean Step Distance (km)",
    color = "Individual"
  ) +
  scale_x_continuous(breaks = 0:23) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 12)
  )

# Based on this, drop all data for C25F14 and C25M36
# Keep only data up to and including week 30 for C25M33
# Could use more finite approach than this in terms of data filtering



# Create final cleaned GPS dataset ---------------------------------------------
# Remove individuals with unusable movement data
# - Remove all fixes for C25F14 and C25M36
# - Retain C25M33 only until end of July 1984

cat_data_clean <- cat_data %>%
  filter(
    !(device.name %in% c("C25F14", "C25M36"))
  ) %>%
  filter(
    !(device.name == "C25M33" &
        date.time.local > as.POSIXct("1984-07-31 23:59:59"))
  )

# Check final individuals and sampling periods
cat_data_clean %>%
  group_by(device.name) %>%
  summarise(
    n = n(),
    start = min(date.time.local),
    end = max(date.time.local),
    .groups = "drop"
  )

# Create spatial object (WGS84)
cat.sf <- st_as_sf(
  cat_data_clean,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

# Create movement object ------------------------------------------------------
cat_df <- as.data.frame(cat_data_clean)

cat.move <- mt_as_move2(
  cat_df,
  coords = c("longitude", "latitude"),
  time = "date.time.local",
  track_id = "device.name"
)

st_crs(cat.move) <- 4326

# Transform to projected CRS for distances -------------------------------------
cat.move.utm <- st_transform(
  cat.move,
  7850
)

# Visualise --------------------------------------------------------------------
# Use WGS84 coordinates for leaflet
cat.move.leaflet <- cat.move %>%
  arrange(device.name, date.time.local)


# Create track lines (LINESTRING) for each individual
tracks_leaflet <- cat.move.leaflet %>%
  group_by(device.name) %>%
  summarise(
    geometry = st_combine(geometry),
    .groups = "drop"
  ) %>%
  st_cast("LINESTRING")


# Colour palette by individual
pal <- colorFactor(
  palette = "Set1",
  domain = cat.move.leaflet$device.name
)


# Plot tracks + locations
leaflet() %>%
  addProviderTiles("Esri.WorldImagery") %>%
  
  # Track lines
  addPolylines(
    data = tracks_leaflet,
    color = ~pal(device.name),
    weight = 2,
    opacity = 0.8
  ) %>%
  
  # GPS locations
  addCircleMarkers(
    data = cat.move.leaflet,
    radius = 3,
    color = ~pal(device.name),
    fill = FALSE,
    weight = 2,
    opacity = 0.8,
    popup = ~paste(
      "Individual:", device.name,
      "<br>Date:", date.time.local
    )
  ) %>%
  
  addLegend(
    position = "topright",
    pal = pal,
    values = cat.move.leaflet$device.name,
    title = "Individual"
  )


# Split tracks by individual --------------------------------------------------
tracks <- split(
  cat.move.utm,
  cat.move.utm$device.name
)


# Calculate step lengths and time lags ----------------------------------------
get_valid_steps <- function(track) {
  
  df <- as.data.frame(track)
  
  coords <- st_coordinates(track)
  
  df <- df %>%
    mutate(
      x = coords[,1],
      y = coords[,2]
    ) %>%
    arrange(date.time.local) %>%
    mutate(
      previous_x = lag(x),
      previous_y = lag(y),
      previous_time = lag(date.time.local),
      
      time_diff = as.numeric(
        difftime(
          date.time.local,
          previous_time,
          units = "hours"
        )
      ),
      
      step_dist = sqrt(
        (x - previous_x)^2 +
          (y - previous_y)^2
      ),
      
      step_km = step_dist / 1000
    ) %>%
    filter(
      is.na(time_diff) |
        (time_diff >= 0.25 & time_diff <= 1.25)
    )
  
  return(df)
}


# Apply to all individuals
valid_steps_all <- bind_rows(
  lapply(tracks, get_valid_steps)
)

# Add movement timing variables
valid_steps_all <- valid_steps_all %>%
  mutate(
    step_hour = lubridate::hour(date.time.local),
    week = lubridate::isoweek(date.time.local),
    year = lubridate::year(date.time.local),
    week_id = paste(year, week, sep = "-"),
    day = as.Date(date.time.local)
  )


# Final check of displacement
valid_steps_all %>%
  group_by(device.name) %>%
  summarise(
    n_steps = n(),
    mean_step_km = mean(step_km, na.rm = TRUE),
    max_step_km = max(step_km, na.rm = TRUE),
    .groups = "drop"
  )

# Summarise by week individual and hour displacement
mean_disp_week <- valid_steps_all %>%
  group_by(device.name, week_id, step_hour) %>%
  summarise(
    mean_step_km = mean(step_km, na.rm = TRUE),
    .groups = "drop"
  )

# Plot hourly displacement per week per individual
ggplot(mean_disp_week, aes(x = step_hour, y = mean_step_km, color = device.name, group = device.name)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~ week_id, ncol = 3) +
  labs(
    title = "Mean Hourly Displacement per Week",
    x = "Hour of Day",
    y = "Mean Step Distance (km)",
    color = "Individual"
  ) +
  scale_x_continuous(breaks = 0:23) +
  theme_minimal()  +
  theme(
    axis.text.x = element_text(size = 6),   # change tick label size
    axis.title.x = element_text(size = 12), # change x-axis title size
    axis.text.y = element_text(size = 8),   # optional: y-axis tick size
    axis.title.y = element_text(size = 12)  # optional: y-axis title size
  )

# Count number of rows
nrow(cat_data_clean)
# Started with 5081
# 5081-3619 = 1462 rows removed through intitial cleaning





# Extra ------------------------------------------------------------------------
# Activity analysis ------------------------------------------------------------
# Mean of means - pooled across individuals
mean_disp_all <- mean_disp %>%
  group_by(step_hour) %>%
  summarize(mean_hourly_disp = mean(mean_step_km, na.rm = TRUE), .groups = "drop")

str(mean_disp_all)

# Plot
ggplot(mean_disp_all, aes(x = step_hour, y = mean_hourly_disp)) +
  geom_line(color = "darkred", linewidth = 1.2) +
  geom_point(color = "darkred", size = 2) +
  labs(
    title = "Mean Hourly Displacement (Equal Weighting per Individual)",
    x = "Hour of Day",
    y = "Mean Displacement (km)"
  ) +
  scale_x_continuous(breaks = 0:23) +
  theme_minimal()

# Smoothed
ggplot(mean_disp_all, aes(x = step_hour, y = mean_hourly_disp)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "darkred", size = 2) +
  geom_smooth(method = "loess", se = TRUE, span = 0.5, color = "black", fill = "grey70", linetype = "dashed") +
  labs(
    title = "Hourly Displacement: Individual Means and Smoothed Trend",
    x = "Hour of Day",
    y = "Mean Displacement (km)"
  ) +
  scale_x_continuous(breaks = 0:23) +
  theme_minimal()

detach("package:ctmm", unload = TRUE)

# With shading for night
ggplot2::ggplot(mean_disp_all, aes(x = step_hour, y = mean_hourly_disp)) +
  # Night shading using annotate
  annotate("rect", xmin = 17, xmax = 23, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.2) +
  annotate("rect", xmin = 0, xmax = 7, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.2) +
  # Data layers
  geom_line(color = "#EA2E08", linewidth = 1) +
  geom_point(color = "#EA2E08", size = 2) +
  geom_smooth(method = "loess", se = TRUE, span = 0.5,
              color = "black", fill = "#009999", linetype = "dashed") +
  labs(
    title = "Hourly Displacement: Individual Means and Smoothed Trend",
    x = "Hour of Day",
    y = "Mean Displacement (km)"
  ) +
  scale_x_continuous(breaks = 0:23, limits = c(0, 23)) +
  theme_minimal()

colnames(valid_steps_all)






## ########################################################################## ##
# SINGLE TIGER -----------------------------------------------------------------
## ########################################################################## ##

cat("\n================ Single tiger  ================\n")

# Load simulated tiger data ----------------------------------------------------
# NOTE THIS IS NOT REAL DATA AND SHOULD NOT BE USED FOR FORMAL ANALYSES!!!
cat("\nNOTE: THIS IS NOT REAL DATA AND SHOULD NOT BE USED FOR FORMAL ANALYSES!!!\n")

# Load the list
files <- list.files("Telemetry_data", full.names = TRUE)
files

# Choose the file from the list
tiger_data <- read.csv(files[2]) 

# View data structure ----------------------------------------------------------
head(tiger_data)

# Initial cleaning -------------------------------------------------------------
# Makes valid R column names but already done by read.csv in this case, so no change
colnames(tiger_data)
names(tiger_data) <- make.names(names(tiger_data), unique = TRUE) 
colnames(tiger_data)

# Make lower case column names
tiger_data <- tiger_data %>%
  rename_with(tolower)

# Reformat column names
colnames(tiger_data)
colnames(tiger_data) <- str_replace_all(colnames(tiger_data), "\\.+", ".")  # Collapse multiple dots
colnames(tiger_data) <- str_replace_all(colnames(tiger_data), "\\.$", "")   # Trim trailing dot
colnames(tiger_data)

# Count number of rows
nrow(tiger_data)

# Keep only complete cases to remove NA lat and long
tiger_data <- tiger_data %>%
  filter(complete.cases(latitude, longitude))

# Count number of rows
nrow(tiger_data)

# Look at the headers and the data structures
head(tiger_data)
summary(tiger_data)


# Format Timestamps ------------------------------------------------------------
tiger_data$date.time.gmt <- as.POSIXct(tiger_data$date.time.gmt, tz = "UTC")
tiger_data$date.time.local <- as.POSIXct(tiger_data$date.time.local, tz = "Asia/Thimphu")


# Remove duplicate timestamps --------------------------------------------------
# Order data by individual and time
tiger_data <- tiger_data %>% arrange(device.id, date.time.local)

# Check for duplicates
tiger_data %>%
  group_by(device.id, date.time.local) %>%
  filter(n() > 1)

tiger_data %>%
  group_by(device.id, date.time.local) %>%
  filter(n() > 1) %>%
  nrow()

# Remove if required
# Only remove if the ID and Time are identical - keeps first occurrence of any duplicate
tiger_data <- tiger_data %>%
  arrange(device.id, date.time.local) %>%
  distinct(device.id, date.time.local, .keep_all = TRUE)

# Check no duplicates remain
tiger_data %>%
  group_by(device.id, date.time.local) %>%
  filter(n() > 1) %>%
  nrow()


# Calculate Metrics ------------------------------------------------------------
tiger_data <- tiger_data %>%
  arrange(device.id, date.time.local) %>%
  group_by(device.id) %>%
  mutate(
    # Time interval in hours
    dt_hrs = as.numeric(
      difftime(date.time.local, lag(date.time.local), units = "hours")
    ),
    
    # Distance between consecutive fixes (metres)
    dist_m = distHaversine(
      cbind(longitude, latitude),
      cbind(lag(longitude), lag(latitude))
    )
  )

# Generate Summary Table -------------------------------------------------------
data.frame(
  First_Fix      = min(tiger_data$date.time.local, na.rm = TRUE),
  Last_Fix       = max(tiger_data$date.time.local, na.rm = TRUE),
  Total_Fixes    = nrow(tiger_data),
  Days_Tracked   = round(as.numeric(diff(range(tiger_data$date.time.local))), 2),
  Median_Int_Hr  = median(tiger_data$dt_hrs, na.rm = TRUE),
  Min_Gap_Hr     = min(tiger_data$dt_hrs, na.rm = TRUE),
  Max_Gap_Hr     = max(tiger_data$dt_hrs, na.rm = TRUE),
  Mean_Dist_m    = mean(tiger_data$dist_m, na.rm = TRUE),
  Median_Dist_m  = median(tiger_data$dist_m, na.rm = TRUE),
  Min_Dist_m     = min(tiger_data$dist_m, na.rm = TRUE),
  Max_Dist_m     = max(tiger_data$dist_m, na.rm = TRUE),
  Missing_Data   = sum(is.na(tiger_data$latitude))
)

# Calculate fixes per year -----------------------------------------------------
fixes_per_year <- tiger_data %>%
  group_by(year = lubridate::year(date.time.local)) %>%
  summarise(Count = n()) %>%
  ungroup()

fixes_per_year

# Visual check of sampling effort over time
ggplot(fixes_per_year, aes(x = factor(year), y = Count)) +
  geom_col(fill = "steelblue") +
  labs(title = "Fixes per Year", x = "Year", y = "Number of Fixes") +
  theme_minimal()

# Plotting ------------------------------------------------------------------
# Time Intervals over Time (Detects Gaps & Duty Cycle shifts)
ggplot(tiger_data, aes(x = date.time.local, y = dt_hrs)) +
  geom_line(color = "grey80") +
  geom_point(alpha = 0.5, size = 1, color = "firebrick") +
  labs(title = "Time Intervals Between Fixes",
       subtitle = "Spikes indicate data gaps or missing fixes",
       x = "Project Timeline", y = "Interval (Hours)") +
  theme_minimal()

# Spatial Trajectory (Quick check for outliers)
# Ensure the year is extracted and treated as a factor (discrete categories)
tiger_data$year <- factor(year(tiger_data$date.time.local))

# Spatial Trajectory colored by Year
ggplot(tiger_data, aes(x = longitude, y = latitude)) +
  geom_path(alpha = 0.3, color = "grey50") +
  geom_point(aes(color = year), alpha = 0.7, size = 1.5) +
  # Using a discrete color scale (Set1, Dark2, etc.)
  scale_color_brewer(palette = "Set1") + 
  coord_quickmap() +
  labs(title = "Spatial Trajectory Check",
       subtitle = "Colored by Year",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

# What does this suggest?
cat("\nQUESTION: What does this suggest?\n")



# Remove data by year ----------------------------------------------------------
tiger_data <- tiger_data %>%
  filter(year != "1936")

# OR filter directly using the timestamp (use if 'year' column has not been created)
tiger_data <- tiger_data %>%
  filter(year(date.time.local) != 1936)

# Check the result
table(year(tiger_data$date.time.local))

# Determine time interval between successive locations
tiger_data <- tiger_data %>%
  arrange(device.id, date.time.local) %>%
  group_by(device.id) %>%
  mutate(
    # Time interval in hours
    dt_hrs = as.numeric(
      difftime(date.time.local, lag(date.time.local), units = "hours")
    ),
    # Distance between consecutive fixes (metres)
    dist_m = distHaversine(
      cbind(longitude, latitude),
      cbind(lag(longitude), lag(latitude))
    )
  )

# Plot Time Intervals over the course of the study
ggplot(filter(tiger_data, !is.na(dt_hrs)), aes(x = date.time.local, y = dt_hrs)) +
  geom_line(color = "grey70", linewidth = 0.3 ) +
  # Fixes within 3 hours of each other considered normal i.e. median of 2 hrs * 1.5 = 3
  geom_point(aes(color = dt_hrs > (median(dt_hrs, na.rm=T) * 1.5)), alpha = 0.6) +
  scale_color_manual(values = c("steelblue", "firebrick"), 
                     name = "Gaps",
                     labels = c("Normal Interval", "Missed Fix/Gap")) +
  labs(title = "Sampling Interval Consistency",
       x = "Date", y = "Interval between fixes (Hours)") +
  theme_minimal()


# Histogram of Intervals (Identifies the primary sampling rate)
ggplot(tiger_data, aes(x = dt_hrs)) +
  geom_histogram(fill = "steelblue", color = "white", binwidth = 2) +
  labs(title = "Distribution of Sampling Intervals",
       x = "Interval (Hours)", y = "Frequency") +
  theme_minimal()

summary(tiger_data$dt_hrs)

# Filter to reasonable time gap
ggplot(
  tiger_data %>% filter(!is.na(dt_hrs)),
  aes(x = dt_hrs)
) +
  geom_histogram(
    fill = "steelblue",
    color = "white",
    binwidth = 1
  ) +
  coord_cartesian(xlim = c(0, 36)) +
  labs(
    title = "Distribution of Sampling Intervals (0–24 h)",
    x = "Interval (Hours)",
    y = "Frequency"
  ) +
  theme_minimal()



# Visualisation of locations ---------------------------------------------------
# Basic plot
plot(
  tiger_data$longitude,
  tiger_data$latitude,
  xlab = "Longitude",
  ylab = "Latitude"
)

# Interactive map 
# Convert to sf object (longitude/latitude decimal degrees)   
gps_sf_ll <- tiger_data %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,        # WGS84 GPS coordinates
    remove = FALSE
  )

# # Colour code by individual
# ids <- unique(gps_sf_ll$device.id)
# 
# pal <- colorFactor(
#   palette = topo.colors(length(ids)),
#   domain  = ids
# )

# Colour by date/time - converts date time to seconds since 1970
# Add the numeric time column to your spatial object
# Ensure chronological order and create fix sequence per individual
gps_sf_ll <- gps_sf_ll %>%
  arrange(device.id, date.time.local) %>%
  group_by(device.id) %>%
  mutate(
    order_num = row_number(),
    time_num = as.numeric(date.time.local)
  ) %>%
  ungroup()

# Option 1: Colour by individual
ids <- unique(gps_sf_ll$device.id)

pal_id <- colorFactor(
  palette = topo.colors(length(ids)),
  domain = ids
)

leaflet(gps_sf_ll) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addCircleMarkers(
    radius = 3,
    color = ~pal_id(device.id),
    stroke = FALSE,
    fillOpacity = 0.8,
    label = ~lapply(
      paste0(
        "<b>ID:</b> ", device.id,
        "<br><b>Fix:</b> ", order_num,
        "<br><b>Time:</b> ",
        format(date.time.local, "%Y-%m-%d %H:%M:%S")
      ),
      htmltools::HTML
    )
  ) %>%
  addLegend(
    "bottomright",
    pal = pal_id,
    values = ~device.id,
    title = "Individual"
  )


# Option 2: Colour by time progression
# Shows the absolute timing of fixes.
# Colours represent when locations occurred across the entire dataset, so useful
# for identifying temporal patterns, gaps, or changes over calendar time. 
# Points closer together in time will be similar colour
pal_time <- colorNumeric(
  palette = "viridis",
  domain = gps_sf_ll$time_num
)

leaflet(gps_sf_ll) %>%
  addTiles() %>%
  addCircleMarkers(
    radius = 3,
    color = ~pal_time(time_num),
    stroke = FALSE,
    fillOpacity = 0.6,
    label = ~lapply(
      paste0(
        "<b>ID:</b> ", device.id,
        "<br><b>Fix:</b> ", order_num,
        "<br><b>Time:</b> ",
        format(date.time.local, "%Y-%m-%d %H:%M:%S")
      ),
      htmltools::HTML
    )
  ) %>%
  addLegend(
    "bottomright",
    pal = pal_time,
    values = ~time_num,
    title = "Time progression"
  )


# Option 3: Colour by movement sequence within each individual
# Shows the relative order of fixes within each
# individual track. Colours represent fix sequence from first to last location,
# making it easier to visualise movement direction and progression regardless of
# the actual date/time.
# Points in sequence will be similar colour despite time difference
pal_order <- colorNumeric(
  palette = "plasma",
  domain = gps_sf_ll$order_num
)

leaflet(gps_sf_ll) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addCircleMarkers(
    radius = 3,
    color = ~pal_order(order_num),
    stroke = FALSE,
    fillOpacity = 0.6,
    label = ~lapply(
      paste0(
        "<b>ID:</b> ", device.id,
        "<br><b>Fix order:</b> ", order_num,
        "<br><b>Time:</b> ",
        format(date.time.local, "%Y-%m-%d %H:%M:%S")
      ),
      htmltools::HTML
    )
  ) %>%
  addLegend(
    "bottomright",
    pal = pal_order,
    values = ~order_num,
    title = "Progress (First → Last)"
  )



# Outlier detection ------------------------------------------------------------
# Review speed and step distributions
# Create utm coordinates
tiger_data_sf <- tiger_data %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(32646)   # GWGS 84 / UTM zone 46N

coords <- st_coordinates(tiger_data_sf)

tiger_data <- tiger_data_sf %>%
  mutate(
    x = coords[, "X"],
    y = coords[, "Y"]
  ) %>%
  st_drop_geometry()

glimpse(tiger_data)


# Calculate speeds between successive fixes
tiger_data <- tiger_data %>%
  arrange(device.id, date.time.local) %>%
  group_by(device.id) %>%
  mutate(
    dist_km  = sqrt((x - lag(x))^2 + (y - lag(y))^2) / 1000,   # convert m to km
    speed_kmh = dist_km / dt_hrs                            # km / hour
  ) %>%
  ungroup()

# Plot distribution of speeds per individual
ggplot(tiger_data, aes(x = speed_kmh)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  facet_wrap(~ device.id, scales = "free_y") +
  labs(x = "Speed (km/h)", y = "Count",
       title = "Distribution of Step Speeds per Individual")

ggplot(tiger_data, aes(x = speed_kmh)) +
  geom_density() +
  theme_bw() +
  facet_wrap(~ device.id, scales = "free_y") +
  labs(x = "Speed (km/h)", y = "Density",
       title = "Speed Density per Individual")

# Plot distribution of speeds combined across individuals
ggplot(tiger_data, aes(x = speed_kmh)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  labs(x = "Speed (km/h)", y = "Count",
       title = "Combined Speed Distribution")

ggplot(tiger_data, aes(x = speed_kmh)) +
  geom_density() +
  theme_bw() +
  labs(x = "Speed (km/h)", y = "Density",
       title = "Combined Speed Density")

# Review speed range
summary(tiger_data$speed_kmh)

# Define acceptable fix interval window around nominal 2-hour schedule
nominal_interval <- median(tiger_data$dt_hrs, na.rm = TRUE)

# Keep observations that are within 2.25 hrs of another
valid_steps <- tiger_data %>%
  filter(!is.na(dt_hrs)) %>%
  filter(dt_hrs <= nominal_interval * 1.125) # Up to 2.25 hrs

cat("Total fixes:", nrow(tiger_data), "\n")
cat("Valid steps:", nrow(valid_steps), "\n")
cat("Excluded steps:", nrow(tiger_data) - nrow(valid_steps), "\n")

# Calculate speed (km/h)
valid_steps <- valid_steps %>%
  mutate(
    speed_kmh = (dist_m / 1000) / dt_hrs
  )

# Summary statistics
cat("\nStep-length summary (m)\n")
print(summary(valid_steps$dist_m))

cat("\nSpeed summary (km/h)\n")
print(summary(valid_steps$speed_kmh))

# Distance thresholds
dist_95 <- quantile(valid_steps$dist_m, 0.95, na.rm = TRUE)
dist_99 <- quantile(valid_steps$dist_m, 0.99, na.rm = TRUE)

# Speed thresholds
speed_95 <- quantile(valid_steps$speed_kmh, 0.95, na.rm = TRUE)
speed_99 <- quantile(valid_steps$speed_kmh, 0.99, na.rm = TRUE)


# STEP LENGTH DISTRIBUTION
ggplot(valid_steps, aes(x = dist_m / 1000)) +
  geom_histogram(
    binwidth = 0.5,
    fill = "darkgreen",
    color = "white"
  ) +
  geom_vline(
    xintercept = dist_95 / 1000,
    colour = "blue",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = dist_99 / 1000,
    colour = "red",
    linewidth = 1
  ) +
  theme_bw() +
  labs(
    title = "Step Length Distribution",
    subtitle = "Only steps ≤ 2.25 h included",
    x = "Distance Between Consecutive Fixes (km)",
    y = "Count"
  )

# Flag long-distance steps
distance_outliers <- valid_steps %>%
  filter(dist_m > dist_99)

cat("95th percentile:", round(dist_95, 0), "m\n")
cat("99th percentile:", round(dist_99, 0), "m\n")
cat("Flagged steps:", nrow(distance_outliers), "\n")


# SPEED DISTRIBUTION
ggplot(valid_steps, aes(x = speed_kmh)) +
  geom_histogram(
    binwidth = 0.25,
    fill = "steelblue",
    color = "white"
  ) +
  geom_vline(
    xintercept = speed_95,
    colour = "blue",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = speed_99,
    colour = "red",
    linewidth = 1
  ) +
  theme_bw() +
  labs(
    title = "Speed Distribution",
    subtitle = "Only steps ≤ 2.25 h included",
    x = "Speed (km/h)",
    y = "Count"
  )

# Flag high-speed steps
speed_outliers <- valid_steps %>%
  filter(speed_kmh > speed_99)

cat("95th percentile:", round(speed_95, 2), "km/h\n")
cat("99th percentile:", round(speed_99, 2), "km/h\n")
cat("Flagged steps:", nrow(speed_outliers), "\n")


# Review top 10 steps and speeds
cat("\nTop 10 longest steps\n")
print(
  valid_steps %>%
    arrange(desc(dist_m)) %>%
    select(device.id, date.time.local, dt_hrs, dist_m, speed_kmh) %>%
    head(10)
)

cat("\nTop 10 fastest steps\n")
print(
  valid_steps %>%
    arrange(desc(speed_kmh)) %>%
    select(device.id, date.time.local, dt_hrs, dist_m, speed_kmh) %>%
    head(10)
)


# Review fix regularity --------------------------------------------------------
# To determine whether successive fixes per ID are within certain amount of time
# of each other
tiger_data <- tiger_data %>%
  arrange(device.id, date.time.local) %>%
  group_by(device.id) %>%
  mutate(regular = dt_hrs > (2 - 0.33) & dt_hrs < (2 + 0.33)
  ) %>%
  ungroup()

table(tiger_data$regular)
prop.table(table(tiger_data$regular))



# Speed filter -----------------------------------------------------------------
# Will identify and remove ~1% of points regardless of whether the animal is a tiger, tortoise, or eagle.
# Does not account for regularity of fix interval
# Easy, biologically interpretable, accounts for fix interval but long gaps will deflate speed estimate
# Does not account for telemetry error but potentially removes fixes owing to telemetry error

cat("\n================ SPEED FILTER ================\n")

summary(tiger_data$speed_kmh)

ggplot(tiger_data, aes(speed_kmh)) +
  geom_histogram(binwidth = 0.5,
                 fill = "steelblue",
                 color = "white") +
  geom_vline(
    xintercept = quantile(tiger_data$speed_kmh,
                          0.99,
                          na.rm = TRUE),
    colour = "red",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = quantile(tiger_data$speed_kmh,
                          0.95,
                          na.rm = TRUE),
    colour = "blue",
    linewidth = 1
  ) +
  theme_bw() +
  labs(
    title = "Speed Distribution",
    x = "Speed (km/h)",
    y = "Count"
  )

speed_threshold <- quantile(
  tiger_data$speed_kmh,
  0.99,
  na.rm = TRUE
)

speed_outliers <- tiger_data %>%
  filter(speed_kmh > speed_threshold)

cat("Speed threshold:", round(speed_threshold,2), "km/h\n")
cat("Flagged locations:", nrow(speed_outliers), "\n")

# This will filter out locations beyond speed threshold
speed_data <- tiger_data %>%
  filter(speed_kmh <= speed_threshold)


## Question ##
cat("\nQUESTION: What are potential issues with filtering simply based on speed between two
successive locations?\n")


# Speed outlier identification and review --------------------------------------
# Calculate speed threshold and flag extreme movement events.
# These are flagged only (not automatically removed) because a high-speed
# movement can be caused by either the preceding or following GPS fix.
# Flagged locations should be visually inspected before removal.

# Inspect speed distribution
summary(tiger_data$speed_kmh)

# Flag locations associated with extreme movement speeds
tiger_data <- tiger_data %>%
  mutate(
    speed_flag = speed_kmh > speed_threshold
  )

# Summary of flagged locations
cat(
  "Flagged locations:",
  sum(tiger_data$speed_flag, na.rm = TRUE),
  "\n"
)


# Extract flagged locations for inspection
speed_outliers <- tiger_data %>%
  filter(speed_flag) %>%
  select(
    device.id,
    date.time.local,
    longitude,
    latitude,
    speed_kmh
  )

speed_outliers


# Visual inspection of flagged locations
ggplot(tiger_data, aes(longitude, latitude)) +
  geom_path(
    aes(group = device.id),
    linewidth = 0.5
  ) +
  geom_point(
    data = filter(tiger_data, speed_flag),
    colour = "red",
    size = 3
  ) +
  theme_bw() +
  labs(
    title = "Flagged High-Speed GPS Locations",
    x = "Longitude",
    y = "Latitude"
  )


# After inspection, remove only confirmed erroneous locations:
# tiger_data <- tiger_data %>%
#   filter(!speed_flag)

# Do not automatically remove all flagged locations without checking,
# because the high-speed step may be caused by the previous fix rather
# than the flagged location itself.


# Create ctmm object -----------------------------------------------------------
library(ctmm)
# Prepare data for ctmm
gps_df <- tiger_data %>%
  dplyr::select(
    device.id,
    date.time.local,
    longitude,
    latitude, 
    dop
  ) %>%
  rename(
    ID = device.id,
    timestamp = date.time.local
  ) %>%
  as.data.frame()

# Create ctmm telemetry object
gps_ctmm <- as.telemetry(
  gps_df,
  timezone = "Asia/Thimphu",
  na.rm = "col",
  keep = TRUE,
  drop = FALSE
)

# Check
class(gps_ctmm[[1]])


# Identify and remove extreme movement outliers using ctmm ---------------------
cat("\nNOTE: ctmm::outlie()` identifies potential outliers by calculating **distance-based
anomalies** (locations unusually far from the animal’s typical movement
distribution) and **speed-based anomalies** (locations requiring implausibly
high movement speeds between fixes, accounting for timestamp precision and
location error), which are then flagged for inspection rather than
automatically removed. Tries to determine the most likely offending locations.
\n")

cat("\nNOTE: Fixed thresholds (e.g., speed > 1 m/s) are species-specific and should
# not be applied universally. A threshold suitable for one species may remove
# valid movements in another.\n")


# Run outlier detection on one individual
OUT <- outlie(gps_ctmm[[1]])

# Visual check of remaining movement pattern
plot(
  OUT,
  units = FALSE
)

cat("\nNOTE: Thicker bluer lines indicate higher speed movements
'Median deviation’ denotes distances from the geometric median, while ‘minimum
speed’ denotes the minimum speed required to explain the location estimate's
displacement as straight-line motion. Both estimates account for telemetry
error and condition on as few data points as possible. The speed estimates
furthermore account for timestamp truncation and assign each timestep's speed
to the most likely offending time, based on its other adjacent speed estimate.\n")

# Conversion of speeds from m/s to km/h for reference
speed_table <- data.frame(
  speed_m_s = seq(0.5, 5, by = 0.5),
  speed_km_h = seq(0.5, 5, by = 0.5) * 3.6
)

speed_table

# Define biological speed threshold
# Example: remove fixes associated with speeds > 1 m/s
speed_threshold <- 1.5

# Identify potential speed outliers
BAD <- OUT$speed > speed_threshold

# Check number of flagged fixes
table(BAD)

# Remove flagged locations from telemetry object
gps_ctmm[[1]] <- gps_ctmm[[1]][!BAD, ]

# Re-run outlier detection to confirm improvement
OUT_clean <- outlie(gps_ctmm[[1]])

plot(
  OUT_clean,
  units = FALSE
)

# Inspect movement metrics
summary(OUT)


# Identify extreme movements --------------------------------------------------
# Flag the top 1% of movements based on displacement distance
# These may represent GPS errors, unrealistic jumps, or data artefacts

BAD <- OUT$distance > quantile(
  OUT$distance,
  0.99,
  na.rm = TRUE
)

# Check number of flagged fixes
table(BAD)


# # Remove flagged fixes from telemetry object
# gps_ctmm[[1]] <- gps_ctmm[[1]][!BAD,]
# 
# 
# # Re-run outlier detection after cleaning
# OUT <- outlie(gps_ctmm[[1]])
# 
# # Visual check of remaining movement pattern
# plot(
#   OUT,
#   units = FALSE
# )


# Optional: identify extreme values using both speed and distance -------------
# This provides a stricter filter by removing fixes that are extreme in either
# displacement OR movement rate

BAD_CTMM <- OUT$speed > quantile(
  OUT$speed,
  0.99,
  na.rm = TRUE
) |
  OUT$distance > quantile(
    OUT$distance,
    0.99,
    na.rm = TRUE
  )

table(BAD_CTMM)





# High speed and sharp turn (v-spike) filter -----------------------------------
calc_bearing <- function(x1,y1,x2,y2){
  atan2(y2-y1,x2-x1)
}

gps_turns <- tiger_data %>%
  arrange(device.id, date.time.local) %>%
  group_by(device.id) %>%
  mutate(
    
    bearing_in =
      calc_bearing(
        lag(x),
        lag(y),
        x,
        y
      ),
    
    bearing_out =
      calc_bearing(
        x,
        y,
        lead(x),
        lead(y)
      ),
    
    turn_angle =
      bearing_out - bearing_in,
    
    turn_angle =
      ifelse(
        turn_angle > pi,
        turn_angle - 2*pi,
        turn_angle
      ),
    
    turn_angle =
      ifelse(
        turn_angle < -pi,
        turn_angle + 2*pi,
        turn_angle
      ),
    
    turn_deg =
      abs(turn_angle * 180/pi)
    
  ) %>%
  ungroup()

# Plot
ggplot(
  gps_turns,
  aes(speed_kmh, turn_deg)
) +
  geom_point(alpha = 0.3) +
  geom_hline(
    yintercept = 150,
    colour = "red",
    linetype = 2
  ) +
  theme_bw() +
  labs(
    title = "Speed vs Turning Angle",
    x = "Speed (km/h)",
    y = "Turning Angle (degrees)"
  )

v_spikes <- gps_turns %>%
  filter(
    speed_kmh > quantile(speed_kmh,
                         0.95,
                         na.rm = TRUE),
    turn_deg > 150
  )

cat("\n================ V-SPIKES ================\n")

cat(
  "Flagged locations:",
  nrow(v_spikes),
  "\n"
)

# Plot
ggplot() +
  geom_path(
    data = gps_turns,
    aes(x = x, y = y),
    colour = "grey80"
  ) +
  geom_point(
    data = v_spikes,
    aes(x = x, y = y),
    colour = "red",
    size = 2
  ) +
  coord_equal() +
  theme_bw() +
  labs(
    title = "Potential GPS V-Spikes",
    subtitle = "Red points are high-speed sharp turns"
  )

# SUMMARY OF ALL METHODS
summary_table <- data.frame(
  Method = c(
    "Speed",
    "Distance",
    "CTMM",
    "V-Spike"
  ),
  Flagged_Locations = c(
    nrow(speed_outliers),
    nrow(distance_outliers),
    sum(BAD_CTMM, na.rm = TRUE),
    nrow(v_spikes)
  )
)

print(summary_table)


# Original data
# 3778 locations
nrow(tiger_data)
# 3778-1415 = 2363



# High speed spikes, out and back ----------------------------------------------
gps_analysis <- tiger_data %>%
  arrange(device.id, date.time.local) %>%
  group_by(device.id) %>%
  mutate(
    # Existing 'Speed In' (speed to get to this point)
    # Ensure this is already calculated in your tiger_data
    speed_in = speed_kmh, 
    
    # 'Speed Out' (speed to get to the next point)
    speed_out = lead(speed_kmh),
    
    # Calculate bearings and turning angles
    bearing_in   = calc_bearing(lag(x), lag(y), x, y),
    bearing_out  = calc_bearing(x, y, lead(x), lead(y)),
    turn_angle   = bearing_out - bearing_in,
    turn_angle   = ifelse(turn_angle > pi, turn_angle - 2*pi, turn_angle),
    turn_angle   = ifelse(turn_angle < -pi, turn_angle + 2*pi, turn_angle),
    turn_deg     = abs(turn_angle * (180 / pi))
  ) %>%
  ungroup()

# --- DEFINE DUAL THRESHOLDS ---
min_speed_thresh <- 3.5  # km/h
min_angle_thresh <- 145  # degrees

# Flagging points where BOTH legs are fast AND the turn is sharp
gps_analysis <- gps_analysis %>%
  mutate(is_spike = case_when(
    speed_in > min_speed_thresh & 
      speed_out > min_speed_thresh & 
      turn_deg > min_angle_thresh ~ "YES",
    TRUE ~ "NO"
  ))

# Review the spikes
v_spikes <- gps_analysis %>% filter(is_spike == "YES")
nrow(v_spikes)

ggplot(gps_analysis, aes(x = speed_in, y = speed_out, color = turn_deg)) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_c() +
  geom_vline(xintercept = min_speed_thresh, linetype = "dashed", color = "red") +
  geom_hline(yintercept = min_speed_thresh, linetype = "dashed", color = "red") +
  labs(title = "Out-and-Back Speed Analysis",
       x = "Speed In (km/h)",
       y = "Speed Out (km/h)",
       color = "Turn Angle") +
  theme_minimal()

# Update the suspicious_points dataframe using the new flags
suspicious_points <- gps_analysis %>%
  filter(is_spike == "YES")

# Map the path and highlight the "Double-Speed" spikes
ggplot() +
  # Draw the full path in grey
  geom_path(data = gps_analysis, aes(x = x, y = y), 
            color = "grey80", alpha = 0.5) +
  
  # Highlight the suspicious "V" apex points in red
  geom_point(data = suspicious_points, aes(x = x, y = y), 
             color = "red", size = 2.5) +
  
  # Label with Speed In / Speed Out
  geom_text(data = suspicious_points, 
            aes(x = x, y = y, 
                label = paste0("In:", round(speed_in, 1), 
                               "\nOut:", round(speed_out, 1))), 
            vjust = -1.2, size = 2.5, color = "darkred", fontface = "bold") +
  
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Map of Suspicious GPS 'V-Spikes'",
    subtitle = paste0("Red points: Speed In/Out > ", min_speed_thresh, 
                      " km/h & Turn > ", min_angle_thresh, "°"),
    x = "Easting (m)", 
    y = "Northing (m)"
  )


# Assessing space use shifts ---------------------------------------------------
# Temporal window
days_per_window <- 7  # Change to 14 or 30 as needed

# Prepare data with the window index
tiger_data <- tiger_data %>%
  arrange(date.time.local) %>%
  mutate(
    days_elapsed = as.numeric(difftime(date.time.local, min(date.time.local), units = "days")),
    time_window = floor(days_elapsed / days_per_window) + 1
  )

# Convert to sf for leaflet
gps_sf_ll <- tiger_data %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

# Change tmeporal window to fcator for more distinct visulaisation
gps_sf_ll$time_window_factor <- as.factor(gps_sf_ll$time_window)

# Create a palette based on the discrete windows
n_windows <- length(unique(gps_sf_ll$time_window))

extended_pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_windows)

pal <- colorFactor(
  palette = extended_pal,
  domain = gps_sf_ll$time_window_factor
)

# Interactive Plot
leaflet(gps_sf_ll) %>%
  # addProviderTiles(providers$Esri.WorldImagery) %>% # Satellite view is often better for ecology
  addTiles() %>%              # Optional: adds toggleable street map
  addCircleMarkers(
    radius = 4,
    color = ~pal(time_window_factor),
    stroke = TRUE,
    weight = 1,
    fillOpacity = 0.9,
    # Updated label to show which window the point belongs to
    label = lapply(paste0(
      "<b>Window: </b>", gps_sf_ll$time_window, "<br>",
      "<b>Time: </b>", format(gps_sf_ll$date.time.local, "%Y-%m-%d %H:%M"), "<br>"
    ), htmltools::HTML),
    labelOptions = labelOptions(direction = "auto")
  ) %>%
  addLegend(
    "bottomright",
    pal    = pal,
    values = ~time_window,
    title  = paste0("Time Window (", days_per_window, " days)"),
    opacity = 1
  ) %>%
  addLayersControl(
    baseGroups = c("Satellite", "Standard Map"),
    options = layersControlOptions(collapsed = FALSE)
  )


# Define which window numbers belong in which group
gps_sf_ll <- gps_sf_ll %>%
  mutate(custom_phase = case_when(
    time_window <= 19 ~ "First Month",
    time_window > 19 & time_window <= 24 ~ "Second Month",
    time_window > 24 ~ "Final Period"
  ))

# Re-factor to keep them in order in the legend
gps_sf_ll$custom_phase <- factor(gps_sf_ll$custom_phase, 
                                 levels = c("First Month", "Second Month", "Final Period"))

pal_custom <- colorFactor(palette = "Dark2", domain = gps_sf_ll$custom_phase)

# Interactive Plot
leaflet(gps_sf_ll) %>%
  addTiles() %>%
  addCircleMarkers(
    radius = 4,
    color = ~pal_custom(custom_phase),
    fillOpacity = 0.8,
    stroke = FALSE,
    label = ~paste0("Phase: ", custom_phase, " | Time: ", format(date.time.local, "%b %d"))
  ) %>%
  addLegend("bottomright", pal = pal_custom, values = ~custom_phase, title = "Temporal Phases")

# Not cleaned data!!!





