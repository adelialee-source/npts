library(ggplot2)
library(sf)
library(dplyr)
library(ggmap)
register_stadiamaps("API-KEY") 

# Getting the location from the data
reshape_data_loc <- function(data) {
  data %>%
    subset(State.Code == 8 & Method.Code != 236 & Method.Code != 238) %>%
    select(Latitude, Longitude)
}

# Load data
AQ_data_PM <- read.csv('daily_88101_2020.csv')
AQ_data_PM_train <- read.csv('daily_88101_2019.csv')
AQ_data_PM_train2 <- read.csv('daily_88101_2018.csv')
AQ_data_PM_train3 <- read.csv('daily_88101_2017.csv')
AQ_data_PM_train4 <- read.csv('daily_88101_2016.csv')
AQ_data_PM_train5 <- read.csv('daily_88101_2015.csv')

reshaped_PM_loc <- reshape_data_loc(AQ_data_PM)
reshaped_PM_train_loc <- reshape_data_loc(AQ_data_PM_train)
reshaped_PM_train2_loc <- reshape_data_loc(AQ_data_PM_train2)
reshaped_PM_train3_loc <- reshape_data_loc(AQ_data_PM_train3)
reshaped_PM_train4_loc <- reshape_data_loc(AQ_data_PM_train4)
reshaped_PM_train5_loc <- reshape_data_loc(AQ_data_PM_train5)


Colorado <- rbind(reshaped_PM_loc, reshaped_PM_train_loc, reshaped_PM_train2_loc, reshaped_PM_train3_loc, reshaped_PM_train4_loc, reshaped_PM_train5_loc)
Colorado <- unique(Colorado)

location_data <- as.data.frame(Colorado)

# Coordinates of the Cameron Peak Fire in 2020 (approximate coordinates)
fire_spot <- data.frame(
  Latitude = 40.4805,
  Longitude = -105.4375,
  value = 100
)

fire_spot$Longitude <- as.numeric(fire_spot$Longitude)
fire_spot$Latitude <- as.numeric(fire_spot$Latitude)

location_data$Longitude <- as.numeric(location_data$Longitude)
location_data$Latitude <- as.numeric(location_data$Latitude)


# Plot
register_stadiamaps(key = "your key code here")

# Set the center location for Colorado
center <- c(lon = -105.5, lat = 39.5)

# Range
bbox <- c(left = -109.5, bottom = 36.5, right = -102, top = 41.5) 
map <- get_stadiamap(bbox = bbox, zoom = 7, maptype = "stamen_toner")

# Plot the map
ggmap(map) +
  geom_point(data = fire_spot, aes(x = Longitude, y = Latitude), color = '#ff7c43', size = 15, shape = 13, stroke = 2) + 
  geom_point(data = location_data, aes(x = Longitude, y = Latitude), color = '#2f4b7c', size = 3) +
  theme_minimal() +
  labs(x = 'Longitude', y = 'Latitude') +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )


