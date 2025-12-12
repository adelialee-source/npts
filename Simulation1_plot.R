# Packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(zoo)
library(janitor)

# Load data
AQ_data_PM <- read.csv('daily_88101_2020.csv')
AQ_data_PM_train <- read.csv('daily_88101_2019.csv')
AQ_data_PM_train2 <- read.csv('daily_88101_2018.csv')
AQ_data_PM_train3 <- read.csv('daily_88101_2017.csv')
AQ_data_PM_train4 <- read.csv('daily_88101_2016.csv')
AQ_data_PM_train5 <- read.csv('daily_88101_2015.csv')

step_1 <- function(data) {
  data %>%
    subset(State.Code == 8 & Method.Code != 236 & Method.Code != 238) %>%
    group_by(County.Code, Date.Local) %>%
    select(County.Code, Date.Local, Arithmetic.Mean) %>%
    mutate(Row = row_number()) %>%
    pivot_wider(names_from = Row, values_from = Arithmetic.Mean) %>%
    rename_with(~ paste0("Obs_", .), starts_with("1") | starts_with("2") | starts_with("3") | starts_with("4") | starts_with("5") | starts_with("6") | starts_with("7") | starts_with("8") | starts_with("9")) %>%
    arrange(County.Code, Date.Local) %>%
    ungroup()
}

# Step 2: Identify the longest row in each County
length_row <- function(data) {
  data %>%
    mutate(length = rowSums(!is.na(select(., starts_with("Obs_")))))
}

longest_row <- function(data){
  data %>% 
    group_by(County.Code) %>%
    mutate(max = max(length))
}

control_max <- function(data){
  Test <- data
  Test <- Test %>%
    mutate(max = case_when(
      County.Code == 1 ~ 3,
      County.Code == 5 ~ 1,
      County.Code == 13 ~ 4,
      County.Code == 31 ~ 14,
      County.Code == 35 ~ 3,
      County.Code == 41 ~ 1,
      County.Code == 45 ~ 4,
      County.Code == 69 ~ 2,
      County.Code == 77 ~ 2,
      County.Code == 101 ~ 1,
      County.Code == 103 ~ 2,
      County.Code == 123 ~ 3,
      TRUE ~ 0
    ))
  return(Test)
}

step_2 <- function(data){
  data <- control_max(data)
  data <- data[data$max != 0,]
  if (ncol(data) < 16){
    data <- cbind(data, NA)
  }
  for(j in 1:nrow(data)){
    max <- data$max[j]
      for(i in 3:(max + 2)){
        if (is.na(data[j, i])) {
          data[j, i] <- mean(as.numeric(data[j, 3:(i-1)]), na.rm = TRUE)
        }
      }
    }
  return(data)
}

step_3 <- function(data){
  data <- data %>% arrange(Date.Local)
  Date <- unique(data$Date.Local)
  County <- unique(data$County.Code)
  County <- County[County != 93]
  Date <- as.matrix(Date)
  current_data <- Date
  
  for(county in County) {
    county_data <- data[data$County.Code == county,]
    county_all_dates <- data.frame(County.Code = county, Date.Local = Date)
    
    merged_data <- merge(county_all_dates, county_data, by = c("County.Code", "Date.Local"), all.x = TRUE)
    merged_data[,-c(1, 2)] <- zoo::na.locf(merged_data[,-c(1, 2)], na.rm = FALSE)
    
    current_data <- cbind(current_data, merged_data[3:(ncol(merged_data)-2)])
  }
  
  clean_data <- current_data %>% remove_empty('cols')
  clean_data[, -1] <- apply(clean_data[, -1], 1, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  })
  return(clean_data)
}


redesign <- function(data){
  Test_1 <- step_1(data)
  Test_2 <- step_2(Test_1)
  Test_3 <- step_3(Test_2)
  
  return(Test_3)
}

redesign_2020 <- redesign(AQ_data_PM)
redesign_2019 <- redesign(AQ_data_PM_train)
redesign_2018 <- redesign(AQ_data_PM_train2)
redesign_2017 <- redesign(AQ_data_PM_train3)
redesign_2016 <- redesign(AQ_data_PM_train4)
redesign_2015 <- redesign(AQ_data_PM_train5)

# Create training data
row_mean_2017 <- rowMeans(redesign_2017[,2:ncol(redesign_2017)], na.rm = T)
row_mean_2017 <- as.matrix(row_mean_2017, ncol = 1)
redesign_2017 <- cbind(redesign_2017, row_mean_2017)
row_mean_2015 <- rowMeans(redesign_2015[,2:ncol(redesign_2015)], na.rm = T)
row_mean_2015 <- as.matrix(row_mean_2015, ncol = 1)
redesign_2015 <- cbind(redesign_2015, row_mean_2015, row_mean_2015, row_mean_2015)
redesign_2015 <- redesign_2015[,1:41]

# Erase 2/29 data to consistency (2016)
k <- which(redesign_2016$current_data == '2016-02-29') #243
redesign_2016 <- redesign_2016[-k,]
redesign_2016 <- redesign_2016[, 1:41]

# Format training data
names(redesign_2015) <- names(redesign_2016) <- names(redesign_2017) <- names(redesign_2018)
training <- rbind(redesign_2015, redesign_2016, redesign_2017, redesign_2018)

# Format for the average data
training$current_data <- as.Date(training$current_data)
training$current_data <- format(training$current_data, format = "%m-%d")

days_avg <- training %>%
  group_by(current_data) %>%
  summarise(across(starts_with('Obs'), ~ mean(. , na.rm = TRUE), .names = "avg_{.col}"))

# Smoothing data
window_size <- 30
days_avg[,-1] <- zoo::rollmean(days_avg[,-1], k = window_size, fill = 8, align = "center")

# Time series
running <- redesign_2020[-k,-1] - days_avg[,-1]
train <- redesign_2019[,-1] - days_avg[,-1]
names(train) <- names(running)
M_t0 <- rbind(train, running)
M_t0 <- as.matrix(M_t0)

# Initial settings
Date <- redesign_2020[-k,1]
training_m <- 483
total_T <- 730 - training_m
thres <- 5.85

# Running test statistics
M_t <- array(0, dim = c(ncol(M_t0), ncol(M_t0), training_m + total_T))
lambda <- rep(0, training_m + total_T)
Gamma_k <- matrix(0, nrow = training_m + total_T, ncol = 1)

for (t in 1:(training_m + total_T)) {
  # Test Statistics
  M_t[,,t] <- M_t0[t,] %*% t(M_t0[t,])
  lambda[t] <- max(eigen(M_t[,,t])$values)
  
  if (t > training_m) {
    # Find nominator V_m only depending on m
    V_m <- 0
    
    for (s in 1:training_m){
      partial <- abs(sum(lambda[1:s]) - (s / training_m) * sum(lambda[1:training_m]))
      V_m <- V_m + partial
    }
    
    # Find denominator depending on k
    sup_D <- 0
    
    for (k in 1:(t - training_m)){
      D_mk <- ( 1 / (training_m+k) ) * (sum(lambda[(training_m+1):(training_m+k)]) - (k / training_m) * sum(lambda[1:training_m]))
      sup_D <- ifelse(sup_D > D_mk, sup_D, D_mk)
    }
    
    Gamma_k[t, 1] <- training_m^2 * (sup_D / V_m)
  }
}


for (t in 1: (training_m + total_T)){
  if(Gamma_k[t, 1] >= thres){
    print(Date[t - 365])
    print(t - 365)
    break
  }
}

Date <- format(Date, format = "%Y - %m-%d")
Date <- as.Date(Date)
data <- data.frame(y = Gamma_k[366:730,1], x = Date)
threshold_crossing <- data$x[t - 365]
monitoring <- data$x[training_m - 366]
real <- data$x[t - 365 - 7]

# Plot
ggplot(data, aes(x = x, y = y)) +
  geom_point(aes(color = ifelse(y > thres, '#ff7c43', '#2f4b7c')), size = 1) +
  # Add a horizontal line at the threshold
  geom_hline(yintercept = thres, linetype = "dashed", color = "#12436D") +
  # Add a vertical line at the threshold crossing point
  geom_vline(xintercept = threshold_crossing, linetype = "dashed", color = "#12436D") +
  geom_vline(xintercept = monitoring, linetype = "dashed", color = "#12436D") +
  geom_vline(xintercept = real, linetype = "solid", color = "#12436D") +
  annotate("text", x = threshold_crossing, y = thres, label = paste("Detection on", threshold_crossing), 
           vjust = 1.6, hjust = -0.1, color = "#12436D") +
  annotate("text", x = monitoring, y = 0, label = paste("Training until", monitoring), 
           vjust = -0.6, hjust = 1.1, color = "#12436D") +
  scale_color_identity() + 
  labs(x = "Day", y = "Test statistic") +
  scale_x_date(
    breaks = seq.Date(from = as.Date('2020-01-01'), to = as.Date('2020-12-31'), by = 'month'),  
    # 1-month interval
    labels = month.abb) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  theme_classic() 

