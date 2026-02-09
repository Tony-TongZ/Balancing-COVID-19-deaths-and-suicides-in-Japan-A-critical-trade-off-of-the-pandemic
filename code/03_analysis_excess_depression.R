############################################
## Social Mood Index (Depression)
## Weekly to Monthly Aggregation
## Excess Depression Estimation
############################################

## ---- 0. Setup ----
set.seed(123)

library(forecast)
library(dplyr)
library(lubridate)

## ---- 1. Read data ----
depression_weekly <- read.csv("data/Social_mood_index_depression.csv",header = FALSE)

mood_weekly <- data.frame(
  CalendarWeek = 1:nrow(depression_weekly),
  depression   = depression_weekly$V2
)

## ---- 2. Convert weekly data to monthly ----
start_date <- as.Date("2015-01-05")  # start of calendar week 1 in 2015

monthly_data <- data.frame(
  Year = integer(),
  Month = integer(),
  depression = numeric()
)

for (i in 1:nrow(mood_weekly)) {
  
  week_start <- start_date + weeks(i - 1)
  week_end   <- week_start + days(6)
  total_days <- as.numeric(week_end - week_start) + 1
  
  current_date <- week_start
  
  while (current_date <= week_end) {
    
    year_now  <- year(current_date)
    month_now <- month(current_date)
    
    days_in_current_month <- days_in_month(current_date) - day(current_date) + 1
    days_used <- min(
      as.numeric(week_end - current_date) + 1,
      days_in_current_month
    )
    
    allocated_value <- mood_weekly$depression[i] * (days_used / total_days)
    
    monthly_data <- monthly_data %>%
      add_row(
        Year = year_now,
        Month = month_now,
        depression = allocated_value
      )
    
    current_date <- current_date + days(days_used)
  }
}

mood_monthly <- monthly_data %>%
  group_by(Year, Month) %>%
  summarise(depression = sum(depression), .groups = "drop") %>%
  mutate(Time = 1:n())

## ---- 3. SARIMA model (pre-pandemic baseline) ----
depression_ts <- ts(
  mood_monthly$depression[1:60],  # Jan 2015 – Dec 2019
  frequency = 12
)

model_sarima <- auto.arima(
  depression_ts,
  d = NA,
  D = NA,
  start.p = 0,
  start.q = 0,
  start.P = 0,
  start.Q = 0,
  seasonal = TRUE,
  stepwise = FALSE
)

forecast_sarima <- forecast(
  model_sarima,
  h = 24,
  level = 95
)

## ---- 4. Counterfactual & excess depression ----
baseline_estimate <- c(
  fitted(model_sarima),
  forecast_sarima$mean
)

mood_monthly <- mood_monthly %>%
  mutate(
    counterfactual = baseline_estimate,
    excess_depression = depression - counterfactual
  )

excess_depression <- mood_monthly$excess_depression[61:84]

## ---- Output ----
excess_depression

############################################################
# End of script
############################################################