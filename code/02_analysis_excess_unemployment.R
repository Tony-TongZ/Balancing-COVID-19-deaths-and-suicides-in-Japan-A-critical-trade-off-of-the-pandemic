############################################################
# Analysis of unemployment by sex and age group
# This script reproduces all numerical results
# used in the manuscript.
# Visualization code is intentionally omitted.
############################################################

## 0. Setup -------------------------------------------------
set.seed(123)

library(forecast)
library(dplyr)
library(tidyr)

## 1. Data input -------------------------------------------
# Input files should be placed in the data/ directory

maledf <- read.csv("data/unemployment_male.csv", header = TRUE)
femaledf <- read.csv("data/unemployment_female.csv", header = TRUE)

## 2. Age-group interpolation ------------------------------
# Original age groups: 15–24, 25–34, ..., 65–69
# Target age groups:   15–19, 20–29, ..., 60–69

spline_interpolation <- function(age, value, new_age) {
  spline_fit <- splinefun(age, value, method = "natural")
  round(spline_fit(new_age), 2)
}

midage <- c((15+24)/2, (25+34)/2, (35+44)/2,
            (45+54)/2, (55+64)/2, (65+69)/2)

midage_new <- c((15+19)/2, (20+29)/2, (30+39)/2,
                (40+49)/2, (50+59)/2, (60+69)/2)

## Construct interpolated datasets
build_interpolated_df <- function(df) {
  out <- data.frame(
    a10 = rep(0, nrow(df)),
    a20 = rep(0, nrow(df)),
    a30 = rep(0, nrow(df)),
    a40 = rep(0, nrow(df)),
    a50 = rep(0, nrow(df)),
    a60 = rep(0, nrow(df))
  )
  
  for (i in seq_len(nrow(out))) {
    out[i, ] <- spline_interpolation(midage, df[i, 2:7], midage_new)
  }
  out
}

maledf_new   <- build_interpolated_df(maledf)
femaledf_new <- build_interpolated_df(femaledf)

## 3. SARIMA estimation function ----------------------------

fit_sarima_excess <- function(series, train_len = 60, h = 24) {
  tsdata <- ts(series[1:train_len], frequency = 12)
  dummy  <- c(rep(0, 36), rep(1, 24))
  
  model <- auto.arima(
    tsdata,
    d = NA, D = NA,
    start.p = 0, start.q = 0,
    start.P = 0, start.Q = 0,
    seasonal = TRUE,
    stepwise = FALSE,
    xreg = dummy
  )
  
  fc <- forecast(
    model,
    h = h,
    level = 95,
    xreg = rep(1, h)
  )
  
  list(
    pre   = c(model$fitted, fc$mean),
    lower = c(fitted(model) - 1.96 * sqrt(model$sigma2), fc$lower),
    upper = c(fitted(model) + 1.96 * sqrt(model$sigma2), fc$upper)
  )
}

## 4. Initialize result objects -----------------------------

init_excess <- function() {
  data.frame(
    a10 = rep(0, 24),
    a20 = rep(0, 24),
    a30 = rep(0, 24),
    a40 = rep(0, 24),
    a50 = rep(0, 24),
    a60 = rep(0, 24)
  )
}

exun_male   <- init_excess()
exun_female <- init_excess()

## 5. Estimation loop ---------------------------------------

for (i in seq_len(ncol(exun_male))) {
  ## Male
  series_m <- maledf_new[, i]
  res_m <- fit_sarima_excess(series_m)
  exun_male[, i] <- series_m[61:84] - res_m$pre[61:84]
  
  ## Female
  series_f <- femaledf_new[, i]
  res_f <- fit_sarima_excess(series_f)
  exun_female[, i] <- series_f[61:84] - res_f$pre[61:84]
}

## 6. Final objects -----------------------------------------
# Objects used in the manuscript:
# - exun_male
# - exun_female

## 7. Save results (optional) -------------------------------
# save(
#   exun_male, exun_female,
#   file = "results/excess_unemployment_results.RData"
# )

############################################################
# End of script
############################################################