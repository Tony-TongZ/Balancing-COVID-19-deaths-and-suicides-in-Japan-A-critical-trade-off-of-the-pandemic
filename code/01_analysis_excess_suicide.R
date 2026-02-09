############################################################
# Analysis of excess suicides by sex and age group
# This script reproduces all numerical results used
# in the manuscript and supplementary materials.
# Exploratory analyses and visualization are omitted.
############################################################

## 0. Setup -------------------------------------------------
set.seed(123)

library(forecast)
library(dplyr)
library(tidyr)

## 1. Data input -------------------------------------------
# Input data files should be placed in the data/ directory

maledf <- read.csv("data/suicide_male.csv", header = TRUE)
femaledf <- read.csv("data/suicide_female.csv", header = TRUE)

## Aggregate age groups (5-year to 10-year groups)
maledf_new <- data.frame(
  au20 = maledf$a10 + maledf$a15,
  a20  = maledf$a20 + maledf$a25,
  a30  = maledf$a30 + maledf$a35,
  a40  = maledf$a40 + maledf$a45,
  a50  = maledf$a50 + maledf$a55,
  a60  = maledf$a60 + maledf$a65,
  aa70 = maledf$a70 + maledf$a75 + maledf$a80 +
    maledf$a85 + maledf$a90 + maledf$a95 + maledf$a100
)
rownames(maledf_new) <- rownames(maledf)

femaledf_new <- data.frame(
  au20 = femaledf$a10 + femaledf$a15,
  a20  = femaledf$a20 + femaledf$a25,
  a30  = femaledf$a30 + femaledf$a35,
  a40  = femaledf$a40 + femaledf$a45,
  a50  = femaledf$a50 + femaledf$a55,
  a60  = femaledf$a60 + femaledf$a65,
  aa70 = femaledf$a70 + femaledf$a75 + femaledf$a80 +
    femaledf$a85 + femaledf$a90 + femaledf$a95 + femaledf$a100
)
rownames(femaledf_new) <- rownames(femaledf)

## 2. SARIMA function --------------------------------------
# Estimate baseline suicides and prediction intervals

fit_sarima_excess <- function(unem, train_len = 60, h = 24) {
  tsunem <- ts(unem[1:train_len], frequency = 12)
  dummy  <- c(rep(0, 36), rep(1, 24))
  
  model <- auto.arima(
    tsunem,
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

## 3. Object initialization --------------------------------
init_baseline <- function(n_time) {
  data.frame(
    au20 = rep(0, n_time),
    a20  = rep(0, n_time),
    a30  = rep(0, n_time),
    a40  = rep(0, n_time),
    a50  = rep(0, n_time),
    a60  = rep(0, n_time),
    aa70 = rep(0, n_time)
  )
}

# Baseline (2015–2020: 66 months)
suibase_male        <- init_baseline(66)
suibase_male_upper  <- init_baseline(66)
suibase_male_lower  <- init_baseline(66)
suibase_female      <- init_baseline(66)
suibase_female_upper<- init_baseline(66)
suibase_female_lower<- init_baseline(66)

# Excess suicides (2020–2021: 24 months)
exsui_male   <- init_baseline(24)
exsui_female <- init_baseline(24)
sui_male_obs <- init_baseline(24)
sui_male_pre <- init_baseline(24)
sui_female_obs <- init_baseline(24)
sui_female_pre <- init_baseline(24)

## 4. Estimation loop --------------------------------------
# Loop over age groups to estimate baseline and excess suicides

for (i in seq_len(ncol(maledf_new))) {
  ## Male
  unem_m <- maledf_new[, i]
  res_m  <- fit_sarima_excess(unem_m)
  
  suibase_male[, i]       <- res_m$pre[1:66]
  suibase_male_upper[, i] <- res_m$upper[1:66]
  suibase_male_lower[, i] <- res_m$lower[1:66]
  exsui_male[, i]         <- unem_m[61:84] - res_m$pre[61:84]
  sui_male_obs[, i] <- unem_m[61:84]
  sui_male_pre[, i] <- res_m$pre[61:84]
  
  ## Female
  unem_f <- femaledf_new[, i]
  res_f  <- fit_sarima_excess(unem_f)
  
  suibase_female[, i]       <- res_f$pre[1:66]
  suibase_female_upper[, i] <- res_f$upper[1:66]
  suibase_female_lower[, i] <- res_f$lower[1:66]
  exsui_female[, i]         <- unem_f[61:84] - res_f$pre[61:84]
  sui_female_obs[, i] <- unem_f[61:84]
  sui_female_pre[, i] <- res_f$pre[61:84]
}

## 5. Final objects ----------------------------------------
# Objects used in the manuscript

suibase       <- cbind(suibase_male, suibase_female)
suibase_upper <- cbind(suibase_male_upper, suibase_female_upper)
suibase_lower <- cbind(suibase_male_lower, suibase_female_lower)

suiobs <- cbind(maledf_new, femaledf_new)

## 6. Save results (optional) -------------------------------
# save(
#   suibase, suibase_upper, suibase_lower,
#   exsui_male, exsui_female, suiobs,
#   file = "results/excess_suicide_results.RData"
# )

############################################################
# End of script
############################################################