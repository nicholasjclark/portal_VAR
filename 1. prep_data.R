# R script originally prepared by Nicholas J Clark (nicholas.j.clark1214@gmail.com)
#### 1. Prep data for modelling ####
library(portalcasting)
library(mgcv)
library(dplyr)

# Load the most recent rodents survey table for control plots
rodents_table <- read.csv('data/rodents_data.csv', as.is = TRUE)

# Calculate means and sds of covariates for later unscaled plotting
rodents_table %>%
  dplyr::mutate(month = lubridate::month(newmoondate),
                year = lubridate::year(newmoondate)) %>%
  dplyr::filter(year > 1995) %>%
  dplyr::select(mintemp) %>%
  dplyr::summarise(mintemp_mean = mean(mintemp, na.rm = TRUE),
                   mintemp_sd = sd(mintemp, na.rm = TRUE)) -> mintemp_stats

rodents_table %>%
  dplyr::mutate(month = lubridate::month(newmoondate),
                year = lubridate::year(newmoondate)) %>%
  dplyr::filter(year > 1995) %>%
  dplyr::select(ndvi) %>%
  dplyr::summarise(ndvi_mean = mean(ndvi, na.rm = TRUE),
                   ndvi_sd = sd(ndvi, na.rm = TRUE)) -> ndvi_stats

# Prep the data for modelling
rodents_table %>%
  dplyr::mutate(month = lubridate::month(newmoondate),
                year = lubridate::year(newmoondate)) %>%
  dplyr::filter(year > 1995) %>%
  # Scale continuous variables for massively improved efficiency of 
  # Stan sampling
  dplyr::mutate(ndvi = as.vector(scale(ndvi)),
                mintemp = as.vector(scale(mintemp)),
                maxtemp = as.vector(scale(maxtemp))) %>%
  dplyr::mutate(ndvi_ma12 = zoo::rollmean(ndvi, k = 12, align = 'right',
                                          na.pad = TRUE)) %>%
  # # Keep the first observation if multiple taken in the same month
  # dplyr::arrange(year, month) %>%
  # dplyr::group_by(month, year) %>%
  # dplyr::slice_head(n = 1) %>%
  tidyr::pivot_longer(cols = colnames(rodents_table)[4:24],
                      names_to = 'series', values_to = 'y') %>%
  dplyr::select(y, series, month, year, 
                newmoonnumber, mintemp:ndvi_ma12) %>%
  dplyr::mutate(time = newmoonnumber - (min(newmoonnumber) - 1))-> model_dat

# Many models will fail if the series of observations is nearly all zeroes. 
# Remove species with < 33 (out of 330) total unique observations (i.e. captures in at least
# 10% of unique trapping sessions)
# as forecasting these is not really useful anyway
model_dat %>%
  dplyr::group_by(series) %>%
  dplyr::summarise(total_obs = length(which(y >= 1))) %>%
  dplyr::filter(total_obs >= 33) %>%
  dplyr::pull(series) -> series_keep

model_dat %>%
  dplyr::filter(series %in% series_keep) %>%
  dplyr::filter(series != 'total') %>%
  dplyr::mutate(series = as.factor(series)) %>%
  dplyr::arrange(time, series) -> model_dat


# Feature engineering
#1. Distributed lag matrices for environmental covariates
# Function to set up a lag matrix for distributed lag nonlinear models
lagard <- function(x, n_lag = 6){
  n <- length(x)
  X <- matrix(NA, n, n_lag)
  for (i in 1:n_lag) X[i:n, i] <- x[i:n - i + 1]
  X
}

# Function to generate predictions for missing real-valued environmental variables
# using a GAM with seasonality and yearly components
approx_gam = function(df, family = gaussian()){
  require(mgcv)
  mod <- bam(y ~ 
               s(month, bs = 'cc', k = 10) +
               s(year, bs = 'bs', m = c(2,1,0),
                 k = 12), data = df, discrete = TRUE)
  preds <- predict(mod, newdata = df, type = 'response')
  
  # Replace any missing values with model-based predictions
  truth <- df$y
  truth[is.na(truth)] <- preds[is.na(truth)]
  truth
}

# Mintemp 6-month lag matrix
unique_times <- sort(unique(model_dat$time))
mintemp <- lagard(approx_gam(model_dat %>%
                               dplyr::select(mintemp, month, year, time) %>%
                               dplyr::arrange(time) %>%
                               dplyr::distinct() %>%
                               dplyr::mutate(y = mintemp)), 6)
mintemp_df <- data.frame(mintemp)
mintemp_df$time <- unique_times
model_dat %>%
  dplyr::select(time, year, month) %>%
  dplyr::left_join(mintemp_df) -> mintemp_df
dim(mintemp_df)[1] == NROW(model_dat)

# Maxtemp 6-month lag matrix
maxtemp <- lagard(approx_gam(model_dat %>%
                               dplyr::select(maxtemp, month, year, time) %>%
                               dplyr::arrange(time) %>%
                               dplyr::distinct() %>%
                               dplyr::mutate(y = maxtemp)), 6)
maxtemp_df <- data.frame(maxtemp)
maxtemp_df$time <- unique_times
model_dat %>%
  dplyr::select(time, year, month) %>%
  dplyr::left_join(maxtemp_df) -> maxtemp_df
dim(maxtemp_df)[1] == NROW(model_dat)

# The lag matrix
lag <- matrix(0:5, nrow(model_dat), 
              6, byrow = TRUE)
dim(lag)[1] == NROW(model_dat)

#2. Create remaining moving average / anomaly versions of environmental covariates
model_dat %>%
  dplyr::left_join(model_dat %>%
                     dplyr::select(time, mintemp, maxtemp, ndvi) %>%
                     dplyr::distinct() %>%
                     dplyr::mutate(mintemp_ma3 = zoo::rollmean(mintemp, k = 3, align = 'right',
                                                               na.pad = TRUE),
                                   maxtemp_ma3 = zoo::rollmean(maxtemp, k = 3, align = 'right',
                                                               na.pad = TRUE))) -> model_dat

# As we now have NAs for the first 11 rows of observations for each lag matrix, 
# as well as NAs for some rows of the moving average covariates,
# filter the data so that no NAs remain for covariates
model_dat %>%
  dplyr::filter(time > 11) %>%
  dplyr::mutate(time = time - 11) -> model_dat

# Impute ndvi_ma12
model_dat %>%
  dplyr::select(-ndvi_ma12) %>%
  dplyr::left_join(
    model_dat %>%
      dplyr::select(time) %>%
      dplyr::distinct() %>%
      dplyr::arrange(time) %>%
      dplyr::bind_cols(
        data.frame(ndvi_ma12 = approx_gam(model_dat %>%
                                            dplyr::select(ndvi_ma12, month, year, time) %>%
                                            dplyr::arrange(time) %>%
                                            dplyr::distinct() %>%
                                            dplyr::mutate(y = ndvi_ma12)))
      )
  ) -> model_dat


mintemp_df %>%
  dplyr::ungroup() %>%
  dplyr::filter(time > 11) %>%
  dplyr::select(-time, -year, -month) %>%
  as.matrix() -> mintemp
dim(mintemp)[1] == NROW(model_dat)

maxtemp_df %>%
  dplyr::ungroup() %>%
  dplyr::filter(time > 11) %>%
  dplyr::select(-time, -year, -month) %>%
  as.matrix() -> maxtemp
dim(maxtemp)[1] == NROW(model_dat)

lag <- tail(lag, NROW(model_dat))
dim(lag)[1] == NROW(model_dat)

# Now create weight matrices that can be used for setting up hierarchical 
# distributed lag terms
weights_dm <- weights_do <- 
  weights_pp <- weights_ol <-
  weights_ot <- weights_pf <- 
  weights_pb <- weights_pe <- weights_rm <-
  matrix(1, ncol = ncol(lag), nrow = nrow(lag))

weights_dm[!(model_dat$series == 'DM'), ] <- 0
weights_do[!(model_dat$series == 'DO'), ] <- 0
weights_ol[!(model_dat$series == 'OL'), ] <- 0
weights_ot[!(model_dat$series == 'OT'), ] <- 0
weights_pb[!(model_dat$series == 'PB'), ] <- 0
weights_pe[!(model_dat$series == 'PE'), ] <- 0
weights_pf[!(model_dat$series == 'PF'), ] <- 0
weights_pp[!(model_dat$series == 'PP'), ] <- 0
weights_rm[!(model_dat$series == 'RM'), ] <- 0


# Create a list to store the full dataset, including lag matrices and 
# moving averages for the environmental covariates 
data_all <- list(lag = lag, 
                 mintemp = mintemp,
                 mintemp_ma3 = model_dat$mintemp_ma3,
                 maxtemp = maxtemp,
                 maxtemp_ma3 = model_dat$maxtemp_ma3,
                 ndvi_ma12 = model_dat$ndvi_ma12,
                 weights_dm = weights_dm,
                 weights_do = weights_do,
                 weights_ol = weights_ol,
                 weights_ot = weights_ot,
                 weights_pb = weights_pb,
                 weights_pe = weights_pe,
                 weights_pf = weights_pf,
                 weights_pp = weights_pp,
                 weights_rm = weights_rm,
                 y = model_dat$y, 
                 month = model_dat$month, 
                 year = model_dat$year, 
                 series = model_dat$series,
                 time = model_dat$time)

# Split data into training and testing; stop training 
# at the end of 2018 so that 2019 can be evaluated. Conditions were 
# challenging in COVID and post-COVID, so evaluation of models may not be 
# as 'fair'
train_inds <- which(model_dat$year < 2019)

data_train <- lapply(seq_along(data_all), function(x){
  if(is.matrix(data_all[[x]])){
    data_all[[x]][train_inds,]
  } else {
    data_all[[x]][train_inds]
  }
})

test_inds <- which(model_dat$year %in% c(2019))

data_test <- lapply(seq_along(data_all), function(x){
  if(is.matrix(data_all[[x]])){
    data_all[[x]][test_inds,]
  } else {
    data_all[[x]][test_inds]
  }
})
names(data_train) <- names(data_test) <- names(data_all)

# Save all objects for forecast modelling
save(model_dat,
     data_all,
     data_train,
     data_test,
     file = 'data/rodents_data_tsobjects.rda')

# Plot some useful descriptors of the raw data
source('Functions/checking_functions.R')
plot_raw_series()
plot_raw_hists()
plot_raw_acfs()

totals <- data.frame(y = data_all$y,
                     time = data_all$time) %>%
  dplyr::group_by(time) %>%
  dplyr::summarise(total = sum(y))

data.frame(time = data_all$time, year = data_all$year) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(min_time = min(time)) %>%
  dplyr::filter(year > 1996) -> year_times

jpeg('Figures/total_series.jpg', width = 5, height = 3.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(0, 1.25, 0, 0))
truth <- totals$total
plot(1, type = "n", bty = 'L',
     xaxt = 'n',
     ylab = '',
     ylim = range(c(truth), na.rm = TRUE),
     xlim = c(0, length(c(truth))),
     xlab = '')
axis(1, at = year_times$min_time, labels = year_times$year, cex.axis = 1,
     tck= -0.025)
lines(x = 1:length(truth), y = truth, lwd = 2, col = "#8F2727")
title(main = 'Full community (9 included species)', cex.main = 1, line = 0.1,
      xpd = NA)
box(bty = 'l', lwd = 2)
title(ylab = 'Total captures', xpd = NA, line = 2.25)
dev.off()


# Training statistics
totals %>%
  dplyr::filter(time < 273) %>%
  summary()

data_all$year[totals %>%
                dplyr::filter(time < 273) %>%
                pull(total) %>%
                which.min()]

totals %>%
  dplyr::filter(time >= 273) %>%
  summary()

# Plot descriptors of covariates; first an STL decomposition of mintemp
mintemp_ts <- data_all$mintemp[which(data_all$series == 'DM'),1] 
mintemp_ts <- ts(mintemp_ts, start = c(1996, 12), frequency = 12)
mintemp_stl <- forecast::mstl(mintemp_ts)
mintemp_stl[,1] <- (mintemp_stl[,1] * mintemp_stats$mintemp_sd) + 
  mintemp_stats$mintemp_mean

jpeg('Figures/mintemp_stl.jpeg', width = 6.5, height = 6.5,
     res = 300, units = 'in')
par(mar=c(2.25, 4, 0.5, 1))
layout(matrix(1:3, nrow = 3))
plot(1, type = "n", bty = 'L',
     xaxt = 'n',
     ylab = 'Minimum temperature (°C)',
     ylim = range(mintemp_stl[,1], na.rm = TRUE),
     xlim = c(0, length(mintemp_stl[,1])),
     xlab = '')
lines(as.vector(mintemp_stl[,1]), lwd = 2)
box(bty = 'l', lwd = 2)
time_axis(labels = FALSE)

plot(1, type = "n", bty = 'L',
     xaxt = 'n',
     ylab = 'Trend component (scaled)',
     ylim = range(as.vector(scale(mintemp_stl[,2])), na.rm = TRUE),
     xlim = c(0, length(mintemp_stl[,1])),
     xlab = '')
lines(as.vector(scale(mintemp_stl[,2])), lwd = 2)
box(bty = 'l', lwd = 2)
time_axis(labels = FALSE)

plot(1, type = "n", bty = 'L',
     xaxt = 'n',
     ylab = 'Seasonality component (scaled)',
     ylim = range(as.vector(scale(mintemp_stl[,3])), na.rm = TRUE),
     xlim = c(0, length(mintemp_stl[,1])),
     xlab = '')
lines(as.vector(scale(mintemp_stl[,3])), lwd = 2)
box(bty = 'l', lwd = 2)
time_axis()
dev.off()

# Now a time series of NDVI
ndvi_ts <- model_dat$ndvi[which(model_dat$series == 'DM')]
ndvi_ts <- (ndvi_ts * ndvi_stats$ndvi_sd) + ndvi_stats$ndvi_mean
ndvi_ma12_ts <- model_dat$ndvi_ma12[which(model_dat$series == 'DM')]

jpeg('Figures/ndvi_ts.jpeg', width = 6.5, height = 5.5,
     res = 300, units = 'in')
par(mar=c(2.25, 4, 0.5, 1))
layout(matrix(1:2, nrow = 2))
plot(1, type = "n", bty = 'L',
     xaxt = 'n',
     ylab = 'NDVI (unitless)',
     ylim = range(ndvi_ts, na.rm = TRUE),
     xlim = c(0, length(ndvi_ts)),
     xlab = '')
lines(ndvi_ts, lwd = 2)
box(bty = 'l', lwd = 2)
time_axis(labels = FALSE)

plot(1, type = "n", bty = 'L',
     xaxt = 'n',
     ylab = 'NDVI moving average (scaled)',
     ylim = range(ndvi_ma12_ts, na.rm = TRUE),
     xlim = c(0, length(ndvi_ma12_ts)),
     xlab = '')
lines(ndvi_ma12_ts, lwd = 2)
box(bty = 'l', lwd = 2)
time_axis()
dev.off()

