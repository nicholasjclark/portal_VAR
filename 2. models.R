library(dplyr)
library(cmdstanr)
library(forecast)
library(mvgam)

# Load the pre-prepared modelling data
load('data/rodents_data_tsobjects.rda')

# Source the bespoke checking / graphical functions
source('Functions/checking_functions.R')

# View some of the raw time series
plot_mvgam_series(data = data_train, series = 'all')
plot_mvgam_series(data = data_train, newdata = data_test, series = 8)

#### Building up the GAM-VAR model ####
# Prior simulation for a baseline version of the GAM-VAR model; note that there will
# likely be many warnings as the predictions from spline prior models could be 
# extremely large
mod1_prior <- mvgam(
  formula = y ~ -1,
  # Species-level hierarchical slopes of NDVI
  # Note we use 'trend' here in place of 'series' to specify
  # the grouping factor. This is because mvgam will allow us to 
  # fit models where multiple time series share the same latent
  # process model, so it is possible that the number of groups 
  # in 'trend' will be fewer than the number of groups in 
  # 'series'. In this case each unique time series will have its
  # own dynamic process, but the effect of NDVI will be estimated
  # jointly
  trend_formula = ~ s(ndvi_ma12, trend, bs = 're') - 1,
  data = data_train,
  # Poisson observation model
  family = poisson(),
  # AR1 dynamics for the process model; note that this process
  # does not have to be centred about zero, so it can capture
  # variation in series-level intercepts without us needing to 
  # explicitly include varying intercept terms
  trend_model = 'AR1',
  burnin = 300,
  samples = 300,
  prior_simulation = TRUE
)

# Plot the prior for the random intercepts and random NDVI slopes, on the log scale
plot(mod1_prior, 're', trend_effects = TRUE)

# These are reasonable given the expected number of captures per session
# Satisfied with our prior model, let's now condition on the observed data
mod1 <- mvgam(
  formula = y ~ -1,
  trend_formula = ~ s(trend, bs = 're') + 
    s(ndvi_ma12, trend, bs = 're') - 1,
  data = data_train,
  newdata = data_test,
  family = poisson(),
  trend_model = AR(),
  algorithm = 'sampling',
  samples = 1600
)

# View the model summary to ensure the MCMC estimator is adequately exploring 
# the joint posterior without any obvious hindrances
summary(mod1)
dir.create('Outputs', showWarnings = FALSE, recursive = TRUE)
save(mod1, file = 'Outputs/mod1.rda')

# We expect some consistent seasonality in captures for the dominant species, so
# distributed lags of minimum temperature make sense. But we need to ensure the 
# computation is sound for this more complex model before adding any dynamic trend
# components
mod2 <- mvgam(
  formula = y ~ -1,
  trend_formula = ~ s(ndvi_ma12, trend, bs = 're') +
    # Hierarchical distributed lags of minimum temperature
    te(mintemp, lag, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')) - 1,
  data = data_train,
  newdata = data_test,
  family = poisson(),
  trend_model = AR(),
  samples = 1600,
  algorithm = 'sampling'
)

# Ensure the MCMC estimator is adequately exploring the joint 
# posterior without any obvious hindrances
summary(mod2)
save(mod2, file = 'Outputs/mod2.rda')

# We expect complex relationships among latent trends, but how 
# can these be captured in a principled way? A VAR(1) dynamic process
# First we need to use suitable containment priors for some key parameters
invgamma_opt_fun <- function(x, lowerbound, upperbound,
                             plower, pupper){
  x <- exp(x)
  val1 <- extraDistr::pinvgamma(lowerbound, x[1], x[2], log.p = TRUE)
  val2 <- extraDistr::pinvgamma(upperbound, x[1], x[2], lower.tail = FALSE, 
                                log.p = TRUE)
  c(val1 - log(plower), val2 - log(pupper))
}

# Now derive a containment for the variance of the NDVI random slopes, which 
# we do not expect to ever be smaller than 0.1 (on the log scale)
opt <- nleqslv::nleqslv(
  c(0, 0),
  invgamma_opt_fun,
  lowerbound = 0.1,
  upperbound = 1,
  plower = 0.01,
  pupper = 0.1,
  control = list(allowSingular = TRUE))

exp(opt$x)

# With prior distributions derived, we can construct the mvgam model with a latent
# VAR1 temporal process
priors <- get_mvgam_priors(
  formula = y ~ -1,
  trend_formula = ~ s(ndvi_ma12, trend, bs = 're') +
    te(mintemp, lag, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')),
  data = data_train,
  family = poisson(),
  trend_model = VAR(cor = TRUE)
)

# Update the prior for the species-level process model variances, with appropriate upper and 
# lower bounds. Here we use brms functionality to easily set priors
priors <- prior(beta(10, 10), class = sigma, lb = 0.2, ub = 1)

# Update the prior for the NDVI random slopes
priors <- c(priors,
            prior(inv_gamma(2.3693353, 0.7311319), class = sigma_raw_trend))
priors

# Ensure the priors are properly incorporated into the model's Stan code
code(mvgam(
  formula = y ~ -1,
  trend_formula = ~ s(ndvi_ma12, trend, bs = 're') +
    te(mintemp, lag, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')),
  data = data_train,
  newdata = data_test,
  family = poisson(),
  trend_model = VAR(cor = TRUE),
  priors = priors,
  run_model = FALSE
))

# Fit the model
modvar <- mvgam(
  formula = y ~ -1,
  trend_formula = ~ s(ndvi_ma12, trend, bs = 're') +
    te(mintemp, lag, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')),
  data = data_train,
  newdata = data_test,
  family = poisson(),
  trend_model = VAR(cor = TRUE),
  priors = priors,
  samples = 1600
)

# Calculate stability and forecast decomposition metrics; save
modvar_metrics <- calc_metrics(modvar, data_test, 'GAMVAR')
save(modvar, modvar_metrics, file = 'Outputs/modvar.rda')

# Now fit the same model to the full set of observed data
modvar_all <- update(modvar,
                     data = data_all,
                     samples = 1600,
                     algorithm = 'sampling')

# No metrics for this model as there is no testing data available
save(modvar_all, file = 'Outputs/modvar_all.rda')

#### Fit Benchmark models ####
# Same GAM linear predictor as GAM-VAR model, but with independent AR1 trends
# (in-text referred to as GAM-AR)
bench1 <- mvgam(
  formula = y ~ -1,
  trend_formula = ~
    s(ndvi_ma12, trend, bs = 're') +
    te(mintemp, lag, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')),
  data = data_train,
  newdata = data_test,
  family = poisson(),
  trend_model = AR(),
  samples = 1600,
  priors = priors
)
save(bench1, file = 'Outputs/bench1.rda')

bench1_all <- update(bench1,
                     data = data_all,
                     samples = 1600,
                     algorithm = 'sampling')
save(bench1_all, file = 'Outputs/bench1_all.rda')

# A slight decrease in complexity again, by removing the hierarchical components
# in the linear predictor and using no pooling to learn these
priors <- c(priors, prior(std_normal(), class = b))
bench1.2 <- mvgam(
  formula = y ~ -1,
  trend_formula = ~
    ndvi_ma12 * trend +
    te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_ol, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_pb, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_pe, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_pf, k = c(3, 4), bs = c('tp', 'cr')) +
    te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')),
  data = data_train,
  newdata = data_test,
  family = poisson(),
  trend_model = AR(),
  samples = 1600,
  priors = priors
)
save(bench1.2, file = 'Outputs/bench1.2.rda')

bench1.2_all <- update(bench1.2,
                       data = data_all,
                       samples = 1600,
                       algorithm = 'sampling')
save(bench1.2_all, file = 'Outputs/bench1.2_all.rda')

# Now the simplest benchmark, with no GAM linear predictor;
# This model uses independent AR1 trends with Poisson observations
# (in-text referred to as AR); note that the observation intercept is now meaningful in 
# this model, so we use a reasonable prior for it
priors <- c(prior(normal(1.5, 1), class = Intercept),
            prior(beta(10, 10), class = sigma, lb = 0.2, ub = 1))
bench2_skeleton <- mvgam(
  formula = y ~ 1,
  data = data_train,
  newdata = data_test,
  family = poisson(),
  trend_model = AR(),
  priors = priors,
  run_model = FALSE
)
code(bench2_skeleton)

# Looks good; run the model
bench2 <- mvgam(
  formula = y ~ 1,
  data = data_train,
  newdata = data_test,
  family = poisson(),
  trend_model = AR(),
  samples = 1600,
  priors = priors
)
save(bench2, file = 'Outputs/bench2.rda')

bench2_all <- update(bench2,
                     data = data_all,
                     samples = 1600,
                     algorithm = 'sampling')
save(bench2_all, file = 'Outputs/bench2_all.rda')

#### Compute exact leave-future-out cross-validation for each model
# across a set of evenly-spaced evaluation periods. When this is finished, we will
# have 6 total exact cross-validations for comparing models ####
evaluation_seq <- round(seq.int(75, 273, length.out = 6), 0)[1:5]
bench1_roll <- lapply(evaluation_seq, function(last_train){
  lfo_exact(object = bench1_all,
            data = data_all,
            last_train = last_train,
            fc_horizon = 12,
            stab_metrics = FALSE,
            samples = 1600)
})

bench1.2_roll <- lapply(evaluation_seq, function(last_train){
  lfo_exact(object = bench1.2_all,
            data = data_all,
            last_train = last_train,
            fc_horizon = 12,
            stab_metrics = FALSE,
            samples = 1600)
})

bench2_roll <- lapply(evaluation_seq, function(last_train){
  lfo_exact(object = bench2_all,
            data = data_all,
            last_train = last_train,
            fc_horizon = 12,
            stab_metrics = FALSE,
            samples = 1600)
})

# When iterating across GAMVAR models, calculate stability and 
# forecast decomposition metrics for each iteration
modvar_roll <- lapply(evaluation_seq, function(last_train){
  lfo_exact(object = modvar_all,
            data = data_all,
            last_train = last_train,
            fc_horizon = 12,
            stab_metrics = TRUE,
            samples = 1600)
})
save(bench1_roll, bench1.2_roll, bench2_roll, modvar_roll,
     evaluation_seq,
     file = 'Outputs/roll_evaluations.rda')

#### Notes ####
# Previous versions used Negative Binomial observation models without the state space
# representation. Suitable priors for the overdispersion parameters (phi) were derived using:
opt <- nleqslv::nleqslv(
  c(0, 0),
  invgamma_opt_fun,
  # Containment of r values at 2 to avoid overly flexible
  # dispersion models that will make it difficult to estimate a trend
  lowerbound = 2,
  upperbound = 600,
  # We only want 1% of prior probability to be below our lower threshold
  # of 2
  plower = 0.01,
  # Allow larger prior probability mass above the upper threshold
  pupper = 0.1,
  control = list(allowSingular = TRUE))
exp(opt$x)