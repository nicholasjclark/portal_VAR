library(dplyr)
library(cmdstanr)
library(forecast)
#remotes::install_github('nicholasjclark/mvgam')
library(mvgam)
setwd("C:/Users/Nick/Google Drive/Academic Work Folder/Ecological forecasting/mv_portalcasting/rodent_evaluation_ms")
source('Functions/checking_functions.R')

# Load the pre-prepared modelling data
load('data/rodents_data_tsobjects.rda')

# View some of the raw time series
plot_mvgam_series(data = data_train, series = 'all')
plot_mvgam_series(data = data_train, newdata = data_test, series = 8)


#### Building up the GAM-VAR model ####
# Prior simulation for a baseline version of the GAM-VAR model
mod1_prior <- mvgam(formula = y ~ 
                      # Series-level hierarchical intercepts
                      s(series, bs = 're') +
                      # Species-level hierarchical slopes of NDVI
                      s(ndvi_ma12, series, bs = 're') - 1,
                    data = data_train,
                    family = 'nb',
                    trend_model = 'None',
                    chains = 1,
                    prior_simulation = TRUE)

# Plot the prior for the random intercepts and slopes, on the log scale
plot(mod1_prior, 're')

# These are reasonable given the expected number of captures per session

# Satisfied with our prior model, let's now condition on the observed data
mod1 <- mvgam(formula = y ~ 
                s(series, bs = 're') +
                s(ndvi_ma12, series, bs = 're') - 1,
              data = data_train,
              newdata = data_test,
              family = 'nb',
              trend_model = 'None',
              use_stan = TRUE)

# View the model summary to ensure the MCMC estimator is adequately exploring 
# the joint posterior without any obvious hindrances
summary(mod1)

dir.create('Outputs', showWarnings = FALSE, recursive = TRUE)
save(mod1, file = 'Outputs/mod1.rda')

# We expect some consistent seasonality in captures for the dominant species, so
# distributed lags of minimum temperature make sense. But we need to ensure the 
# computation is sound for this more complex model before adding any dynamic trend
# components
mod2 <- mvgam(formula = y ~ 
                s(series, bs = 're') +
                s(ndvi_ma12, series, bs = 're') +
                te(mintemp, lag, k = c(3, 4), bs = c('tp', 'cr')) +
                te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
                te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
                te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
                te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')) - 1,
              data = data_train,
              newdata = data_test,
              family = 'nb',
              trend_model = 'None',
              use_stan = TRUE)

# Ensure the MCMC estimator is adequately exploring the joint 
# posterior without any obvious hindrances
summary(mod2)
save(mod2, file = 'Outputs/mod2.rda')


# Model 2 fits without issue and indicates there is evidence for seasonality in captures
# for most species. Now for model expansion: we expect complex relationships among latent trends, but how 
# can these be captured in a principled way? A latent VAR(1) dynamic process is one option

# First we need to use suitable containment priors for some key parameters; most notably
# the overdispersion parameters for the Negative Binomial sampling distribution. This code
# is modified from relevant code in the brms R package: https://github.com/paul-buerkner/brms
invgamma_opt_fun <- function(x, lowerbound, upperbound,
                             plower, pupper){
  x <- exp(x)
  val1 <- extraDistr::pinvgamma(lowerbound, x[1], x[2], log.p = TRUE)
  val2 <- extraDistr::pinvgamma(upperbound, x[1], x[2], lower.tail = FALSE, log.p = TRUE)
  c(val1 - log(plower), val2 - log(pupper))
}

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

# Now a containment for the variance of the NDVI random slopes, which 
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

# With prior distributions derived, we can construct the necessary data 
# objects for conditioning the GAM-VAR model
modvar_skeleton <- mvgam(formula = y ~ 
                           s(ndvi_ma12, series, bs = 're') +
                           te(mintemp, lag, k = c(3, 4), bs = c('tp', 'cr')) +
                           te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
                           te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
                           te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
                           te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')) - 1,
                         data = data_train,
                         newdata = data_test,
                         family = 'nb',
                         trend_model = 'AR1',
                         use_stan = TRUE,
                         run_model = FALSE)
model_data <- modvar_skeleton$model_data

# The model file has already been created externally (saved in the Stan directory). This 
# was necessary becase mvgam does not yet automatically create models with latent VAR processes.
# Compile the model into C++ code using Cmdstanr routines
cmd_mod <- cmdstan_model('Stan/var_gam.stan',
                         stanc_options = list('O1',
                                              'canonicalize=deprecations,braces,parentheses'))

# Condition the model; takes ~ 1.5 hours to sample on a 
# Intel(R) Core(TM) i5-8500 CPU @ 3.00GHz with 32GB RAM
fit <- cmd_mod$sample(data = model_data,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 500,
                      iter_sampling = 500,
                      refresh = 100,
                      init = modvar_skeleton$inits,
                      max_treedepth = 11)

# Read in the cmdstan files and convert to a stanfit object
out_gam_mod <- mvgam:::read_csv_as_stanfit(fit$output_files(),
                                           variables = c('b',
                                                         'mu_raw', 
                                                         'sigma_raw', 
                                                         'sigma',
                                                         'trend',
                                                         'ypred',
                                                         'mus',
                                                         'beta_var',
                                                         'rho', 
                                                         'r'))

# Convert the resulting object to mvgam class for easier 
# manipulations and summaries
out_gam_mod <- mvgam:::repair_stanfit(out_gam_mod)
modvar <- modvar_skeleton
modvar$max_treedepth <- 11
modvar$fit_engine <- 'stan'
class(modvar) <- 'mvgam'
modvar$model_output <- out_gam_mod
modvar$trend_model = 'AR1'
modvar$model_file <- cmd_mod$print()

modvar$resids <- mvgam:::get_mvgam_resids(object = list(
  model_output = out_gam_mod,
  fit_engine = 'stan',
  family = 'Negative Binomial',
  obs_data = modvar$obs_data,
  ytimes = modvar$ytimes,
  n_cores = 4))
names(modvar$resids) <- levels(data_train$series)
save(modvar, file = 'Outputs/modvar.rda')

# A final check of computational faithfulness and convergence
summary(modvar)

#### Fit Benchmark models ####
# Same GAM linear predictor as GAM-VAR model, but with independent AR1 trends
# (in-text referred to as GAM-AR)
priors <- get_mvgam_priors(formula = y ~ 
                             s(ndvi_ma12, series, bs = 're') +
                             te(mintemp, lag, k = c(3, 4), bs = c('tp', 'cr')) +
                             te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
                             te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
                             te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
                             te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')) - 1,
                           data = data_train,
                           family = 'nb',
                           trend_model = 'AR1',
                           use_stan = TRUE)
priors$prior[1] <- 'ar1 ~ normal(0.5, 0.25);'
priors$prior[2] <- 'sigma ~ beta(8, 12);'
priors$prior[3] <- 'r_inv ~ gamma(0.5408871, 6.8770244);'
bench1 <- mvgam(formula = y ~ 
           s(ndvi_ma12, series, bs = 're') +
           te(mintemp, lag, k = c(3, 4), bs = c('tp', 'cr')) +
           te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
           te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
           te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
           te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')) - 1,
         data = data_train,
         newdata = data_test,
         family = 'nb',
         priors = priors,
         trend_model = 'AR1',
         use_stan = TRUE)
save(bench1, file = 'Outputs/bench1.rda')

# Now a simpler benchmark with no GAM linear predictor;
# This model uses independent AR1 trends with Negative Binomial observations
# (in-text referred to as AR)
priors <- get_mvgam_priors(formula = y ~ 1,
                 data = data_train,
                 family = 'nb',
                 trend_model = 'AR1',
                 use_stan = TRUE)
priors$prior[1] <- 'ar1 ~ normal(0.5, 0.25);'
priors$prior[2] <- 'sigma ~ beta(8, 12);'
priors$prior[3] <- 'r_inv ~ gamma(0.5408871, 6.8770244);'
bench2 <- mvgam(formula = y ~ 1,
                data = data_train,
                newdata = data_test,
                family = 'nb',
                trend_model = 'AR1',
                use_stan = TRUE,
                priors = priors)
save(bench2, file = 'Outputs/bench2.rda')

