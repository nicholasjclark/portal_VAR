#### R script originally prepared by Nicholas J Clark (nicholas.j.clark1214@gmail.com) ####

#### Checking functions for bespoke model evaluation ####

#### Function to use exact leave-future-out cross-validation to evaluate models
# on a small subset of evaluation splits (5 splits) ####
lfo_exact = function(object, data, last_train, fc_horizon = 12,
                     stab_metrics = TRUE, model_name = 'GAMVAR', ...){
  
  # Ensure data and object will work with newer versions of mvgam
  data <- mvgam:::validate_series_time(data,
                                       trend_model = attr(object$model_data, 
                                                          'trend_model'))
  object$share_obs_params <- FALSE
  
  # Split data into training and testing at the appropriate time point
  data_split <- mvgam:::cv_split(data, last_train = last_train,
                                 fc_horizon = fc_horizon)
  
  # Update the model using the new data splits
  split_mod <- update(object,
                      data = data_split$data_train,
                      newdata = data_split$data_test,
                      lfo = !stab_metrics,
                      ...)
  
  # Extract forecast distributions
  fc <- forecast(split_mod)
  
  # Calculate forecast scores
  fc_score_energy <- data.frame(score = 
                                  score(fc, score = 'energy',
                                        log = TRUE)$all_series$score[1:fc_horizon])
  fc_score_energy$eval_timepoint <- last_train
  
  fc_score_var <- data.frame(score = score(fc, score = 'variogram', 
                                           log = TRUE)$all_series$score[1:fc_horizon])
  fc_score_var$eval_timepoint <- last_train
  
  cmbn_score <- fc_score_var
  cmbn_score$score <- log(fc_score_var$score * fc_score_energy$score)
  cmbn_score$score_type <- 'combination'
  
  if(stab_metrics){
    # Calculate stability and forecast composition metrics
    fc_metrics <- calc_metrics(object = split_mod,
                               data_test = data_split$data_test,
                               model_name = model_name)
  } else {
    fc_metrics <- NULL
  }

  # Return scores
  return(list(energy_score = fc_score_energy, 
              var_score = fc_score_var,
              cmbn_score = cmbn_score,
              fc_metrics = fc_metrics))
}

#### Functions to plot uncertainty contributions ####
plot_unc_props = function(unc_props){
  ggplot(unc_props, aes(x = horizon, 
                        y = `Proportion of variance`, 
                        fill = Process)) +
    geom_area() +
    scale_fill_viridis(discrete = TRUE) +
    facet_wrap(~ series) +
    labs(x = 'Forecast horizon (lunar months)',
         y = 'Proportion of forecast variance') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(legend.title = element_blank())
}

plot_sp_unc_time = function(fc_uncertainty_props,
                            series){
  s_name <- levels(data_all$series)[series]
  unc_props <- fc_uncertainty_props %>%
    dplyr::filter(series == s_name,
                  horizon < 13)
  train_labs <- paste0('T = ', 
                       unique(unc_props$end_train))
  names(train_labs) <- unique(unc_props$end_train)
  ggplot(unc_props, aes(x = horizon, 
                        y = `Proportion of variance`, 
                        fill = Process)) +
    geom_area() +
    scale_fill_viridis(discrete = TRUE) +
    facet_wrap(~ end_train,
               labeller = 
                 labeller(end_train = train_labs)) +
    labs(x = 'Forecast horizon (lunar months)',
         y = 'Proportion of forecast variance',
         title = species_names[series]) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(plot.title = element_text(face = "italic"),
          legend.title = element_blank(),
          strip.text = element_text(face = "italic"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.ticks = element_line()) +
    geom_hline(aes(yintercept=-Inf)) + 
    geom_vline(aes(xintercept=-Inf))
}

#### Wrapper function to calculate all stability / forecast composition metrics ####
calc_metrics = function(object, data_test, model_name){
  
  # Last training time
  end_train <- max(object$obs_data$time)
  
  # Calculate stability metrics and their uncertainties
  stab_metrics <- stability_metrics(object) %>%
    dplyr::mutate(end_train = end_train,
                  model_name = model_name)
  
  # Extract trend parameter estimates in the correct format
  trend_estimates <- trend_pars(object, data_test)
  
  # Calculate proportional contributions to forecast uncertainty
  
  # 1. No process error but including uncertainty in the NDVI
  # GAM params (ignoring other effects)
  coef_fix <- (1:length(coef(object$trend_mgcv_model)))[
    !grepl('ndvi', names(coef(object$trend_mgcv_model)))]
  gam_ndvi <- gamvar_unc(object,
                         data_test,
                         trend_pars = trend_estimates,
                         coefs_fix = coef_fix,
                         trend_fix = TRUE,
                         process_error = FALSE)
  
  # 2. Now including all params in the GAM
  gam_all <- gamvar_unc(object,
                        data_test,
                        trend_pars = trend_estimates,
                        trend_fix = TRUE,
                        process_error = FALSE)
  
  # 3. Now including uncertainty in AR / VAR params
  gamvar_noerror <- gamvar_unc(object,
                               data_test,
                               trend_pars = trend_estimates,
                               trend_fix = FALSE,
                               process_error = FALSE)
  
  # 4. Now including process error
  gamvar_all <- gamvar_unc(object,
                           data_test,
                           trend_pars = trend_estimates,
                           trend_fix = FALSE,
                           process_error = TRUE)
  
  # Calculate proportions for all series
  unc_props <- do.call(rbind, 
                       lapply(seq_len(nlevels(object$obs_data$series)),
                              function(i){
                                varmat <- rbind(gam_ndvi[,i],
                                                gam_all[,i],
                                                gamvar_noerror[,i],
                                                gamvar_all[,i])
                                
                                plot_propvar(varmat = varmat,
                                             series = i,
                                             varnames = c('NDVI effect',
                                                          'Mintemp smooth',
                                                          'Interactions',
                                                          'Process error'))                         
                                
                              }))
  
  unc_props %>%
    dplyr::group_by(series, horizon) %>%
    dplyr::mutate(n = sum(prop_var)) %>%
    dplyr::mutate(`Proportion of variance` = prop_var / n) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(series = levels(modvar$obs_data$series)[series]) %>%
    dplyr::mutate(Process = sub('_', ' ', process)) %>%
    dplyr::mutate(Process = factor(Process, levels = c(
      'Process error',
      'Interactions',
      'Mintemp smooth',
      'NDVI effect')))  %>%
    dplyr::mutate(end_train = end_train,
                  model_name = model_name) -> unc_props
  
  return(list(unc_props = unc_props,
              stab_metrics = stab_metrics))
}

#### Function to calculate, and optionally plot, 
# proportions of forecast variance ####
plot_propvar = function(varmat,
                        varnames,
                        series,
                        legend = TRUE,
                        draw = FALSE){
  
  # Tidy the varmat
  varmat_norm <- apply(varmat, 2, function(x) {x / max(x)})
  do.call(rbind, lapply(rev(1:NROW(varmat)), function(i){
    if(i > 1){
      data.frame(prop_var = pmax(0, varmat_norm[i,] - varmat_norm[i - 1]),
                 process = gsub(' ', '_', varnames)[i],
                 series = series,
                 horizon = 1:NCOL(varmat_norm))
    } else {
      data.frame(prop_var = varmat_norm[i,],
                 process = gsub(' ', '_', varnames)[i],
                 series = series,
                 horizon = 1:NCOL(varmat_norm))
    }
    
  })) -> out
  
  if(!draw){
    return(out)
  } else {
    pred_rel <- apply(varmat, 2, function(x) {x / max(x)})
    plot(1:NCOL(pred_rel),
         pred_rel[1,],
         ylim = c(0,1),
         type = 'n',
         bty = 'l',
         ylab = "Proportion of variance",
         xlab = "Forecast horizon",
         xaxs = "i",
         yaxs = "i")
    cols <- c('black',
              viridis::viridis(NROW(varmat))[2:NROW(varmat)])
    
    ciEnvelope(1:NCOL(pred_rel),
               rep(0, ncol(pred_rel)),
               pred_rel[1,],
               col = cols[1])
    ciEnvelope(1:NCOL(pred_rel),
               pred_rel[1,],
               pred_rel[2,],
               col = cols[2])
    if(NROW(varmat) > 2){
      ciEnvelope(1:NCOL(pred_rel),
                 pred_rel[2,],
                 pred_rel[3,],
                 col = cols[3])
    }
    if(NROW(varmat) > 3){
      ciEnvelope(1:NCOL(pred_rel),
                 pred_rel[3,],
                 pred_rel[4,],
                 col = cols[4])
    }
    if(NROW(varmat) > 4){
      ciEnvelope(1:NCOL(pred_rel),
                 pred_rel[4,],
                 pred_rel[5,],
                 col = cols[5])
    }
    if(legend){
      legend("topright",
             legend = rev(varnames),
             col = rev(cols),
             lty = 1, lwd = 6, bg = 'white',
             box.col = 'white')
    }
    
    box(bty = 'l', lwd = 2)
    return(out)
  }
}

#### Function to extract posterior estimates of trend parameters
# for use in mvgam internals ####
trend_pars = function(object, data_test){
  
  # Forecast horizon
  data_test <- mvgam:::validate_series_time(data_test, 
                                            name = 'data_test',
                                            trend_model = attr(object$model_data, 'trend_model'))
  data_train <- mvgam:::validate_series_time(object$obs_data, 
                                            name = 'obs_data',
                                            trend_model = attr(object$model_data, 'trend_model'))
  
  h <- max(data_test$index..time..index - max(data_train$index..time..index))
  ending_time <- max(data_train$index..time..index)
  
  # Trend parameters
  trend_pars <- mvgam:::extract_trend_pars(object,
                                           keep_all_estimates = FALSE,
                                           ending_time = ending_time)
  return(list(h = h,
              trend_pars = trend_pars))
}

#### Function to compute variance of forecasts under specific conditions
# i.e. using only single draw for the process or setting particular
# GAM coefs to zero ####
gamvar_unc = function(object,
                      data_test,
                      trend_pars,
                      driver_fix = TRUE,
                      coefs_fix = NULL,
                      trend_fix = TRUE,
                      process_error = TRUE){
  
  # Validate training and test data
  data_test <- mvgam:::validate_series_time(data_test, 
                                            name = 'data_test',
                                            trend_model = attr(object$model_data, 'trend_model'))
  data_train <- mvgam:::validate_series_time(object$obs_data, 
                                             name = 'obs_data',
                                             trend_model = attr(object$model_data, 'trend_model'))
  
  # Create trend Xp matrix
  make_trend_Xp = function(object, data_train, data_test){
    n_series <- nlevels(object$obs_data$series)
    Xp_trend <- mvgam:::trend_Xp_matrix(newdata = mvgam:::sort_data(data_test),
                                        trend_map = object$trend_map,
                                        mgcv_model = object$trend_mgcv_model)
    Xp_trend_last <- mvgam:::trend_Xp_matrix(newdata = data_train,
                                             trend_map = object$trend_map,
                                             mgcv_model = object$trend_mgcv_model)
    
    # Ensure the last three observed values are used, in case the obs_data
    # was not supplied in order
    data.frame(time = data_train$index..time..index,
               series = object$obs_data$series,
               row_id = 1:length(data_train$index..time..index)) %>%
      dplyr::arrange(time, series) %>%
      dplyr::pull(row_id) -> sorted_inds
    
    linpred_order <- vector(length = 3 * n_series)
    last_rows <- tail(sort(sorted_inds), 3 * n_series)
    for(i in seq_along(last_rows)){
      linpred_order[i] <- which(sorted_inds == last_rows[i])
    }
    
    # Deal with any offsets
    if(!all(attr(Xp_trend_last, 'model.offset') == 0)){
      offset_vec <- attr(Xp_trend_last, 'model.offset')
      offset_last <- offset_vec[linpred_order]
      offset_last[is.na(offset_last)] <- 0
      full_offset <- c(offset_last, attr(Xp_trend, 'model.offset'))
    } else {
      full_offset <- 0
    }
    
    # Bind the last 3 linpred rows with the forecast linpred rows
    Xp_trend <- rbind(Xp_trend_last[linpred_order, , drop = FALSE],
                      Xp_trend)
    attr(Xp_trend, 'model.offset') <- full_offset
    
    return(Xp_trend)
  }
  
  # Posterior GAM betas
  betas_trend <- mvgam:::mcmc_chains(object$model_output, 'b_trend')
  
  if(!is.null(coefs_fix)){
    # Fix specified coefs to zero so particular GAM terms can
    # be dropped from the predictions
    betas_trend[,coefs_fix] <- 0
  }
  
  # Propagate the trend for 250 draws of betas, using the fixed
  # seed to ensure the VAR process is the same in each draw
  gam_draws <- lapply(seq_len(min(250, NROW(betas_trend))), 
                                 function(i){
    # Sample the driver estimates if specified
    if(driver_fix){
      Xp_trend <- make_trend_Xp(object, data_train, data_test)
    } else {
      # This assumes data_test is a list of possible test data sets
      ind <- sample(1:length(data_test), 1)
      Xp_trend <- make_trend_Xp(object, data_train, data_test[[ind]])
    }
    propagate_trends(h = trend_pars$h,
                     trend_pars = trend_pars$trend_pars,
                     data_test = data_test,
                     Xp_trend = Xp_trend,
                     betas_trend = betas_trend[i,],
                     seed = trend_fix,
                     samp_index = ifelse(trend_fix, 1, i),
                     process_error = process_error)
  })
  
  # Calculate variances of predictions across horizons for each series
  n_series <- nlevels(object$obs_data$series)
  gam_vars <- matrix(NA, nrow = NROW(gam_draws[[1]]),
                     ncol = n_series)
  for(i in 1:n_series){
    gam_vars[,i] <- apply(do.call(cbind,
                                  lapply(gam_draws, `[`,,i)), 1, var)
  }
  
  return(gam_vars)
}


#### Function to propagate dynamic processes while allowing for certain
# components to be ignored or fixed ####
propagate_trends = function(h,
                            trend_pars,
                            data_test,
                            Xp_trend = NULL,
                            betas_trend = NULL,
                            samp_index = 1,
                            seed = FALSE,
                            process_error = TRUE){
  
  if(!'last_lvs' %in% names(trend_pars)){
    trend_pars$last_lvs <- trend_pars$last_trends
  }
  
  # One realisation of trend parameters
  trend_pars <- mvgam:::extract_general_trend_pars(trend_pars = trend_pars,
                                                   samp_index = samp_index)
  
  # Reconstruct the A and Sigma matrices
  if('A' %in% names(trend_pars)){
    Amat <- matrix(trend_pars$A, nrow = length(trend_pars$last_lvs),
                   ncol = length(trend_pars$last_lvs),
                   byrow = TRUE)
    ar1 <- rlang::missing_arg()
  } else if('ar1' %in% names(trend_pars)){
    ar1 <- trend_pars$ar1
    Amat <- rlang::missing_arg()
  } else {
    ar1 <- rep(1, length(trend_pars$last_lvs))
    Amat <- rlang::missing_arg()
  }
  
  if('ar2' %in% names(trend_pars)){
    ar2 <- trend_pars$ar2
  } else {
    ar2 <- rep(0, length(trend_pars$last_lvs))
  }
  
  if('ar3' %in% names(trend_pars)){
    ar3 <- trend_pars$ar3
  } else {
    ar3 <- rep(0, length(trend_pars$last_lvs))
  }
  
  if('Sigma' %in% names(trend_pars)){
    Sigmamat <- matrix(trend_pars$Sigma, nrow = length(trend_pars$last_lvs),
                       ncol = length(trend_pars$last_lvs),
                       byrow = TRUE)
  } else if('sigma' %in% names(trend_pars)){
    Sigmamat <- matrix(0, nrow = length(trend_pars$last_lvs),
                       ncol = length(trend_pars$last_lvs),
                       byrow = TRUE)
    diag(Sigmamat) <- trend_pars$sigma
  } else {
    Sigmamat <- matrix(0, nrow = length(trend_pars$last_lvs),
                       ncol = length(trend_pars$last_lvs),
                       byrow = TRUE)
    diag(Sigmamat) <- 1 / trend_pars$tau
  }
  
  if(!process_error){
    Sigmamat <- matrix(0, ncol = ncol(Sigmamat),
                       nrow = nrow(Sigmamat))
    diag(Sigmamat) <- .Machine$double.eps
  }
  
  
  # Reconstruct the last trend matrix
  last_trendmat <- do.call(cbind,(lapply(trend_pars$last_lvs,
                                         function(x) tail(x, 3))))
  
  # If this is a moving average model, reconstruct theta matrix and
  # last error matrix
  if('theta' %in% names(trend_pars)){
    thetamat <- matrix(trend_pars$theta,
                       nrow = length(trend_pars$last_lvs),
                       ncol = length(trend_pars$last_lvs),
                       byrow = TRUE)
    errormat <- rbind(rep(0, length(trend_pars$last_lvs)),
                      rep(0, length(trend_pars$last_lvs)),
                      tail(trend_pars$error, length(trend_pars$last_lvs)))
    
  } else {
    thetamat <- rlang::missing_arg()
    errormat <- rlang::missing_arg()
  }
  
  # Prep VARMA parameters
  varma_params <- mvgam:::prep_varma_params(A = Amat,
                                            ar1 = ar1,
                                            ar2 = ar2,
                                            ar3 = ar3,
                                            Sigma = Sigmamat,
                                            last_trends = last_trendmat,
                                            last_errors = errormat,
                                            theta = thetamat,
                                            Xp_trend = Xp_trend,
                                            betas_trend = betas_trend,
                                            h = h)
  
  # Simulate one realisation of VAR process for h timesteps ahead
  if(seed){
    set.seed(1)
  } else {
    set.seed(NULL)
  }
  trend_realisation <- mvgam:::sim_varma(A = varma_params$A,
                                         A2 = varma_params$A2,
                                         A3 = varma_params$A3,
                                         drift = varma_params$drift,
                                         theta = varma_params$theta,
                                         Sigma = varma_params$Sigma,
                                         last_trends = varma_params$last_trends,
                                         last_errors = varma_params$last_errors,
                                         Xp_trend = varma_params$Xp_trend,
                                         betas_trend = varma_params$betas_trend,
                                         h = varma_params$h)
  
  return(trend_realisation)
}


#### Function to compute measures of stability
# These measures of stability assess how systems respond to
# environmental fluctuations, not how variable the systems are in general ####
stability_metrics = function(object){
  
  # Take posterior draws of the interaction matrix
  B_post <- as.matrix(object, variable = 'A', regex = TRUE)
  
  # Take posterior draws of Sigma
  Sigma_post <- as.matrix(object, variable = 'Sigma', regex = TRUE)
  
  metrics <- do.call(rbind, lapply(
    seq_len(min(1000, NROW(B_post))),
    function(i){
      
      B <- matrix(B_post[i,],
                  nrow = nlevels(object$obs_data$series),
                  ncol = nlevels(object$obs_data$series))
      p <- dim(B)[1]
      
      # If we want to get the variance of the stationary distribution (Sigma_inf)
      Sigma <- matrix(Sigma_post[i,],
                      nrow = nlevels(object$obs_data$series),
                      ncol = nlevels(object$obs_data$series))
      vecS_inf <- solve(diag(p * p) - kronecker(B, B)) %*% as.vector(Sigma)
      Sigma_inf <- matrix(vecS_inf, nrow = p)

      # The difference in volume between Sigma_inf and Sigma is:
      # det(Sigma_inf - Sigma) = det(Sigma_inf) * det(B) ^ 2
      # according to Ives et al 2003 (eqn 24)
      
      # We can take partial derivatives to determine which elements of 
      # Sigma_inf contribute
      # most to rates of change in the proportion of Sigma_inf that is due to 
      # unmodelled environmental variation (process error)
      int_env <- det(Sigma_inf) * t(solve(Sigma_inf))
      
      # Proportion of interspecific covariance to
      # to overall environmental variation contribution (i.e. how important are
      # correlated errors for controlling the shape of the stationary forecast
      # distribution?)
      dat <- data.frame(intersp_env_cont = mean(abs(int_env[lower.tri(int_env)])) /
        (mean(abs(diag(int_env))) + mean(abs(int_env[lower.tri(int_env)]))))
      
      # Proportion of volume of Sigma_inf attributable to species interactions,
      # measuring the degree to which species interactions increase
      # the variance of the stationary distribution (Sigma_inf) relative
      # to the variance of the process error (Sigma)
      # lower values = more stability
      dat$sp_prop = abs(det(B)) ^ 2
      
      # Ives et al 2003 suggest to scale this by the number of series for more direct
      # comparisons among different studies
      dat$sp_prop_adj <- abs(det(B)) ^ (2 / p)
      
      # Sensitivity of the species interaction proportion to particular
      # interactions is also calculated using partial derivatives 
      # (note the use of 2 here because we squared det(B) in the above eqn)
      int_sens <- 2 * det(B) * t(solve(B))
      
      # Proportion of interspecific contributions to
      # to overall interaction contribution
      dat$intersp_interact_cont <- mean(abs(int_sens[lower.tri(int_sens)])) /
        (mean(abs(diag(int_sens))) + mean(abs(int_sens[lower.tri(int_sens)])))
      
      # Reactivity, measuring the degree to which the system moves
      # away from a stable equilibrium following a perturbation
      # values > 0 suggest the system is reactive, whereby a
      # perturbation of the system in one period can be amplified in the next period
      # Following Neubert et al 2009 Ecology (Detecting reactivity)
      dat$reactivity <- log(max(svd(B)$d))
      
      # Return rate of transition distribution to the stationary distribution
      # Asymptotic return rate of the mean
      # lower values = more stability
      dat$mean_returnrate <- max(abs(eigen(B)$values))
      
      # Asymptotic return rate of the variance
      # lower values = more stability
      dat$var_returnrate <- max(abs(eigen(B %x% B)$values))
      dat
    }))
  return(metrics)
}

#### Generic function to plot a line histogram ####
# Most of this was adapted from Michael Betancourt's case study code:
# https://betanalpha.github.io/assets/case_studies/taylor_models.html
plot_line_hist <- function(s, min_val, max_val, 
                           delta = 0.02, 
                           line_col = 'black',
                           poly_col = 'black') {
  bins <- seq(min_val, max_val, delta)
  B <- length(bins) - 1
  idx <- rep(1:B, each=2)
  x <- sapply(1:length(idx),
              function(b) if(b %% 2 == 1) bins[idx[b]] else bins[idx[b] + 1])
  x <- c(min_val - delta, x, max_val + delta)
  
  counts <- hist(s, breaks=bins, plot=FALSE)$density
  y <- counts[idx]
  y <- c(0, y, 0)
  
  lines(x, y, lwd = 2, col = 'white')
  lines(x, y, lwd = 1.5, col = line_col)
  polygon(x, y, col = poly_col, border = NA)
}

#### Species names for easier plotting ####
species_names <- c('Dipodomys merriami',
                   'Dipodomys ordii',
                   'Onychomys leucogaster',
                   'Onychomys torridus',
                   'Chaetodipus baileyi',
                   'Peromyscus eremicus',
                   'Perognathus flavus',
                   'Chaetodipus penicillatus',
                   'Reithrodontomys megalotis')

abbrev_names <- c('D. merriami',
                   'D. ordii',
                   'O. leucogaster',
                   'O. torridus',
                   'C. baileyi',
                   'P. eremicus',
                   'P. flavus',
                   'C. penicillatus',
                   'R. megalotis')

#### For formatting time series x axes ####
data.frame(time = data_all$time, year = data_all$year) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(min_time = min(time)) %>%
  dplyr::filter(year > 1996) -> year_times

time_axis = function(labels = TRUE){
  if(labels){
    axis(1, at = year_times$min_time, labels = year_times$year, cex.axis = 1,
         tck= -0.05)
  } else {
    axis(1, at = year_times$min_time, labels = NA, cex.axis = 1,
         tck= -0.05)
  }

}

#### Plot raw data descriptions ####
plot_raw_series = function(filepath = 'Figures/raw_series.jpg'){
  raw_series <- lapply(seq_along(levels(data_all$series)), function(x){
    data_all$y[which(data_all$series == levels(data_all$series)[x])]
  })
  
  jpeg(filepath, width = 6.25, height = 4.25,
       res = 300, units = 'in')
  par(mar=c(2, 2, 1, 1),
      oma = c(0, 1.25, 0, 0))
  layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
  for(x in 1:9){
    truth <- raw_series[[x]]
    plot(1, type = "n", bty = 'L',
         xaxt = 'n',
         ylab = 'No. captures',
         ylim = range(c(truth), na.rm = TRUE),
         xlim = c(0, length(c(truth))),
         xlab = '')
    if(x > 6){
      time_axis()
    } else {
      time_axis(labels = FALSE)
    }
    lines(x = 1:length(truth), y = truth, lwd = 1.25, col = "#8F2727")
    title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.1,
          xpd = NA)
    box(bty = 'l', lwd = 2)
    
    if(x == 4){
      title(ylab = 'Number of captures', xpd = NA, line = 2.25)
    }
    
  }
  dev.off()
}

plot_raw_hists = function(filepath = 'Figures/raw_hists.jpg'){
  raw_series <- lapply(seq_along(levels(data_all$series)), function(x){
    data_all$y[which(data_all$series == levels(data_all$series)[x])]
  })
  
  jpeg(filepath, width = 6.25, height = 4.25,
       res = 300, units = 'in')
  par(mar=c(2, 2, 1, 1),
      oma = c(1.25, 1.25, 0, 0))
  layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
  for(x in 1:9){
    truth <- raw_series[[x]]
    hist(truth, border = 'white', col = "#8F2727", breaks = 30,
         main = '', freq = TRUE)
    title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.1,
          xpd = NA)
    box(bty = 'l', lwd = 2)
    if(x == 4){
      title(ylab = 'Frequency', xpd = NA, line = 2.25)
    }
    if(x == 8){
      title(xlab = 'Number of captures', xpd = NA, line = 2.25)
    }
  }
  dev.off()
}

plot_raw_acfs = function(filepath = 'Figures/raw_acfs.jpg'){
  raw_series <- lapply(seq_along(levels(data_all$series)), function(x){
    data_all$y[which(data_all$series == levels(data_all$series)[x])]
  })
  
  jpeg(filepath, width = 6.25, height = 4.25,
       res = 300, units = 'in')
  par(mar=c(2, 2, 1, 1),
      oma = c(1.25, 1.25, 0, 0))
  layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
  for(x in 1:9){
    truth <- raw_series[[x]]
    acf(c(truth),
        na.action = na.pass, bty = 'L',
        lwd = 2.5, ci.col = 'black', col = "#8F2727",
        main = '', ylab = 'Autocorrelation')
    acf1 <- acf(c(truth), plot = F,
                na.action = na.pass)
    clim <- qnorm((1 + .95)/2)/sqrt(acf1$n.used)
    abline(h = clim,  col = '#FFFFFF', lwd = 1.85)
    abline(h = clim,  col = 'black', lwd = 1.5, lty = 'dashed')
    abline(h = -clim,  col = '#FFFFFF', lwd = 1.85)
    abline(h = -clim,  col = 'black', lwd = 1.5, lty = 'dashed')
    box(bty = 'L', lwd = 2)
    
    title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.1,
          xpd = NA)
    if(x == 4){
      title(ylab = 'Autocorrelation function', xpd = NA, line = 2.25)
    }
    if(x == 8){
      title(xlab = 'Temporal lag (lunar months)', xpd = NA, line = 2.25)
    }
  }
  dev.off()
}


#### Get constrained forecast distribution ####
get_fc_constrained = function(object, 
                             bound = 196,
                             newdata){
  
  # Extract forecasts for all species
  data_train <- object$obs_data
  ends <- seq(0, dim(MCMCvis::MCMCchains(object$model_output, 'ypred'))[2],
              length.out = NCOL(object$ytimes) + 1)
  starts <- ends + 1
  starts <- c(1, starts[-c(1, (NCOL(object$ytimes)+1))])
  ends <- ends[-1]
  
  all_fcs <- lapply(seq_along(object$resids), function(series){
    if(object$fit_engine == 'stan'){
      preds <- MCMCvis::MCMCchains(object$model_output, 'ypred')[,seq(series,
                                                                      dim(MCMCvis::MCMCchains(object$model_output,
                                                                                              'trend'))[2],
                                                                      by = NCOL(object$ytimes))]
    } else {
      preds <- MCMCvis::MCMCchains(object$model_output, 'ypred')[,starts[series]:ends[series]]
    }
    preds
  })
  
  # Adjust forecasts based on the supplied upper bound, which by default
  # is set as the usual number of traps available across all control plots
  # (196)
  totals <- Reduce('+', all_fcs)
  totals[totals < bound] <- 1
  offsets <- matrix(bound, nrow = NROW(totals), ncol = NCOL(totals))
  offsets[totals < bound] <- 1
  
  preds_cons <- lapply(seq_along(all_fcs), function(x){
    # If the total predicted trappings is above the upper bound,
    # re-scale species-level forecasts based on their empirical proportions
    # of the total so that the total now matches the upper bound
    floor((all_fcs[[x]] / totals) * offsets)
  })
  names(preds_cons) <- names(object$resids)
  
  return(preds_cons)
}

#### Variogram score probabilistic forecast evaluation ####
evaluate_variogram = function(object, newdata, bound = 196,
                              weights = NULL){
  
  truth_data <- data.frame(y = newdata$y,
                          time = newdata$time,
                          series = newdata$series)
  
  # Get truths (out of sample) into correct format
  truths <- do.call(rbind, lapply(seq_along(object$resids), function(series){
    s_name <- levels(object$obs_data$series)[series]
    truth_data %>%
      dplyr::filter(series == s_name) %>%
      dplyr::select(time, y) %>%
      dplyr::distinct() %>%
      dplyr::arrange(time) %>%
      dplyr::pull(y)
  }))
  
  # Now produce a list of forecasts, one for each time point the fc_horizon
  # Note, this assumes the data in newdata have already been fed to the model
  # as missing observations!!!!
  fc_length <- length(unique(newdata$time))
  fcs <- get_fc_constrained(object, bound = bound)
  fcs <- lapply(seq_along(fcs), function(series){
    preds <- fcs[[series]]
    preds[, tail(1:NCOL(preds), fc_length)]
  })

  fcs_per_horizon <- lapply(seq_len(fc_length), function(horizon){
    do.call(rbind, lapply(seq_along(fcs), function(fc){
      fcs[[fc]][,horizon]
    }))
  })
  
  # Compute the variogram score, using the median pairwise difference
  # from the forecast distribution (scoringRules::vs_sample uses the 
  # mean, which is not appropriate for skewed distributions)
  variogram_score = function(obs, fc_sample){
    
    # log the truth and fc so that we aren't attributing huge amount of weight
    #   # to the more abundant species
      obs <- log(obs + 0.001)
      fc_sample <- log(fc_sample + 0.001)

    out <- matrix(NA, length(obs), length(obs))
    for(i in 1:length(obs)){
      for(j in 1:length(obs)){
        if(i == j){
          out[i,j] <- 0
        } else {
          v_fc <- quantile(abs(fc_sample[i,] -  fc_sample[j,]) ^ 0.5, 0.5)
          v_dat <- abs(obs[i] - obs[j]) ^ 0.5
          out[i,j] <- 2 * ((v_dat - v_fc) ^ 2)
        }
      }
    }
    # Divide by two as we have (inefficiently) computed each pairwise
    # comparison twice
    sum(out / 2)
  }
  
  # Calculate variogram score on all observations in fc_horizon
  vg_scores <- unlist(lapply(seq_len(fc_length), function(horizon){
    variogram_score(obs = truths[,horizon], 
                    fc_sample = fcs_per_horizon[[horizon]])   
  }))
  
  # Calculate scaled DRPS for weighting the variogram score by sharpness
  drps_scores <- rowSums(do.call(cbind, lapply(seq_along(fcs), function(series){
    mvgam:::drps_mcmc_object(fc = log(t(fcs[[series]]) + 0.001), 
                             truth = log(truths[series,] + 0.001))[,1]
  }))) 
  vg_scores * drps_scores
}

#### Function to fit a LOESS trendline to exact leave-future-out forecast scores ####
loess_scores = function(scores){
  predict(loess(y ~ time, span = 0.95, data = data.frame(y = log(scores),
                                                         time = seq_len(length(scores)))),
          newdata = data.frame(time = seq_len(length(scores))))
}

#### Function for rolling evaluations with weighted variogram / DRPS scores ####
roll_evaluation = function(object, 
                           evaluation_seq, 
                           n_samples, 
                           n_cores){
  # logged variogram evaluations
  var_roll <- roll_eval_mvgam(object, 
                              evaluation_seq = evaluation_seq,
                              n_samples = n_samples, 
                              n_cores = n_cores,
                              fc_horizon = 12,
                              score = 'variogram',
                              log = TRUE)
  
  # logged DRPS evaluations (per series)
  drps_roll <- roll_eval_mvgam(object, 
                               evaluation_seq = evaluation_seq,
                               n_samples = n_samples, 
                               n_cores = n_cores,
                               fc_horizon = 12,
                               score = 'drps',
                               log = TRUE)
  
  # Calculate summed DRPS scores per horizon
  drps_sum_scores <- data.frame(score = rowSums(do.call(cbind, lapply(seq_along(drps_roll$series_evals), 
                                                 function(series){
                                                   drps_roll$series_evals[[series]]$all_scores[,1]
                                                 }))),
                                horizon = as.factor(var_roll$all_scores$eval_horizon))
  # Variogram scores per horizon
  var_scores <- data.frame(score = var_roll$all_scores$score,
                           horizon = as.factor(var_roll$all_scores$eval_horizon))
  # Equally-weighted ensemble scores per horizon
  combn_scores <- data.frame(score = var_scores$score * drps_sum_scores$score,
                             horizon = as.factor(var_roll$all_scores$eval_horizon))
  return(list(drps_sum_scores = drps_sum_scores,
              var_scores = var_scores,
              combn_scores = combn_scores,
              coverage_90 = drps_roll$interval_coverage))
}

#### Plot mintemp trends against residuals ####
plot_mintemp_resids = function(object, series){
  ts_dat <- data_train$mintemp_ma3[which(data_train$series == 
                                              levels(model_dat$series)[series])]
  layout(matrix(c(1,1,2,3), nrow = 2, byrow = T))
  lims <- sort(c(max(object$resids[[series]][1:500,], na.rm = TRUE), 
            max(object$resids[[series]][1:500,], na.rm = TRUE) * -1 * 
              sign(max(object$resids[[series]][1:500,], na.rm = TRUE))))
  plot(as.vector(scale(ts_dat)),
       type = 'l', lwd = 2, col = 'darkred', bty = 'l',
       ylab = 'Value',
       main = levels(data_train$series)[series],
       ylim = lims,
       xaxt = 'n',
       xlab = '')
  time_axis()
  lincoefs <- vector(length = 500)
  for(i in 1:500){
    lines(object$resids[[series]][i,],
          col = '#D3D3D330',
          lwd = 0.75)
    lincoefs[i] <- coef(lm(as.vector(scale(ts_dat)) ~ object$resids[[series]][i,]))[2]
  }
  lines(as.vector(scale(ts_dat)),
        col = 'white', lwd = 3.5)
  lines(as.vector(scale(ts_dat)),
        col = 'darkred', lwd = 3)
  box(bty = 'l', lwd = 2)
  text(x = -5,
       y = lims[2],
       adj = 0,
       labels = 'Mintemp moving average',
       col = 'darkred',
       xpd = TRUE)
  text(x = -5,
       y = lims[2] - 1.25,
       adj = 0,
       labels = 'Residual draws',
       col = 'grey50',
       xpd = TRUE)
  lims <- c(max(abs(lincoefs)), 
            max(abs(lincoefs)) * -1 * sign(max(abs(lincoefs))))
  hist(lincoefs, breaks = seq(min(lims),
                              max(lims),
                              length.out = 20), 
       col = 'darkred',
       border = 'white', xlab = bquote(beta[mintemp]), 
       freq = FALSE,
       main = '',
       lwd = 2,
       ylab = 'Density')
  abline(v = 0, lwd = 3.5, col = 'white')
  abline(v = 0, lwd = 3, col = 'black')
  
  acfs <- vector(length = 500)
  for(i in 1:500){
    acfs[i] <- acf(object$resids[[series]][i,], lag.max = 12, plot = FALSE,
                   na.action = na.pass)$acf[[13]]
  }
  lims <- c(max(abs(acfs)), max(abs(acfs)) * -1 * sign(max(abs(acfs))))
  hist(acfs, breaks = seq(min(lims),
                          max(lims), 
                          length.out = 20), 
       col = 'darkred',
       border = 'white', xlab = 'Lag-12 ACF', 
       freq = FALSE,
       main = '',
       lwd = 2,
       ylab = '')
  abline(v = 0, lwd = 3.5, col = 'white')
  abline(v = 0, lwd = 3, col = 'black')
  layout(1)
}

#### Plot NDVI retrodictive residuals ####
ndvi_retro_residual <- function(series, object, n_bins = 20, 
                                xlabel, show_xlabs){
  
  # Plotting colours
  c_light <- c("#DCBCBC")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  
  # Pull out observed NDVI_ma12 values for this species
  
  obs_covar <- object$obs_data$ndvi_ma12[which(object$obs_data$series == 
                                                 levels(object$obs_data$series)[series])]
  
  # Pull out trend for this species
  resids <- object$resids[[series]]
  trend <- MCMCvis::MCMCchains(object$model_output, 'trend')[,seq(series,
                                                                  dim(MCMCvis::MCMCchains(object$model_output, 'trend'))[2],
                                                                  by = NCOL(object$ytimes))]
  trend <- trend[,1:length(obs_covar)]
  
  # Set up bins for trend component
  breaks <- seq(min(obs_covar), max(obs_covar), 
                (max(obs_covar) - min(obs_covar)) / n_bins)
  idx <- rep(1:n_bins, each=2)
  xs <- sapply(1:length(idx), function(b)
    if(b %% 2 == 1) breaks[idx[b]] else breaks[idx[b] + 1])
  
  S <- length(obs_covar)
  
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- matrix(NA, nrow=9, ncol=n_bins)
  
  for (b in 1:n_bins) {
    bin_idx <- which(breaks[b] <= obs_covar & obs_covar < breaks[b + 1])
    
    if (length(bin_idx)) {
      cred[,b] <- quantile(trend[,bin_idx] + resids[,bin_idx], 
                           probs = probs, na.rm = TRUE)
    }
    if ( is.na(cred[1, b]) ) cred[,b] = rep(0, length(probs))
  }
  
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9, n]))
  
  ymin <- min(pad_cred[1,], na.rm=T)
  ymax <- max(pad_cred[9,], na.rm=T)
  
  # Generate the plot
  if(xlabel){
    xlab <- 'NDVI moving average'
  } else {
    xlab <- ''
  }
  plot(1, type="n", main="",
       xlim=c(min(obs_covar), max(obs_covar)), 
       xlab=xlab,
       ylim=c(ymin, ymax), ylab="DS residual",
       xaxt = 'n',
       bty = 'l')
  if(show_xlabs){
    axis(1)
  } else {
    axis(1, labels = NA)
  }
  
  polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  
  for (b in 1:n_bins) {
    if (length(which(breaks[b] <= obs_covar & obs_covar < breaks[b + 1]))) {
      lines(xs[(2 * b - 1):(2 * b)], pad_cred[5,(2 * b - 1):(2 * b)], col=c_dark, lwd=2)
    }
  }
  box(bty = 'l', lwd = 2)
}


#### Plot NDVI trends against residuals ####
plot_ndvi_resids = function(object, series){
  ts_dat <- data_train$ndvi_ma12[which(data_train$series == 
                                            levels(model_dat$series)[series])]
  layout(matrix(c(1,1,2,3), nrow = 2, byrow = T))
  lims <- sort(c(max(object$resids[[series]][1:500,], na.rm = TRUE), 
            max(object$resids[[series]][1:500,], na.rm = TRUE) * -1 * 
              sign(max(object$resids[[series]][1:500,], na.rm = TRUE))))
  plot(as.vector(scale(ts_dat)),
       type = 'l', lwd = 2, col = 'darkred', bty = 'l',
       ylab = 'Value',
       main = levels(data_train$series)[series],
       ylim = lims,
       xaxt = 'n',
       xlab = '')
  time_axis()
  lincoefs <- vector(length = 500)
  for(i in 1:500){
    lines(object$resids[[series]][i,],
          col = '#D3D3D330',
          lwd = 0.75)
    lincoefs[i] <- coef(lm(as.vector(scale(ts_dat)) ~ object$resids[[series]][i,]))[2]
  }
  lines(as.vector(scale(ts_dat)),
        col = 'white', lwd = 3.5)
  lines(as.vector(scale(ts_dat)),
        col = 'darkred', lwd = 3)
  box(bty = 'l', lwd = 2)
  text(x = -5,
       y = lims[2],
       adj = 0,
       labels = 'NDVI moving average',
       col = 'darkred',
       xpd = TRUE)
  text(x = -5,
       y = lims[2] - 1.25,
       adj = 0,
       labels = 'Residual draws',
       col = 'grey50',
       xpd = TRUE)
  lims <- c(max(abs(lincoefs)), 
            max(abs(lincoefs)) * -1 * sign(max(abs(lincoefs))))
  hist(lincoefs, breaks = seq(min(lims),
                              max(lims),
                              length.out = 20), 
       col = 'darkred',
       border = 'white', xlab = bquote(beta[SOI]), 
       freq = FALSE,
       main = '',
       lwd = 2,
       ylab = 'Density')
  abline(v = 0, lwd = 3.5, col = 'white')
  abline(v = 0, lwd = 3, col = 'black')
  
  acfs <- vector(length = 500)
  for(i in 1:500){
    acfs[i] <- acf(object$resids[[series]][i,], lag.max = 12, plot = FALSE,
                   na.action = na.pass)$acf[[13]]
  }
  lims <- c(max(abs(acfs)), max(abs(acfs)) * -1 * sign(max(abs(acfs))))
  hist(acfs, breaks = seq(min(lims),
                          max(lims), 
                          length.out = 20), 
       col = 'darkred',
       border = 'white', xlab = 'Lag-12 ACF', 
       freq = FALSE,
       main = '',
       lwd = 2,
       ylab = '')
  abline(v = 0, lwd = 3.5, col = 'white')
  abline(v = 0, lwd = 3, col = 'black')
  layout(1)
}

#### Calculate DRPS score ####
drps_score <- function(truth, fc, interval_width = 0.9){
  nsum <- 1000
  Fy = ecdf(fc)
  ysum <- 0:nsum
  indicator <- ifelse(ysum - truth >= 0, 1, 0)
  score <- sum((indicator - Fy(ysum))^2)
  
  # Is value within empirical interval?
  interval <- quantile(fc, probs = c((1-interval_width)/2, (interval_width + (1-interval_width)/2)),
                       na.rm = TRUE)
  in_interval <- ifelse(truth <= interval[2] & truth >= interval[1], 1, 0)
  return(c(score, in_interval))
}

# Wrapper to operate on all observations in fc_horizon
drps_mcmc_object <- function(truth, fc, interval_width = 0.9){
  indices_keep <- which(!is.na(truth))
  if(length(indices_keep) == 0){
    scores = data.frame('drps' = rep(NA, length(truth)),
                        'interval' = rep(NA, length(truth)))
  } else {
    scores <- matrix(NA, nrow = length(truth), ncol = 2)
    for(i in indices_keep){
      scores[i,] <- drps_score(truth = as.vector(truth)[i],
                               fc = fc[,i], interval_width)
    }
  }
  scores
}


#### NDVI posterior contrast plot ####
plot_ndvi_contrast = function(object, series = 1, 
                              breaks = 30,
                              title = TRUE, 
                              xlabel = TRUE,
                              xlimits,
                              xaxis = TRUE,
                              line_labs = TRUE,
                              show_xlabs = TRUE){
  
  # Simulate NDVI from a green year vs 
  # a brown year to estimate a contrast
  
  # Create fake NDVI series with green summer
  scenario_1 <- rep(0.5, 200)
  
  # Set all other covariates to zero
  newdata <- data_train
  newdata <- lapply(seq_along(newdata), function(x){
    if(is.matrix(newdata[[x]])){
      matrix(0, ncol = NCOL(newdata[[x]]),
             nrow = 200)
    } else {
      rep(0, 200)
    }
  })
  names(newdata) <- names(data_train)
  newdata$series <- factor(rep(levels(data_train$series)[series], 200),
                           levels = levels(data_train$series))
  newdata$ndvi_ma12 <- scenario_1
  newdata$lag <- matrix(0:5, 200, 6, byrow = TRUE)
  
  # Generate posterior predictions for the new data
  Xp <- mvgam:::trend_Xp_matrix(newdata = newdata,
                                trend_map = object$trend_map,
                                series = 'all',
                                mgcv_model = object$trend_mgcv_model)
  trend_betas <- as.matrix(object, variable = 'trend_betas')
  test <- matrix(NA, nrow = NROW(trend_betas), ncol = NROW(Xp))
  for(i in 1:NROW(trend_betas)){
    test[i,] <- exp(Xp %*% trend_betas[i, ])
  }
  #test <- exp(suppressWarnings(predict(object, newdata = newdata, type = 'link')))
  
  # Repeat for a brown year
  scenario_2 <- rep(-0.5, 200)
  newdata$ndvi_ma12 <- scenario_2
  Xp <- mvgam:::trend_Xp_matrix(newdata = newdata,
                                trend_map = object$trend_map,
                                series = 'all',
                                mgcv_model = object$trend_mgcv_model)
  test2 <- matrix(NA, nrow = NROW(trend_betas), ncol = NROW(Xp))
  for(i in 1:NROW(trend_betas)){
    test2[i,] <- exp(Xp %*% trend_betas[i, ])
  }
  #test2 <- exp(suppressWarnings(predict(object, newdata = newdata,  type = 'link')))
  
  contrast_distribution <- as.vector(test - test2)
  
  # Plot histogram
  if(missing(xlimits)){
    hist_data <- hist(contrast_distribution,
                      breaks = breaks,
                      plot = F)
  } else {
    break_seq <- seq(xlimits[1], xlimits[2], 
                     by = (xlimits[2] - xlimits[1]) / breaks)
    
    if(!0 %in% break_seq){
      break_seq <- break_seq - (min(abs(break_seq)))
    }
    
    if(any(contrast_distribution < break_seq[1])){
      contrast_distribution <- contrast_distribution[!c(contrast_distribution < break_seq[1])]
    }
    
    if(any(contrast_distribution > tail(break_seq, 1))){
      contrast_distribution <- contrast_distribution[!c(contrast_distribution > tail(break_seq, 1))]
    }

    hist_data <- hist(contrast_distribution,
                      breaks = break_seq,
                      plot = F)
  }
  
  
  if(missing(xlimits)){
    xlimits <- c(-1 * max(abs(contrast_distribution),
                          na.rm = TRUE), max(abs(contrast_distribution),
                                             na.rm = TRUE))
  }
  ylimits <- range(hist_data$density)

  if(xlabel){
    xlab <- paste0('High - low NDVI contrast for ', 
                  levels(data_train$series)[series])
  } else {
    xlab <- ''
  }
  
  plot(1, type = "n", bty = 'n',
       ylab = '',
       xlab = xlab,
       xaxt = 'n',
       yaxt = 'n',
       xlim = xlimits,
       ylim = ylimits)
    hist(contrast_distribution, border = 'white',
         breaks = hist_data$breaks,
         freq = FALSE,
         lwd = 2, col = "#8F2727", add = TRUE)

  hist_data$density[which(hist_data$breaks>= 0)] <- 0
  plot(hist_data, freq = F, add = TRUE,
       border = 'white',
       col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5])
  
  text(x = xlimits[2] * 3/5,
       y = ylimits[2] * 4/5,
       labels = paste0(round(length(which(contrast_distribution > 0)) / 
                               length(contrast_distribution), 2) * 100, 
                       '%'),
       col = "#8F2727",
       pos = 4)
  text(x = xlimits[1] * 3/5,
       y = ylimits[2] * 4/5,
       labels = paste0(round( 1 - (length(which(contrast_distribution > 0)) / 
                               length(contrast_distribution)), 2) * 100, 
                       '%'),
       col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5],
       pos = 2)
  abline(v = 0, lwd = 3, col = 'white')
  abline(v = 0, lwd = 2.5, col = 'black')
  
  if(show_xlabs){
    axis(side = 1, lwd = 2, cex.axis = 0.8, tck= -0.08)
  } else {
    axis(side = 1, lwd = 2, labels = NA, tck= -0.08)
  }
  
}

#### Plot residuals over time ####
plot_resids_time = function(object, series){
  series_residuals <- object$resids[[series]]
  
  probs <- c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  cred <- sapply(1:NCOL(series_residuals),
                 function(n) quantile(series_residuals[,n],
                                      probs = probs, na.rm = TRUE))
  c_light <- c("#DCBCBC")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  
  plot(1, type = "n", bty = 'L',
       xaxt = 'n',
       xlab = '',
       ylab = '',
       xlim = c(0, NCOL(series_residuals)),
       ylim = range(cred, na.rm = TRUE))
  pred_vals <- 1:NCOL(cred)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(pred_vals, cred[5,], col = c_dark, lwd = 1)
}

#### Plot ACF functions for series residuals ####
plot_mod_acfs = function(object, series){

  probs <- c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  c_light <- c("#DCBCBC")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  
  series_residuals <- object$resids[[series]]
  acf1 <- acf(series_residuals[1,], plot = F,
              na.action = na.pass)
  resid_acf <- matrix(NA, nrow = NROW(series_residuals),
                      ncol = length(acf1$acf[,,1]))
  for(i in 1:NROW(series_residuals)){
    resid_acf[i, ] <- acf(series_residuals[i,], plot = F,
                          na.action = na.pass)$acf[,,1]
  }
  
  sorted_x <- seq(1:NCOL(resid_acf))
  N <- length(sorted_x)
  idx <- rep(1:N, each = 2)
  repped_x <- rep(sorted_x, each = 2)
  
  x <- sapply(1:length(idx),
              function(k) if(k %% 2 == 0)
                repped_x[k] + min(diff(sorted_x))/2 else
                  repped_x[k] - min(diff(sorted_x))/2)
  cred <- sapply(1:NCOL(resid_acf),
                 function(n) quantile(resid_acf[,n],
                                      probs = probs, na.rm = T))
  cred <- cred[, -1]
  clim <- qnorm((1 + .95)/2)/sqrt(acf1$n.used)
  plot(1, type = "n", bty = 'L',
       xlab = '',
       ylab = '',
       xlim = c(1, N-1),
       xaxt = 'n',
       ylim = range(c(cred,
                      -clim - 0.05,
                      clim + 0.05)))
  N <- N - 1
  rect(xleft = x[seq(1, N*2, by = 2)],
       xright = x[seq(2, N*2, by = 2)],
       ytop =  cred[9,],
       ybottom =  cred[1,],
       col = c_light,
       border = 'transparent')
  rect(xleft = x[seq(1, N*2, by = 2)],
       xright = x[seq(2, N*2, by = 2)],
       ytop =  cred[8,],
       ybottom =  cred[2,],
       col = c_light_highlight,
       border = 'transparent')
  rect(xleft = x[seq(1, N*2, by = 2)],
       xright = x[seq(2, N*2, by = 2)],
       ytop =  cred[7,],
       ybottom =  cred[3,],
       col = c_mid,
       border = 'transparent')
  rect(xleft = x[seq(1, N*2, by = 2)],
       xright = x[seq(2, N*2, by = 2)],
       ytop =  cred[6,],
       ybottom =  cred[4,],
       col = c_mid_highlight,
       border = 'transparent')
  
  for (k in 1:N) {
    lines(x = c(x[seq(1, N*2, by = 2)][k],x[seq(2, N*2, by = 2)][k]),
          y = c(cred[5,k], cred[5,k]),
          col = c_dark, lwd = 2)
  }
  abline(h = clim,  col = '#FFFFFF60', lwd = 2.85)
  abline(h = clim,  col = 'black', lwd = 2.5, lty = 'dashed')
  abline(h = -clim,  col = '#FFFFFF60', lwd = 2.85)
  abline(h = -clim, col = 'black', lwd = 2.5, lty = 'dashed')
}

#### Plot residual QQ norm functions ####
plot_mod_qq = function(object, series){
  
  probs <- c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  c_light <- c("#DCBCBC")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  
  series_residuals <- object$resids[[series]]
  coords <- qqnorm(series_residuals[1,], plot.it = F)
  resid_coords_y <- matrix(NA, nrow = NROW(series_residuals), ncol = length(coords$y))
  for(i in 1:NROW(series_residuals)){
    if(all(is.na(series_residuals[i,]))){
      resid_coords_y[i,] <- rep(NA, length(coords$y))
    } else {
      norm_coords <- qqnorm(series_residuals[i,], plot.it = FALSE)
      coords_y <- norm_coords$y
      coords_y[abs(coords_y) > 3.75] <- NA
      resid_coords_y[i,] <- coords_y[order(norm_coords$x)]
    }
  }
  
  cred <- sapply(1:NCOL(resid_coords_y),
                 function(n) quantile(resid_coords_y[,n],
                                      probs = probs,
                                      na.rm = TRUE))
  pred_vals <- coords$x[order(coords$x)]
  pred_vals <- pred_vals[complete.cases(cred[1,])]
  plot(x = pred_vals,
       y = cred[5,][complete.cases(cred[1,])],
       bty = 'L',
       xlab = '',
       xaxt = 'n',
       ylab = '',
       pch = 16,
       col = 'white',
       yaxt = 'n',
       cex = 1,
       xlim = c(-3.25, 3.25),
       ylim = c(-3.25, 3.25),
       tck = -0.04)

  polygon(c(pred_vals, rev(pred_vals)), c(cred[1,][complete.cases(cred[1,])],
                                          rev(cred[9,][complete.cases(cred[1,])])),
          col = c_light, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[2,][complete.cases(cred[1,])],
                                          rev(cred[8,][complete.cases(cred[1,])])),
          col = c_light_highlight, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[3,][complete.cases(cred[1,])],
                                          rev(cred[7,][complete.cases(cred[1,])])),
          col = c_mid, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[4,][complete.cases(cred[1,])],
                                          rev(cred[6,][complete.cases(cred[1,])])),
          col = c_mid_highlight, border = NA)
  lines(pred_vals, cred[5,][complete.cases(cred[1,])], col = c_dark, lwd = 1.5)
  qqline(cred[5,][complete.cases(cred[1,])], col = '#FFFFFF60', lwd = 2)
  qqline(cred[5,][complete.cases(cred[1,])], col = 'black', lwd = 1.5)
  box(bty = 'l', lwd = 2)
}

#### Plot standardised trend estimates for two series together ####
plot_trend_comp = function(object, series1, series2, hide_xlabels = FALSE){
  
  # Extract trend estimates for series 1
  data_train <- object$obs_data
  ends <- seq(0, dim(MCMCvis::MCMCchains(object$model_output, 'ypred'))[2],
              length.out = NCOL(object$ytimes) + 1)
  starts <- ends + 1
  starts <- c(1, starts[-c(1, (NCOL(object$ytimes)+1))])
  ends <- ends[-1]
  
  if(object$fit_engine == 'stan'){
    preds1 <- MCMCvis::MCMCchains(object$model_output, 'trend')[,seq(series1,
                                                                    dim(MCMCvis::MCMCchains(object$model_output,
                                                                                            'trend'))[2],
                                                                    by = NCOL(object$ytimes))]
  } else {
    preds1 <- MCMCvis::MCMCchains(object$model_output, 'trend')[,starts[series1]:ends[series1]]
  }
  
  preds1 <- (preds1 - mean(preds1)) / sd(preds1)
  
  # Extract trend estimates for series 2
  if(object$fit_engine == 'stan'){
    preds2 <- MCMCvis::MCMCchains(object$model_output, 'trend')[,seq(series2,
                                                                    dim(MCMCvis::MCMCchains(object$model_output,
                                                                                            'trend'))[2],
                                                                    by = NCOL(object$ytimes))]
  } else {
    preds2 <- MCMCvis::MCMCchains(object$model_output, 'trend')[,starts[series2]:ends[series2]]
  }
  
  preds2 <- (preds2 - mean(preds2)) / sd(preds2)
  
  # Plot quantiles
  probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  cred1 <- sapply(1:NCOL(preds1),
                 function(n) quantile(preds1[,n],
                                      probs = probs, na.rm = TRUE))
  cred2 <- sapply(1:NCOL(preds2),
                  function(n) quantile(preds2[,n],
                                       probs = probs, na.rm = TRUE))
  
  c_light <- c("#DCBCBC")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  pred_vals <- 1:NCOL(preds1)
  

  if(hide_xlabels){
    plot(1, type = "n", bty = 'L',
         ylab = 'Latent states (z-scored)',
         xlab = '',
         xaxt = 'n',
         xlim = c(0, length(pred_vals)),
         ylim = range(c(cred1, cred2)))
    
    box(bty = 'l', lwd = 2)
    
    polygon(c(pred_vals, rev(pred_vals)), c(cred1[2,], rev(cred1[8,])),
            col = c_light_highlight, border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred1[3,], rev(cred1[7,])),
            col = c_mid, border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred1[4,], rev(cred1[6,])),
            col = c_mid_highlight, border = NA)
    lines(pred_vals, cred1[5,], col = c_dark, lwd = 2.5)
    
    blues <- RColorBrewer::brewer.pal(n = 5, 'Blues')
    
    polygon(c(pred_vals, rev(pred_vals)), c(cred2[2,], rev(cred2[8,])),
            col = adjustcolor(blues[2], alpha.f = 0.4), border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred2[3,], rev(cred2[7,])),
            col = adjustcolor(blues[3], alpha.f = 0.6), border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred2[4,], rev(cred2[6,])),
            col = adjustcolor(blues[4], alpha.f = 0.8), border = NA)
    lines(pred_vals, cred2[5,], col = adjustcolor(blues[5], alpha.f = 0.9), lwd = 2.5)
  } else {
    plot(1, type = "n", bty = 'L',
         ylab = 'Latent states (z-scored)',
         xlab = '',
         xaxt = 'n',
         xlim = c(0, length(pred_vals)),
         ylim = range(c(cred1, cred2)))
    
    box(bty = 'l', lwd = 2)
    
    polygon(c(pred_vals, rev(pred_vals)), c(cred1[2,], rev(cred1[8,])),
            col = c_light_highlight, border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred1[3,], rev(cred1[7,])),
            col = c_mid, border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred1[4,], rev(cred1[6,])),
            col = c_mid_highlight, border = NA)
    lines(pred_vals, cred1[5,], col = c_dark, lwd = 2.5)
    
    blues <- RColorBrewer::brewer.pal(n = 5, 'Blues')
    
    polygon(c(pred_vals, rev(pred_vals)), c(cred2[2,], rev(cred2[8,])),
            col = adjustcolor(blues[2], alpha.f = 0.4), border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred2[3,], rev(cred2[7,])),
            col = adjustcolor(blues[3], alpha.f = 0.6), border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred2[4,], rev(cred2[6,])),
            col = adjustcolor(blues[4], alpha.f = 0.8), border = NA)
    lines(pred_vals, cred2[5,], col = adjustcolor(blues[5], alpha.f = 0.9), lwd = 2.5)
    time_axis()
  }

}

#### Plot constrained forecasts and calculate individual series DRPS ####
plot_fc_constrained = function(fcs, series, newdata, ylim,
                              ylab, xlab, hide_xlabels = FALSE,
                              colours = 'reds',
                              return_forecasts = FALSE,
                              realisations = FALSE, ...){
  
  preds <- fcs[[series]]
  s_name <- names(fcs)[series]
  
  if(!missing(newdata)){
    data_test <- newdata
  }
  
  # Plot quantiles of the forecast distribution
  preds_last <- preds[1,]
  probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  cred <- sapply(1:NCOL(preds),
                 function(n) quantile(preds[,n],
                                      probs = probs, na.rm = TRUE))
  if(colours == 'reds'){
    c_light <- c("#DCBCBC")
    c_light_highlight <- c("#C79999")
    c_mid <- c("#B97C7C")
    c_mid_highlight <- c("#A25050")
    c_dark <- c("#8F2727")
    c_dark_highlight <- c("#7C0000")
  } else {
    cols <- RColorBrewer::brewer.pal(n = 6, 'Blues')
    c_light <- cols[1]
    c_light_highlight <- cols[2]
    c_mid <- cols[3]
    c_mid_highlight <- cols[4]
    c_dark <- cols[5]
    c_dark_highlight <- cols[6]
  }
  
  if(missing(ylim)){
    ytrain <- data.frame(series = data_train$series,
                         time = data_train$time,
                         y = data_train$y) %>%
      dplyr::filter(series == s_name) %>%
      dplyr::select(time, y) %>%
      dplyr::distinct() %>%
      dplyr::arrange(time) %>%
      dplyr::pull(y)
    ylim <- c(min(cred, min(ytrain, na.rm = TRUE)),
              max(cred, max(ytrain, na.rm = TRUE)) + 2)
  }
  
  if(missing(ylab)){
    ylab <- paste0('Predicitons for ', levels(data_train$series)[series])
  }
  
  if(missing(xlab)){
    xlab <- 'Time'
  }
  
  pred_vals <- seq(1:length(preds_last))
  if(hide_xlabels){
    plot(1, type = "n", bty = 'L',
         xlab = '',
         xaxt = 'n',
         ylab = ylab,
         xlim = c(0, length(preds_last)),
         ylim = ylim, ...)
  } else {
    plot(1, type = "n", bty = 'L',
         xlab = xlab,
         ylab = ylab,
         xlim = c(0, length(preds_last)),
         ylim = ylim, ...)
  }
  
  if(realisations){
    for(i in 1:n_realisations){
      lines(x = pred_vals,
            y = preds[i,],
            col = 'white',
            lwd = 2.5)
      lines(x = pred_vals,
            y = preds[i,],
            col = sample(c("#DCBCBC",
                           "#C79999",
                           "#B97C7C",
                           "#A25050",
                           "#7C0000"), 1),
            lwd = 2.25)
    }
  } else {
    polygon(c(pred_vals, rev(pred_vals)), c(cred[1,], rev(cred[9,])),
            col = c_light, border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred[2,], rev(cred[8,])),
            col = c_light_highlight, border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred[3,], rev(cred[7,])),
            col = c_mid, border = NA)
    polygon(c(pred_vals, rev(pred_vals)), c(cred[4,], rev(cred[6,])),
            col = c_mid_highlight, border = NA)
    lines(pred_vals, cred[5,], col = c_dark, lwd = 2.5)
  }
  box(bty = 'L', lwd = 2)
  
  if(!missing(newdata)){
    
    if(class(data_train)[1] == 'list'){
      data_train <- data.frame(series = data_train$series,
                               y = data_train$y,
                               time = data_train$time)
      data_test <- data.frame(series = data_test$series,
                              y = data_test$y,
                              time = data_test$time)
    }
    
    last_train <- (NROW(data_train) / length(fcs))
    
    # Plot training and testing points
    points(dplyr::bind_rows(data_train, data_test) %>%
             dplyr::filter(series == s_name) %>%
             dplyr::select(time, y) %>%
             dplyr::distinct() %>%
             dplyr::arrange(time) %>%
             dplyr::pull(y), pch = 16, col = "white", cex = 0.8)
    points(dplyr::bind_rows(data_train, data_test) %>%
             dplyr::filter(series == s_name) %>%
             dplyr::select(time, y) %>%
             dplyr::distinct() %>%
             dplyr::arrange(time) %>%
             dplyr::pull(y), pch = 16, col = "black", cex = 0.65)
    abline(v = last_train, col = '#FFFFFF60', lwd = 2.85)
    abline(v = last_train, col = 'black', lwd = 2.5, lty = 'dashed')
    
    # Calculate out of sample DRPS and print the score
    drps_score <- function(truth, fc, interval_width = 0.9){
      nsum <- 1000
      Fy = ecdf(fc)
      ysum <- 0:nsum
      indicator <- ifelse(ysum - truth >= 0, 1, 0)
      score <- sum((indicator - Fy(ysum))^2)
      
      # Is value within empirical interval?
      interval <- quantile(fc, probs = c((1-interval_width)/2, (interval_width + (1-interval_width)/2)),
                           na.rm = TRUE)
      in_interval <- ifelse(truth <= interval[2] & truth >= interval[1], 1, 0)
      return(c(score, in_interval))
    }
    
    # Wrapper to operate on all observations in fc_horizon
    drps_mcmc_object <- function(truth, fc, interval_width = 0.9){
      indices_keep <- which(!is.na(truth))
      if(length(indices_keep) == 0){
        scores = data.frame('drps' = rep(NA, length(truth)),
                            'interval' = rep(NA, length(truth)))
      } else {
        scores <- matrix(NA, nrow = length(truth), ncol = 2)
        for(i in indices_keep){
          scores[i,] <- drps_score(truth = as.vector(truth)[i],
                                   fc = fc[,i], interval_width)
        }
      }
      scores
    }
    
    truth <- as.matrix(data_test %>%
                         dplyr::filter(series == s_name) %>%
                         dplyr::select(time, y) %>%
                         dplyr::distinct() %>%
                         dplyr::arrange(time) %>%
                         dplyr::pull(y))
    last_train <- length(data_train %>%
                           dplyr::filter(series == s_name) %>%
                           dplyr::select(time, y) %>%
                           dplyr::distinct() %>%
                           dplyr::arrange(time) %>%
                           dplyr::pull(y))
    
    fc <- preds[,(last_train+1):NCOL(preds)]
    
    
    if(all(is.na(truth))){
      message('No non-missing values in data_test$y; cannot calculate DRPS')
      message()
    } else {
      message('Out of sample DRPS:')
      print(sum(drps_mcmc_object(as.vector(truth),
                                 fc)[,1], na.rm = TRUE))
      message()
    }
    
  } else {
    if(class(data_train)[1] == 'list'){
      data_train <- data.frame(series = data_train$series,
                               y = data_train$y,
                               time = data_train$time)
    }
    
    points(data_train %>%
             dplyr::filter(series == s_name) %>%
             dplyr::select(time, y) %>%
             dplyr::distinct() %>%
             dplyr::arrange(time) %>%
             dplyr::pull(y),pch = 16, col = "white", cex = 0.8)
    points(data_train %>%
             dplyr::filter(series == s_name) %>%
             dplyr::select(time, y) %>%
             dplyr::distinct() %>%
             dplyr::arrange(time) %>%
             dplyr::pull(y),pch = 16, col = "black", cex = 0.65)
  }
  
  if(return_forecasts){
    if(!missing(newdata)){
      return(preds[,(last_train+1):NCOL(preds)])
    } else {
      return(preds)
    }
    
  }
}

#### Plot VAR coefficients in a more interpretable manner; ####
plot_var_coefs = function(object, filepath = 'Figures/VAR_betas.jpeg'){
  coef_estimates <- MCMCvis::MCMCchains(object$model_output, 'A')
  
  plot_coef_hist = function(coef_distribution, main,
                            show_xlabs = TRUE,
                            linecol = 'black',
                            borders = FALSE){
    coef_distribution[coef_distribution < -1.25] <- -1.25
    coef_distribution[coef_distribution > 1.25] <- 1.25
    hist_data <- hist(coef_distribution,
                      breaks = seq(-1.25, 1.25, length.out = 75),
                      plot = F)
    
    xlimits <- 1.1 * c(-1 * max(abs(coef_distribution),
                                na.rm = TRUE), max(abs(coef_distribution),
                                                   na.rm = TRUE))
    ylimits <- range(hist_data$density)
    
    plot(1, type = "n", bty = 'n',
         ylab = '',
         xaxt = 'n',
         yaxt = 'n',
         xlim = c(-1, 1),
         ylim = ylimits,
         main = main)
    
    if(borders){
      rect(xleft = -1, xright = 1, ybottom = -100, ytop = 100, col = 'grey80',
           border = NA)
    }
    
    line <- par(lwd = 0.4)
    hist(coef_distribution, border = 'white',
         breaks = hist_data$breaks,
         freq = FALSE,
         lwd = 0.5, col = "#8F2727", add = TRUE)
    
    hist_data$density[which(hist_data$breaks>= 0)] <- 0
    plot(hist_data, freq = F, add = TRUE,
         border = 'white',
         col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5],
         lwd = 0.5)
    if(borders){
      abline(v = 0, lwd = 1.5, col = 'grey80')
    } else {
      abline(v = 0, lwd = 1.5, col = 'white')
    }
    
    abline(v = 0, lwd = 1.2, col = linecol)
    if(show_xlabs){
      axis(side = 1, at = seq(-1, 1, by = 0.25), labels = c(-1, NA, -0.5, NA, 0, NA, 0.5, NA, 1),
           lwd = 1.2, cex.axis = 0.8, col = linecol, tck= -0.08)
    } else {
      axis(side = 1, at = seq(-1, 1, by = 0.25), lwd = 1.2, labels = NA, col = linecol, tck= -0.08)
    }
    
  }
  
  jpeg(filepath, width = 8, height = 5.5,
       res = 300, units = 'in')
  par(mgp=c(3,0.05,0), mar=c(0.9, 0, 0, 0),
      oma = c(0,0,0.7,0))
  layout(matrix(1:81, ncol = 9, nrow = 9, byrow = TRUE))
  for(x in 1:81){
    plot_coef_hist(coef_estimates[,x], 
                   main = '',
                   show_xlabs = x == 73,
                   linecol = 'black',
                   borders = x %in% seq(1, 81, by = 10))
    if(x %in% 1:9){
      title(main = levels(data_all$series)[which(1:9 == x)], 
            cex.main = 0.8, line = 0.1,
            xpd = NA)
    }
    
    if(x %in% seq(1, 81, by = 9)){
      title(ylab = levels(data_all$series)[which(seq(1, 81, by = 9) == x)], 
            cex.lab = 0.8, line = -1,
            font.lab = 2)
    }
    
  }
  dev.off()
  
}

#### Plot VAR covariances in a more interpretable manner; ####
plot_var_cors = function(object, filepath = 'Figures/VAR_cors.jpeg'){
  
  # Extract covariance estimates and convert to correlation matrices
  cov_estimates <- MCMCvis::MCMCchains(object$model_output, 'Sigma')
  cor_estimates <- matrix(NA, nrow = NROW(cov_estimates),
                          ncol = NCOL(cov_estimates))
  for(i in 1:NROW(cor_estimates)){
    cor_estimates[i, ] <- as.vector(cov2cor(matrix(cov_estimates[i,],
                                                   nrow = 9, ncol = 9)))
  }
  
  
  plot_coef_hist = function(coef_distribution, main,
                            show_xlabs = TRUE,
                            linecol = 'black',
                            borders = FALSE){
    coef_distribution[coef_distribution < -1] <- -1
    coef_distribution[coef_distribution > 1] <- 1
    hist_data <- hist(coef_distribution,
                      breaks = seq(-1, 1, length.out = 100),
                      plot = F)
    
    xlimits <- 1.1 * c(-1 * max(abs(coef_distribution),
                                na.rm = TRUE), max(abs(coef_distribution),
                                                   na.rm = TRUE))
    ylimits <- range(hist_data$density)
    
    plot(1, type = "n", bty = 'n',
         ylab = '',
         xaxt = 'n',
         yaxt = 'n',
         xlim = c(-1, 1),
         ylim = ylimits,
         main = main)
    
    if(!borders){
      line <- par(lwd = 0.4)
      hist(coef_distribution, border = 'white',
           breaks = hist_data$breaks,
           freq = FALSE,
           lwd = 0.5, col = "#8F2727", add = TRUE)
      
      hist_data$density[which(hist_data$breaks>= 0)] <- 0
      plot(hist_data, freq = F, add = TRUE,
           border = 'white',
           col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5],
           lwd = 0.5)
      if(borders){
        abline(v = 0, lwd = 1.5, col = 'grey80')
      } else {
        abline(v = 0, lwd = 1.5, col = 'white')
      }
      
      abline(v = 0, lwd = 1.2, col = linecol)
      if(show_xlabs){
        axis(side = 1, at = seq(-1, 1, by = 0.25), 
             labels = c(-1, NA, -0.5, NA, 0, NA, 0.5, NA, 1),
             lwd = 1.2, cex.axis = 0.8, col = linecol, tck= -0.08)
      } else {
        axis(side = 1, at = seq(-1, 1, by = 0.25), lwd = 1.2, 
             labels = NA, col = linecol, tck= -0.08)
      }
      
    }
    
  }
  
  jpeg(filepath, width = 8, height = 5.5,
       res = 300, units = 'in')
  par(mgp=c(3,0.05,0), mar=c(0.9, 0, 0, 0),
      oma = c(0,0,0.7,0))
  layout(matrix(1:81, ncol = 9, nrow = 9))
  for(x in 1:81){
    plot_coef_hist(cor_estimates[,x], 
                   main = '',
                   show_xlabs = x == 9,
                   linecol = 'black',
                   borders = x %in% seq(1, 81, by = 10))
    if(x %in% seq(1, 81, by = 9)){
      title(main = levels(data_all$series)[which(seq(1, 81, by = 9) == x)], 
            cex.main = 0.8, line = 0.1,
            xpd = NA)
    }
    
    if(x %in% 1:9){
      title(ylab = levels(data_all$series)[which(1:9 == x)], 
            cex.lab = 0.8, line = -1,
            font.lab = 2)
    }
    
  }
  dev.off()
  
}


#### Mintemp distributed lag plot ####
plot_mintemp_conditional = function(object, series, xlabel = TRUE,
                                 ylabel = TRUE){
  newdata <- lapply(seq_along(data_all), function(x){
    if(is.matrix(data_all[[x]])){
      matrix(0, nrow = 12, ncol = NCOL(data_all[[x]]))
    } else {
      rep(0, 12)
    }
  })
  names(newdata) <- names(data_all)
  newdata$series <- rep(levels(data_all$series)[series], 12)
  newdata$lag <- matrix(0:5, 12, 6, byrow = TRUE)
  newdata$mintemp <- data_all$mintemp[which(data_all$series == 'DM')[2:13],]
  if(levels(data_all$series)[series] == 'DM'){
    newdata$weights_dm <- matrix(1, 12, 6)
  }
  if(levels(data_all$series)[series] == 'DO'){
    newdata$weights_do <- matrix(1, 12, 6)
  }
  if(levels(data_all$series)[series] == 'OT'){
    newdata$weights_ot <- matrix(1, 12, 6)
  }
  if(levels(data_all$series)[series] == 'PP'){
    newdata$weights_pp <- matrix(1, 12, 6)
  }
  Xp <- mvgam:::trend_Xp_matrix(newdata = newdata,
                                trend_map = object$trend_map,
                                series = 'all',
                                mgcv_model = object$trend_mgcv_model)
  trend_betas <- as.matrix(object, variable = 'trend_betas')
  preds <- matrix(NA, nrow = NROW(trend_betas), ncol = NROW(Xp))
  for(i in 1:NROW(trend_betas)){
    preds[i,] <- Xp %*% trend_betas[i, ]
  }
  preds <- (preds - mean(preds, na.rm = TRUE)) / sd(preds, na.rm = TRUE)
  probs <- c(0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85)
  cred <- sapply(1:NCOL(preds),
                 function(n) quantile(preds[,n],
                                      probs = probs, na.rm = TRUE))
  cred <- cred - mean(cred[5,])
  
  c_light <- c("#DCBCBC")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  
  pred_vals <- 1:12
  if(ylabel){
    plot(1, type = "n", bty = 'L',
         xlab = '',
         xaxt = 'n',
         ylab = '',
         xlim = c(1, 12),
         ylim = c(-1.1 * max(abs(cred)), 1.1 * max(abs(cred))))
  } else {
    plot(1, type = "n", bty = 'L',
         xlab = '',
         xaxt = 'n',
         ylab = '',
         yaxt = 'n',
         xlim = c(1, 12),
         ylim = c(-1.1 * max(abs(cred)), 1.1 * max(abs(cred))))
    axis(side = 2, labels = NA)
  }

  polygon(c(pred_vals, rev(pred_vals)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(pred_vals, cred[5,], col = c_dark, lwd = 2.5)
  abline(h = 0, lwd = 3, col = 'white')
  abline(h = 0, lwd = 2.5)
  box(bty = 'L', lwd = 2)
  axis(side = 1, at = seq(1:12),
         labels = NA, cex.axis = 1,
         tck= -0.05)
  
if(xlabel){
  text(x=1:12, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=month.abb, srt=45, adj=1, xpd=NA)
}
  
  title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.25,
        xpd = NA)
  
}

#### VAR impulse response plots ####
# Calculate the forecast error variance decomposition, which is 
# based upon the orthogonalised impulse response coefficient 
# matrices Ψh and allow the user to analyse the contribution of 
# variable j to the h-step forecast error variance of variable k. 
# the result is a percentage figure
calc_vardecomps = function(object){
  beta_vars <- MCMCvis::MCMCchains(object$model_output, 'A')
  sigmas <- MCMCvis::MCMCchains(object$model_output, 'Sigma')
  all_decomps <- lapply(seq_len(NROW(beta_vars)), function(draw){
    
    # Get necessary VAR parameters into a simple list format
    x <- list(K = 9,
              A = matrix(beta_vars[draw,], 
                         nrow = 9, 
                         ncol = 9,
                         byrow = TRUE),
              Sigma = matrix(sigmas[draw,], 
                             nrow = 9, 
                             ncol = 9,
                             byrow = TRUE),
              p = 1)
    
    # Calculate the orthogonal variance decomposition
    ortho_vardecomp(x, h = 12, orthog = TRUE)
    
  })
  return(all_decomps)
}


# Calculate credible intervals of decompositions
get_decomp_creds = function(var_decomps, series1 = 1, series2 = 1){
  
  
  responses <- do.call(rbind, lapply(seq_len(length(var_decomps)),
                                     function(draw){
                                       var_decomps[[draw]][[series1]][, series2]
                                     }))

  probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  cred <- sapply(1:NCOL(responses),
                 function(n) quantile(responses[,n],
                                      probs = probs, na.rm = TRUE))
  cred
}

# Function to compute Generalized Impulse Response functions
# code modified from R code generously provided by Clinton Watkins:
# https://www.clintonwatkins.com/post/2021-generalised-impulse-response-function-r/#:~:text=This%20function%20is%20a%20replacement,bootstrap%20confidence%20intervals%20for%20IRFs.
gen_irf = function(x, h = 6, cumulative = TRUE, orthog = FALSE){
  
  spec_names <- species_names
  impulse <- spec_names
  response <- spec_names
  
  # Create arrays to hold calculations     
  # [1:nlags, 1:nvariables, shocked variable ] 
  IRF_o = array(data = 0, dim = c(h,x$K,x$K),
                 dimnames = list(NULL,spec_names,spec_names))       
  IRF_g = array(data = 0, dim = c(h,x$K,x$K),
                 dimnames = list(NULL,spec_names,spec_names))
  IRF_g1 = array(data = 0, dim = c(h,x$K,x$K))
  
  # Estimation of orthogonalised and generalised IRFs
  if(orthog){
    var_ma <- var_psi(x, h)
  } else {
    var_ma <- var_phi(x, h)
  }

  sigma.u <- x$Sigma
  P <- t(chol(sigma.u))
  sig_jj <- diag(sigma.u)
  
  for (jj in 1:x$K){
    indx_ <- matrix(0,x$K,1)
    indx_[jj,1] <- 1
    
    for (kk in 1:h){  #kk counts the lag
      IRF_o[kk, ,jj] <- var_ma[, ,kk]%*%P%*%indx_  # Peseran-Shin eqn 7 (OIRF)
      IRF_g1[kk, ,jj] <- var_ma[, ,kk]%*%sigma.u%*%indx_
      IRF_g[kk, ,jj] <- sig_jj[jj]^(-0.5)*IRF_g1[kk, ,jj]  # Peseran-Shin eqn 10 (GIRF)
      
    }
  }
  
  if(orthog==TRUE){
    irf <- IRF_o
  } else if(orthog==FALSE) {
    irf <- IRF_g
  } else {
    stop("\nError! Orthogonalised or generalised IRF?\n")
  }
  
  idx <- length(impulse)
  irs <- list()
  for (ii in 1:idx) {
    irs[[ii]] <- matrix(irf[1:(h), response, impulse[ii]], nrow = h)
    colnames(irs[[ii]]) <- response
    if (cumulative) {
      if (length(response) > 1) 
        irs[[ii]] <- apply(irs[[ii]], 2, cumsum)
      if (length(response) == 1) {
        tmp <- matrix(cumsum(irs[[ii]]))
        colnames(tmp) <- response
        irs[[ii]] <- tmp
      }
    }
  }
  names(irs) <- impulse
  result <- irs
  return(result)
  
}

# Convert a VAR A matrix to its moving average representation
var_phi = function(x, h = 10){
  h <- abs(as.integer(h))
  K <- x$K
  p <- x$p
  A <- as.array(x$A)
  if(h >= p){
    As <- array(0, dim = c(K, K, h + 1))
    for(i in (p + 1):(h + 1)){
      As[, , i] <- matrix(0, nrow = K, ncol = K)
    }
  } else {
    As <- array(0, dim = c(K, K, p))
  }
  As[, , 1] <- A
  Phi <- array(0, dim=c(K, K, h + 1))
  Phi[, ,1] <- diag(K)
  Phi[, , 2] <- Phi[, , 1] %*% As[, , 1]
  if (h > 1) {
    for (i in 3:(h + 1)) {
      tmp1 <- Phi[, , 1] %*% As[, , i-1]
      tmp2 <- matrix(0, nrow = K, ncol = K)
      idx <- (i - 2):1
      for (j in 1:(i - 2)) {
        tmp2 <- tmp2 + Phi[, , j+1] %*% As[, , idx[j]]
      }
      Phi[, , i] <- tmp1 + tmp2
    }
  }
  return(Phi)
}

# Convert a VAR A matrix to its orthogonalised moving average representation
var_psi = function(x, h=10){
    h <- abs(as.integer(h))
    Phi <- var_phi(x, h = h)
    Psi <- array(0, dim=dim(Phi))
    sigma.u <- x$Sigma
    P <- t(chol(sigma.u))
    dim3 <- dim(Phi)[3]
    for(i in 1:dim3){
      Psi[, , i] <- Phi[, , i] %*% P
    }
    return(Psi)
  }

# Forecast error decomposition
var_fecov = function(x, h) {
  sigma.u <- x$Sigma
  Sigma.yh <- array(NA, dim = c(x$K, x$K, h))
  Sigma.yh[, , 1] <- sigma.u
  Phi <- var_phi(x, h = h)
  if (h > 1) {
    for (i in 2:h) {
      temp <- matrix(0, nrow = x$K, ncol = x$K)
      for (j in 2:i) {
        temp <- temp + Phi[, , j] %*% sigma.u %*% t(Phi[, , j])
      }
      Sigma.yh[, , i] <- temp + Sigma.yh[, , 1]
    }
  }
  return(Sigma.yh)
}

ortho_vardecomp = function(x, h = 10, ...){

    h <- abs(as.integer(h))
    K <- x$K
    p <- x$p
    ynames <- species_names
    msey <- var_fecov(x, h = h)
    Psi <- var_psi(x, h = h)
    mse <- matrix(NA, nrow = h, ncol = K)
    Omega <- array(0, dim = c(h, K, K))
    for(i in 1 : h){
      mse[i, ] <- diag(msey[, , i])
      temp <- matrix(0, K, K)
      for(l in 1 : K){
        for(m in 1 : K){
          for(j in 1 : i){
            temp[l, m] <- temp[l, m] + Psi[l , m, j]^2
          }
        }
      }
      temp <- temp / mse[i, ]
      for(j in 1 : K){
        Omega[i, ,j] <- temp[j, ]
      }
    }
    result <- list()
    for(i in 1 : K){
      result[[i]] <- matrix(Omega[, , i], nrow = h, ncol = K)
      colnames(result[[i]]) <- ynames
    }
    names(result) <- ynames
    return(result)
  }

# Calculate impulse responses for each posterior draw
calc_irfs = function(object){
  beta_vars <- MCMCvis::MCMCchains(object$model_output, 'A')
  sigmas <- MCMCvis::MCMCchains(object$model_output, 'Sigma')
  all_irfs <- lapply(seq_len(NROW(beta_vars)), function(draw){
    
    # Get necessary VAR parameters into a simple list format
    x <- list(K = 9,
              A = matrix(beta_vars[draw,], 
                         nrow = 9, 
                         ncol = 9,
                         byrow = TRUE),
              Sigma = matrix(sigmas[draw,], 
                             nrow = 9, 
                             ncol = 9,
                             byrow = TRUE),
              p = 1)
    
    # Calculate the orthogonal irf for 12 steps ahead, initialising 
    # the shocks in a way that they can be directly compared
    gen_irf(x, h = 12, cumulative = FALSE, orthog = FALSE)
    
  })
  return(all_irfs)
}

plot_impulse_responses = function(all_irfs, series, filepath){
  
  # Extract IRFs for the specific series
  impulse_responses <- lapply(seq_along(all_irfs), function(j){
    all_irfs[[j]][series]
  })
  
  # Plot imp responses for all species apart from the impacted one
  c_light <- c("#DCBCBC")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  
  jpeg(filepath, width = 6.25, height = 4.25,
       res = 300, units = 'in')
  par(mar=c(2, 2, 1, 1),
      oma = c(1.25, 1.25, 0, 0))
  layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
  for(x in 1:9){
    if(x == series){
      plot(x = 1:12, y = 1:12, type = "n", bty = 'L',
           ylim = c(-2, 2),
           yaxt = 'n',
           xaxt = 'n',
           ylab = '',
           xlab = '')
      abline(h = 0, lwd = 2.5)
      box(bty = 'l', lwd = 2)
    } else {
      responses <- do.call(rbind, lapply(seq_along(impulse_responses), function(j){
        impulse_responses[[j]][[1]][,x]
      }))
      
      probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
      cred <- sapply(1:NCOL(responses),
                     function(n) quantile(responses[,n],
                                          probs = probs, na.rm = TRUE))
      pred_vals <- 1:12
      plot(1, type = "n", bty = 'L',
           xlab = '',
           xaxt = 'n',
           ylab = '',
           xlim = c(1, 12),
           ylim = c(-1.1 * max(abs(cred)), 1.1 * max(abs(cred))))
      polygon(c(pred_vals, rev(pred_vals)), c(cred[1,], rev(cred[9,])),
              col = c_light, border = NA)
      polygon(c(pred_vals, rev(pred_vals)), c(cred[2,], rev(cred[8,])),
              col = c_light_highlight, border = NA)
      polygon(c(pred_vals, rev(pred_vals)), c(cred[3,], rev(cred[7,])),
              col = c_mid, border = NA)
      polygon(c(pred_vals, rev(pred_vals)), c(cred[4,], rev(cred[6,])),
              col = c_mid_highlight, border = NA)
      lines(pred_vals, cred[5,], col = c_dark, lwd = 2.5)
      abline(h = 0, lwd = 3, col = 'white')
      abline(h = 0, lwd = 2.5)
      box(bty = 'L', lwd = 2)
    }
    if(x > 6){
      axis(1, cex.axis = 1,
           tck= -0.05)
    } else {
      axis(1, labels = NA, cex.axis = 1,
           tck= -0.05)
    }
    title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.35,
          xpd = NA)
    if(x == 4){
      title(ylab = bquote(Impulse~response~to~excess~captures~of~italic(.(species_names[series]))), 
            xpd = NA, line = 2.25)
    }
    if(x == 8){
      title(xlab = 'Forecast horizon (lunar months)', xpd = NA, line = 2.25)
    }
  }
  dev.off()
  
}

