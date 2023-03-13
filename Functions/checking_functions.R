#### Checking functions for bespoke model evaluation ####

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
  
  ymax <- max(y) + 1
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
      axis(1, at = seq(2, 324,
                       by = 12), labels = seq(1997, 2023), cex.axis = 1,
           tck= -0.05)
    } else {
      axis(1, at = seq(2, 324,
                       by = 12), labels = NA, cex.axis = 1,
           tck= -0.05)
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
      title(xlab = 'Temporal lag (months)', xpd = NA, line = 2.25)
    }
  }
  dev.off()
}


#### Get constrained forecast distribution ####
get_fc_constained = function(object, 
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
evaluate_variogram = function(object, newdata, bound = 196){
  
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
  fcs <- get_fc_constained(object, bound = bound)
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
    out <- matrix(NA, length(obs), length(obs))
    for(i in 1:length(obs)){
      for(j in 1:length(obs)){
        if(i == j){
          out[i,j] <- 0
        } else {
          v_fc <- quantile(abs(fc_sample[i,] - fc_sample[j,]) ^ 0.5, 0.5)
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
  unlist(lapply(seq_len(fc_length), function(horizon){
    variogram_score(obs = truths[,horizon], 
                    fc_sample = fcs_per_horizon[[horizon]])   
  }))
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
  axis(1, at = seq(0, 336,
                   by = 12), labels = seq(1995, 2023), cex.axis = 1)
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
    axis(side = 1, lwd = 2, cex.axis = 0.8, tck = -0.08)
  } else {
    axis(side = 1, lwd = 2, labels = NA, tck = -0.08)
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
  axis(1, at = seq(0, 336,
                   by = 12), labels = seq(1995, 2023), cex.axis = 1)
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
  scenario_1 <- rep(0.75, 200)
  
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
  test <- exp(suppressWarnings(predict(object, newdata = newdata, type = 'link')))
  
  # Repeat for a brown year
  scenario_2 <- rep(-0.75, 200)
  newdata$ndvi_ma12 <- scenario_2
  test2 <- exp(suppressWarnings(predict(object, newdata = newdata,  type = 'link')))
  
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
  pred_vals <- 1:NCOL(leftover)
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
  box(bty = 'l', lwd = 2)

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
         ylab = 'Posterior trends',
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
         ylab = 'Posterior trends (standardised)',
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
    axis(1, at = seq(0, 336,
                     by = 12), labels = seq(1995, 2023), cex.axis = 1)
  }

}

#### Plot constrained forecasts and calculate individual series DRPS ####
plot_fc_constained = function(fcs, series, newdata, ylim,
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
# a cell represents the effect of the column on the row, so cell
# 3,6 represents the effect of species 6's trend at t-1 on species 3's
# trend at t
plot_var_coefs = function(object, filepath = 'Figures/VAR_betas.jpeg'){
  coef_estimates <- MCMCvis::MCMCchains(object$model_output, 'beta_var')
  
  plot_coef_hist = function(coef_distribution, main,
                            show_xlabs = TRUE,
                            linecol = 'black',
                            borders = FALSE){
    coef_distribution[coef_distribution < -1.25] <- -1.25
    coef_distribution[coef_distribution > 1.25] <- 1.25
    hist_data <- hist(coef_distribution,
                      breaks = seq(-1.25, 1.25, length.out = 100),
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
      rect(xleft = -1, xright = 1, ybottom = -10, ytop = 10, col = 'grey80',
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
      title(main = levels(data_all$series)[x], cex.main = 0.8, line = 0.1,
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
  preds <- predict(object, newdata = newdata, type = 'link')
  preds <- (preds - mean(preds, na.rm = TRUE)) / sd(preds, na.rm = TRUE)
  probs <- c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  cred <- sapply(1:NCOL(preds),
                 function(n) quantile(preds[,n],
                                      probs = probs, na.rm = TRUE))
  
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
# Calculate expected impulse variations when shocking each species
# separately; take proportions of these over a horizon of 6 to 
# calculate variance decompositions
calculate_var_decomps = function(object){
  beta_vars <- MCMCvis::MCMCchains(object$model_output, 'beta_var')
  
  var_decomps <- lapply(1:1000, function(x){
    all_imps <- lapply(1:9, function(series){
      imp_states <- rep(1, 9)
      imp_states[series] <- 1 + log(3)
      pulses <- iterate_var1(var_coefs = matrix(beta_vars[x,], 
                                                nrow = 9, 
                                                ncol = 9, 
                                                byrow = TRUE),
                             inits = imp_states)
      pulses <- rbind(rep(1, 9), pulses)
      abs(apply(pulses, 2, diff))
    })
    
    total_pulses <- Reduce('+', all_imps)
    lapply(1:9, function(series){
      all_imps[[series]] / total_pulses
    })
    
  })
  var_decomps
}

# Calculate credible intervals of decompositions
get_decomp_creds = function(var_decomps, series1 = 1, series2 = 1){
  responses <- do.call(rbind, lapply(var_decomps, function(matrix){
    matrix[[series]][, series2]
  }))
  probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  cred <- sapply(1:NCOL(responses),
                 function(n) quantile(responses[,n],
                                      probs = probs, na.rm = TRUE))
  cred
}

# Plot impulse responses
iterate_var1 = function(var_coefs, inits){
  states <- matrix(NA, nrow = 7, ncol = 9)
  states[1,] <- inits
  for(t in 2:7){
    states[t,] <- var_coefs %*% (states[t-1,] + rnorm(9, 0, 0.1))
  }
  states[2:7,]
}

plot_impulse_responses = function(object, series, filepath){
  beta_vars <- MCMCvis::MCMCchains(object$model_output, 'beta_var')
  
  impulse_responses <- lapply(1:1000, function(x){
    # Base prediction for 6 months ahead (all innovations are zero)
    base_preds <- iterate_var1(var_coefs = matrix(beta_vars[x,], 
                                                  nrow = 9, 
                                                  ncol = 9, 
                                                  byrow = TRUE),
                               inits = rep(1, 9))
    
    # Impulse response prediction (innovation for one species is three more captures
    # than expected from the GAM model)
    imp_states <- rep(1, 9)
    imp_states[series] <- 1 + log(3)
    imp_preds <- iterate_var1(var_coefs = matrix(beta_vars[x,], 
                                                 nrow = 9, 
                                                 ncol = 9, 
                                                 byrow = TRUE),
                              inits = imp_states)
    
    # Calculate deviations
    imp_preds - base_preds
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
      plot(x = 1:6, y = 1:6, type = "n", bty = 'L',
           ylim = c(-2, 2),
           yaxt = 'n',
           xaxt = 'n',
           ylab = '',
           xlab = '')
      abline(h = 0, lwd = 2.5)
      box(bty = 'l', lwd = 2)
    } else {
      responses <- do.call(rbind, lapply(impulse_responses, function(matrix){
        matrix[,x]
      }))
      
      probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
      cred <- sapply(1:NCOL(responses),
                     function(n) quantile(responses[,n],
                                          probs = probs, na.rm = TRUE))
      pred_vals <- 1:6
      plot(1, type = "n", bty = 'L',
           xlab = '',
           xaxt = 'n',
           ylab = '',
           xlim = c(1, 6),
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
      title(xlab = 'Forecast horizon (months)', xpd = NA, line = 2.25)
    }
  }
  dev.off()
  
}

