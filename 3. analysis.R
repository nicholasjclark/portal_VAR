#### Analyse candidate models and produce model-based visualisations ####
library(mvgam)
library(scoringRules)
setwd("C:/Users/Nick/Google Drive/Academic Work Folder/Ecological forecasting/mv_portalcasting/rodent_evaluation_ms")

# Load the pre-prepared modelling data and visualisation functions
load('data/rodents_data_tsobjects.rda')
source('Functions/checking_functions.R')

# Load all relevant models
load('Outputs/bench1.rda')
load('Outputs/bench2.rda')
load('Outputs/modvar.rda')
load('Outputs/bench1_all.rda')
load('Outputs/bench2_all.rda')
load('Outputs/modvar_all.rda')
load('Outputs/modvar.rda')

# Load the exact leave-future-out cross-validation results and plot them
load('Outputs/roll_evaluations.rda')
tot_scores <- matrix(NA, nrow = 3, ncol = 6)
jpeg('Figures/variogram_scores.jpg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
layout(matrix(1:6, ncol = 3, byrow = TRUE))
par(mar = c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
for(i in 1:6){
  if(i == 6){
    # If this is the primary model fitted, we need to extract 
    # forecast scores
    bench1_scores <- exp(log(score(forecast(bench1), score = 'variogram',
                                   log = TRUE)$all_series$score[1:12] *
                               score(forecast(bench1), score = 'energy',
                                     log = TRUE)$all_series$score[1:12]))
    bench2_scores <- exp(log(score(forecast(bench2), score = 'variogram',
                                   log = TRUE)$all_series$score[1:12] *
                               score(forecast(bench2), score = 'energy',
                                     log = TRUE)$all_series$score[1:12]))
    modvar_scores <- exp(log(score(forecast(modvar), score = 'variogram',
                                   log = TRUE)$all_series$score[1:12] *
                               score(forecast(modvar), score = 'energy',
                                     log = TRUE)$all_series$score[1:12]))
    
  } else {
    # Other scores are already calculated
    bench1_scores <- exp(bench1_roll[[i]]$cmbn_score$score)
    bench2_scores <- exp(bench2_roll[[i]]$cmbn_score$score)
    modvar_scores <- exp(modvar_roll[[i]]$cmbn_score$score)
  }
  
  # Compute total scores
  tot_scores[1,i] <- sum(bench1_scores, na.rm = TRUE)
  tot_scores[2,i] <- sum(bench2_scores, na.rm = TRUE)
  tot_scores[3,i] <- sum(modvar_scores, na.rm = TRUE)
  
  # Create loess smooth lines
  bench1_lines <- loess_scores(bench1_scores)
  bench2_lines <- loess_scores(bench2_scores)
  modvar_lines <- loess_scores(modvar_scores)
  
  plot(1, type = "n", bty = 'L',
       ylab = '',
       ylim = c(1.6, 8),
       xlim = c(1, length(bench1_scores)),
       xlab = '',
       xaxt = 'n',
       yaxt = 'n',
       main = bquote(italic(T)~' = '~.(ifelse(i <= 5, evaluation_seq[i],
                                              273))))
  
  if(i > 3){
    axis(side = 1)
  } else {
    axis(side = 1, labels = NA)
  }
  
  if(i %in% c(1, 4)){
    axis(side = 2)
  } else {
    axis(side = 2, labels = NA)
  }
  
  lines(bench1_lines, lwd = 3.25, 
        col = 'white')
  lines(bench1_lines, lwd = 2.9, 
        col = "#8F2727")
  lines(bench2_lines, lwd = 3.25, 
        col = 'white')
  lines(bench2_lines, lwd = 2.9, 
        col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5])
  lines(modvar_lines, lwd = 3.25, 
        col = 'white')
  lines(modvar_lines, lwd = 2.9)
  
  points(log(bench1_scores), pch = 16, cex = 0.9, 
         col = 'white')
  points(log(bench1_scores), pch = 16, cex = 0.7, 
         col = "#8F2727")
  points(log(bench2_scores), pch = 16, cex = 0.9, 
         col = 'white')
  points(log(bench2_scores), pch = 16, cex = 0.7, 
         col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5])
  points(log(modvar_scores), pch = 16, cex = 0.9, 
         col = 'white')
  points(log(modvar_scores), pch = 16, cex = 0.7)
  
  box(bty = 'l', lwd = 2)
  if(i == 6){
    text(x = 3.7, y = bench2_lines[3] - 0.4, labels = 'AR',
         col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1,
         adj = 1)
    text(x = 4.5, y = bench1_lines[3] + 0.4, labels = 'GAM-AR',
         col = "#8F2727", cex = 1, adj = 1)
    text(x = 5, y = modvar_lines[3] - 0.5, labels = 'GAM-VAR',
         cex = 1, adj = 1)
  }
}
mtext('Forecast horizon (lunar months)', side = 1, outer = TRUE,
      line = 0.25, cex = 0.9)
mtext('Weighted variogram score', side = 2, outer = TRUE,
      line = 0.3, cex = 0.9)
dev.off()

# GAM-VAR provides lowest overall scores in four of six cv folds
apply(tot_scores, 2, which.min)

# GAM-VAR is superior
# Inspect process error estimates for the various models and compare predictions
jpeg('Figures/trend_sigmas.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 0, 0, 0))
modvar_all_sigmas <- MCMCvis::MCMCchains(modvar$model_output,
                                     'Sigma')
names_org <- expand.grid(as.character(1:9),
                         as.character(1:9))

layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  bench1_sigmas <- MCMCvis::MCMCchains(bench1$model_output, 
                                       paste0('sigma[', x, ']'),
                                       ISB = FALSE)
  bench2_sigmas <- MCMCvis::MCMCchains(bench2$model_output, 
                                       paste0('sigma[', x, ']'),
                                       ISB = FALSE)
  modvar_sigmas <- sqrt(modvar_all_sigmas[,  
                                          which(names_org[,1] == x & 
                                                  names_org[,2] == x)])
  break_seq <- seq(0, 1, length.out = 50)
  hist_data <- hist(bench1_sigmas,
                    breaks = break_seq,
                    plot = F)
  hist_data2 <- hist(bench2_sigmas,
                     breaks = break_seq,
                     plot = F)
  hist_data3 <- hist(modvar_sigmas,
                     breaks = break_seq,
                     plot = F)
  ylimits <- range(c(hist_data$density,
                     hist_data2$density,
                     hist_data3$density))
  plot(1, type = "n", bty = 'n',
       ylab = '',
       xlab = '',
       xaxt = 'n',
       yaxt = 'n',
       xlim = c(0.2, 1),
       ylim = ylimits * 1.2)
  
  plot_line_hist(s = bench1_sigmas, 
                 min_val = min(break_seq),
                 max_val = max(break_seq),
                 delta = ((max(break_seq) - min(break_seq)) / 
                            length(break_seq)),
                 line_col = "#8F2727",
                 poly_col = scales::alpha("#8F2727",
                               0.35))
  plot_line_hist(s = bench2_sigmas, 
                 min_val = min(break_seq),
                 max_val = max(break_seq),
                 delta = ((max(break_seq) - min(break_seq)) / 
                            length(break_seq)),
                 line_col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5],
                 poly_col = scales::alpha(RColorBrewer::brewer.pal(n = 5, 'Blues')[5],
                                     0.35))
  plot_line_hist(s = modvar_sigmas, 
                 min_val = min(break_seq),
                 max_val = max(break_seq),
                 delta = ((max(break_seq) - min(break_seq)) / 
                            length(break_seq)),
                 line_col = 'black',
                 poly_col = scales::alpha("black",
                                          0.35))

  if(x > 6 ){
    axis(side = 1, lwd = 2, cex.axis = 0.8, tck= -0.08)
  } else {
    axis(side = 1, lwd = 2, labels = NA, tck= -0.08)
  }
  
  title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.35,
        xpd = NA)
  if(x == 8){
    title(xlab = expression(paste("Process standard deviation (sqrt(", Sigma['[i,i]'], "))")),
          xpd = NA, line = 2.25)
  }
}
dev.off()


# Look at a few unconstrained forecasts and calculate univariate
# Discrete Rank Probability Scores
plot(bench1, 'forecast', series = 1, newdata = data_test)
plot(bench2, 'forecast', series = 1, newdata = data_test)
plot(modvar, 'forecast', series = 1, newdata = data_test)

plot(bench1, 'forecast', series = 2, newdata = data_test)
plot(bench2, 'forecast', series = 2, newdata = data_test)
plot(modvar, 'forecast', series = 2, newdata = data_test)

plot(bench1, 'forecast', series = 3, newdata = data_test)
plot(bench2, 'forecast', series = 3, newdata = data_test)
plot(modvar, 'forecast', series = 3, newdata = data_test)

plot(bench1, 'forecast', series = 4, newdata = data_test)
plot(bench2, 'forecast', series = 4, newdata = data_test)
plot(modvar, 'forecast', series = 4, newdata = data_test)

plot(bench1, 'forecast', series = 5, newdata = data_test)
plot(bench2, 'forecast', series = 5, newdata = data_test)
plot(modvar, 'forecast', series = 5, newdata = data_test)

plot(bench1, 'forecast', series = 6, newdata = data_test)
plot(bench2, 'forecast', series = 6, newdata = data_test)
plot(modvar, 'forecast', series = 6, newdata = data_test)

plot(bench1, 'forecast', series = 7, newdata = data_test)
plot(bench2, 'forecast', series = 7, newdata = data_test)
plot(modvar, 'forecast', series = 7, newdata = data_test)

plot(bench1, 'forecast', series = 8, newdata = data_test)
plot(bench2, 'forecast', series = 8, newdata = data_test)
plot(modvar, 'forecast', series = 8, newdata = data_test)

plot(bench1, 'forecast', series = 9, newdata = data_test)
plot(bench2, 'forecast', series = 9, newdata = data_test)
plot(modvar, 'forecast', series = 9, newdata = data_test)

#### Analyse and evaluate the GAM-VAR model ####
# Plot residuals for all series to look for any unmodelled 
# systematic variation
jpeg('Figures/resids_time.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  plot_resids_time(object = modvar_all, series = x)
  
  if(x > 6){
    time_axis()
  } else {
    time_axis(labels = FALSE)
  }
  
  title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.35,
        xpd = NA)
  if(x == 4){
    title(ylab = 'Randomized quantile residuals', xpd = NA, line = 2.25)
  }
  
}
dev.off()

# Plot residual ACF functions for all series
jpeg('Figures/resid_acfs.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  plot_mod_acfs(object = modvar_all, series = x)
  if(x > 6){
    axis(1, at = seq(1, 24, by = 2), tck = -0.05, cex.axis = 1)
  } else {
    axis(1, at = seq(1, 24, by = 2), labels = NA, tck = -0.05, cex.axis = 1)
  }
  
  title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.35,
        xpd = NA)
  if(x == 4){
    title(ylab = 'Residual autocorrelation function', xpd = NA, line = 2.25)
  }
  if(x == 8){
    title(xlab = 'Lag (lunar months)',
          xpd = NA, line = 2.25)
  }
}
dev.off()

# Plot residual QQ normal functions for all series
jpeg('Figures/resid_qqnorms.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  plot_mod_qq(object = modvar_all, series = x)
  if(x > 6){
    axis(1, tck = -0.05, cex.axis = 1)
  } else {
    axis(1, labels = NA, tck = -0.05, cex.axis = 1)
  }
  
  title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.35,
        xpd = NA)
  if(x %in% c(1,4,7)){
    axis(side = 2)
  } else {
    axis(side = 2, labels = NA)
  }
  if(x == 4){
    title(ylab = 'Sample quantiles', xpd = NA, line = 2.25)
  }
  if(x == 8){
    title(xlab = 'Theoretical quantiles',
          xpd = NA, line = 2.25)
  }
}
dev.off()

# Calculate expert-adjusted forecast distributions, constrained by 
# the total number of available traps (196)
cons_fcs <- get_fc_constrained(modvar)

# Plot the forecasts, which will also show the Discrete Rank Probability Score
# for the out of sample period
jpeg('Figures/all_fcs.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  plot_fc_constrained(cons_fcs, series = x, newdata = data_test, 
                      hide_xlabels = TRUE, ylab = '')
  title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.35,
        xpd = NA)
  if(x %in% c(1, 4, 7)){
    title(ylab = 'Posterior predictions', xpd = NA, line = 2.25)
  }
  
  if(x > 6){
    time_axis()
  } else {
    time_axis(labels = FALSE)
  }
}
dev.off()

# Plot the NDVI random slope distributions
plot(modvar_all, 're', trend_effects = TRUE)

# Plot the distributed lag posterior medians as bivariate heatmaps
jpeg('Figures/dist_lags.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:6, ncol = 3, nrow = 2, byrow = TRUE))
plot_mvgam_smooth(modvar_all, smooth = 1, trend_effects = TRUE)
plot_mvgam_smooth(modvar_all, smooth = 2, trend_effects = TRUE)
plot_mvgam_smooth(modvar_all, smooth = 3, trend_effects = TRUE)
plot_mvgam_smooth(modvar_all, smooth = 4, trend_effects = TRUE)
plot_mvgam_smooth(modvar_all, smooth = 5, trend_effects = TRUE)
mtext('Minimum temperature (scaled)', side = 1, outer = TRUE,
      cex = 0.8)
mtext('Lag (in months)', side = 2, outer = TRUE,
      cex = 0.8)
dev.off()

# Plot high vs low NDVI posterior contrast distributions
jpeg('Figures/NDVI_contrasts.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 0, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  plot_ndvi_contrast(object = modvar_all, series = x,
                     xlimits = c(-2, 2),
                     xlabel = FALSE,
                     show_xlabs = x > 6)
  title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.35,
        xpd = NA)
  if(x == 8){
    title(xlab = 'Expected change in captures with higher NDVI', xpd = NA, line = 2.25)
  }
}
dev.off()

# Plot mintemp conditional curves
jpeg('Figures/mintemp_conditionals.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(1.5, 1.5, 1, 1),
    oma = c(1.75, 1.75, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  plot_mintemp_conditional(object = modvar_all, series = x,
                           xlabel = x > 6,
                           ylabel = x %in% c(1,4,7))
  if(x == 4){
    title(ylab = 'Conditional minimum temperature effect (scaled)',
          xpd = NA, line = 2.25)
  }
}
dev.off()


# Plot AR and VAR coefficients
plot_var_coefs(object = modvar_all)

# Plot variance-covariance estimates (as correlation matrices)
plot_var_cors(object = modvar_all)

# Plot some of the interesting forecast comparisons
rstan::stan_hist(modvar_all$model_output, c('A[2,4]',
                                        'A[4,2]'))
jpeg('Figures/DO_PF_trends.jpeg', width = 6.5, height = 6.5,
     res = 300, units = 'in')
par(mar=c(2.25, 4, 0.5, 1))
layout(matrix(1:3, nrow = 3, byrow = T))
plot_trend_comp(object = modvar, series1 = 2, series2 = 7,
                hide_xlabels = TRUE)
time_axis(labels = FALSE)
abline(v = max(data_train$time), col = '#FFFFFF60', lwd = 2.85)
abline(v = max(data_train$time), col = 'black', lwd = 2.5, lty = 'dashed')
text(x = max(data_train$time) + 14, y = 2.05, labels = 'Testing')
arrows(x0 = max(data_train$time) + 6, x1 = max(data_train$time) + 24,
       y0 = 1.72, y1 = 1.72, lwd = 1, length = 0.08)
text(x = max(data_train$time) - 15, y = 2.05, labels = 'Training')
arrows(x0 = max(data_train$time) - 6, x1 = max(data_train$time) - 24,
       y0 = 1.72, y1 = 1.72, lwd = 1, length = 0.08)
text(x = 65, y = -1.45, labels = expression(italic(Dipodomys~ordii)),
     col = "#8F2727", cex = 1)
text(x = 48, y = 1.8, labels = expression(italic(Perognathus~flavus)),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1)
plot_fc_constrained(cons_fcs, series = 2, newdata = data_test, 
                   hide_xlabels = TRUE,
                   ylab = 'Posterior predictions')
time_axis(labels = FALSE)
plot_fc_constrained(cons_fcs, series = 7, newdata = data_test, 
                   hide_xlabels = TRUE,
                   colours = 'blues',
                   ylab = 'Posterior predictions')
time_axis()
dev.off()

# Compare the above for the three models
jpeg('Figures/DO_PF_trends_allmods.jpeg', width = 6.5, height = 6.5,
     res = 300, units = 'in')
par(mar=c(2.25, 4, 0.75, 1))
layout(matrix(1:3, nrow = 3, byrow = T))
plot_trend_comp(object = modvar, series1 = 2, series2 = 7,
                hide_xlabels = TRUE)
title('GAM-VAR', adj = 0)
time_axis(labels = FALSE)
abline(v = max(data_train$time), col = '#FFFFFF60', lwd = 2.85)
abline(v = max(data_train$time), col = 'black', lwd = 2.5, lty = 'dashed')
text(x = max(data_train$time) + 14, y = 2.05, labels = 'Testing')
arrows(x0 = max(data_train$time) + 6, x1 = max(data_train$time) + 24,
       y0 = 1.72, y1 = 1.72, lwd = 1, length = 0.08)
text(x = max(data_train$time) - 15, y = 2.05, labels = 'Training')
arrows(x0 = max(data_train$time) - 6, x1 = max(data_train$time) - 24,
       y0 = 1.72, y1 = 1.72, lwd = 1, length = 0.08)
text(x = 65, y = -1.25, labels = expression(italic(Dipodomys~ordii)),
     col = "#8F2727", cex = 1)
text(x = 48, y = 1.8, labels = expression(italic(Perognathus~flavus)),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1)
plot_trend_comp(object = bench1, series1 = 2, series2 = 7,
                hide_xlabels = TRUE)
title('GAM-AR', adj = 0)
time_axis(labels = FALSE)
abline(v = max(data_train$time), col = '#FFFFFF60', lwd = 2.85)
abline(v = max(data_train$time), col = 'black', lwd = 2.5, lty = 'dashed')
plot_trend_comp(object = bench2, series1 = 2, series2 = 7,
                hide_xlabels = TRUE)
title('AR', adj = 0)
time_axis()
abline(v = max(data_train$time), col = '#FFFFFF60', lwd = 2.85)
abline(v = max(data_train$time), col = 'black', lwd = 2.5, lty = 'dashed')
dev.off()

# Calculate Generalized Impulse Response Functions
all_irfs <- calc_irfs(modvar_all)

# Plot impulse response functions for select species
plot_impulse_responses(all_irfs = all_irfs, 
                       series = 1, 
                       filepath = 'Figures/DM_imp_response.jpg')
plot_impulse_responses(all_irfs = all_irfs, 
                       series = 2, 
                       filepath = 'Figures/DO_imp_response.jpg')
plot_impulse_responses(all_irfs = all_irfs,  
                       series = 5, 
                       filepath = 'Figures/PB_imp_response.jpg')
plot_impulse_responses(all_irfs = all_irfs,  
                       series = 8, 
                       filepath = 'Figures/PP_imp_response.jpg')


# Plot variance decompositions for select species
var_decomps <- calc_vardecomps(object = modvar_all)

# Inspect variance decompositions for a few species in more detail
jpeg('Figures/var_decomps.jpg', width = 6.25, height = 3.2,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:2, ncol = 2, nrow = 1, byrow = TRUE))
series = 7
all_creds <- matrix(NA, nrow = 9, ncol = 12)
for(i in 1:9){
  all_creds[i,] <- get_decomp_creds(var_decomps, series, i)[5,1:12]
}

pred_vals <- 1:12
plot(1, type = "n", bty = 'L',
     xlab = 'Forecast horizon (lunar months)',
     xaxt = 'n',
     ylab = '',
     xlim = c(1, 12),
     ylim = range(all_creds))
cols <- rep('grey80', 9)
cols[4] <- "#8F2727"
cols[8] <- RColorBrewer::brewer.pal(n = 5, 'Blues')[5]
for(i in 1:9){
  lines(pred_vals, all_creds[i, ], col = 'white', lwd = 1.5)
  lines(pred_vals, all_creds[i, ], col = cols[i], lwd = 1)
}
lines(pred_vals, all_creds[series, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[series, ], col = 'black', lwd = 3)

lines(pred_vals, all_creds[3, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[3, ], col = cols[4], lwd = 3)
lines(pred_vals, all_creds[7, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[7, ], col = cols[8], lwd = 3)

text(x = 5.10, y = .34, labels = bquote(italic(.(abbrev_names[3]))),
     col = "#8F2727", cex = 1)
text(x = 3.50, y = .12, labels = bquote(italic(.(abbrev_names[7]))),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1)

box(bty = 'L', lwd = 2)
title(main = bquote(italic(.(species_names[series]))), cex.main = 1, line = 0.35,
      xpd = NA)
axis(1, cex.axis = 1,
     tck= -0.05)

# Series 6
series = 6
all_creds <- matrix(NA, nrow = 9, ncol = 12)
for(i in 1:9){
  all_creds[i,] <- get_decomp_creds(var_decomps, series, i)[5,1:12]
}

pred_vals <- 1:12
plot(1, type = "n", bty = 'L',
     xlab = 'Forecast horizon (lunar months)',
     xaxt = 'n',
     yaxt = 'n',
     ylab = '',
     xlim = c(1, 12),
     ylim = c(0, 1))
cols <- rep('grey80', 9)
cols[4] <- "#8F2727"
cols[5] <- RColorBrewer::brewer.pal(n = 5, 'Blues')[5]
for(i in 1:9){
  lines(pred_vals, all_creds[i, ], col = 'white', lwd = 1.5)
  lines(pred_vals, all_creds[i, ], col = cols[i], lwd = 1)
}
lines(pred_vals, all_creds[series, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[series, ], col = 'black', lwd = 3)

lines(pred_vals, all_creds[3, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[3, ], col = cols[4], lwd = 3)
lines(pred_vals, all_creds[6, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[6, ], col = cols[5], lwd = 3)

text(x = 3.1, y = .155, labels = bquote(italic(.(abbrev_names[3]))),
     col = "#8F2727", cex = 1)
text(x = 4, y = .082, labels = bquote(italic(.(abbrev_names[6]))),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1)

box(bty = 'L', lwd = 2)
title(main = bquote(italic(.(species_names[series]))), cex.main = 1, line = 0.35,
      xpd = NA)
axis(1, cex.axis = 1,
     tck= -0.05)
axis(2, labels = NA)
mtext('Forecast horizon (lunar months)', side = 1, outer = TRUE,
      line = 0.25, cex = 0.9)
mtext('Proportional contribution to trend variance', side = 2, outer = TRUE,
      line = 0.3, cex = 0.9)
dev.off()


#### A final figure to illustrate changes in community composition over time ####
points_circle <- function(n = 196, alpha = 2, 
                          geometry = 'planar',
                          med_counts) {
  
  radius <- function(k,n,b) {
    ifelse(
      k > n-b,
      1,
      sqrt(k-1/2)/sqrt(n-(b+1)/2)
    )
  }
  b <- round(alpha*sqrt(n))
  phi <- (sqrt(5)+1)/2
  r <- radius(1:n,n,b)
  theta <- 1:n * ifelse(geometry[1] == 'geodesic', 360*phi, 2*pi/phi^2)
  
  points <- data.frame(x = r*cos(theta),
                       y = r*sin(theta)) %>%
    dplyr::arrange(x, y)
  
  point_vals <- sort(factor(sample(med_counts$series, size = n, replace = TRUE,
                                   prob = med_counts$med_count),
                            levels = levels(data_train$series)))
  cols <- RColorBrewer::brewer.pal(9, 'Paired')
  
  plot(points, pch = 16, bty = 'none', xaxt = 'n', yaxt = 'n',
       ylab = '', xlab = '',
       col =   cols[as.numeric(point_vals)])
}

# Grab posterior median predictions for specific timepoints and plot
jpeg('Figures/med_community_comps.jpg', width = 6.25, height = 2.85,
     res = 300, units = 'in')
par(mar=c(1, 1, 1.5, 1),
    oma = c(0, 0, 0, 0))
layout(matrix(1:8, nrow = 2, byrow = TRUE))
modvar_mus <- MCMCvis::MCMCchains(modvar_all$model_output, 'mus')
for(i in seq(1999, 2022, by = 3)){
  print(i)
  index <- min(unique(data_all$time[which(data_all$year == i & data_all$month == 7)]))
  med_counts = data.frame(series = levels(data_all$series),
                          med_count = unlist(lapply(seq_len(length(levels(data_all$series))), function(x){
                            exp(quantile(modvar_mus[,seq(x,dim(modvar_mus)[2],
                                                         by = length(levels(data_all$series)))][,index], probs = 0.5))
                          })))
  points_circle(med_counts = med_counts, n = 250)
  title(main = paste0(i-2, ' - ', i), cex.main = 0.85, line = 0.95,
        xpd = NA)
  
  if(i == 1999){
    mtext(bquote(italic(.("D. mirriami"))*phantom(.(" dominates"))), 
          col=RColorBrewer::brewer.pal(9, 'Paired')[1], adj=0.5, cex = 0.5,
          line = 0.1)
    mtext(bquote(phantom(italic(.("D. mirriami")))*" dominates"),  adj=0.5, cex = 0.5,
          line = 0.1)
  }
  
  if(i == 2002){
    mtext(bquote(italic(.("C. baileyi"))*phantom(.(" begins reshuffling community"))), 
          col=RColorBrewer::brewer.pal(9, 'Paired')[5], adj=0.5, cex = 0.5,
          line = 0.1)
    mtext(bquote(phantom(italic(.("C. baileyi")))*" begins reshuffling community"),  adj=0.5, cex = 0.5,
          line = 0.1)
  }
  if(i == 2005){
    mtext(bquote(Shift~to~"'weaker'"~competitors),  adj=0.5, cex = 0.5,
          line = 0.1)
    mtext(bquote("("*phantom(italic(.("D. ordii")))*phantom(.("; "))*phantom(italic(.("C. penicillatus")))*phantom(.(")"))), 
          adj=0.5, cex = 0.5,
          line = -0.4)
    mtext(bquote(phantom(.("("))*italic(.("D. ordii"))*phantom(.("; "))*phantom(italic(.("C. penicillatus")))*phantom(.(")"))), 
          col=RColorBrewer::brewer.pal(9, 'Paired')[2],
          adj=0.5, cex = 0.5,
          line = -0.4)
    mtext(bquote(phantom(.("("))*phantom(italic(.("D. ordii")))*"; "*phantom(italic(.("C. penicillatus")))*phantom(.(")"))), 
          adj=0.5, cex = 0.5,
          line = -0.4)
    mtext(bquote(phantom(.("("))*phantom(italic(.("D. ordii")))*phantom(.("; "))*italic(.("C. penicillatus"))*phantom(.(")"))), 
          col=RColorBrewer::brewer.pal(9, 'Paired')[8],
          adj=0.5, cex = 0.5,
          line = -0.4)
    mtext(bquote(phantom(.("("))*phantom(italic(.("D. ordii")))*phantom(.("; "))*phantom(italic(.("C. penicillatus")))*")"), 
          adj=0.5, cex = 0.5,
          line = -0.4)
  }
  if(i == 2008){
    mtext(bquote("High diversity pre-drought"), adj=0.5, cex = 0.5,
          line = 0.1)
  }
  
  if(i == 2011){
    mtext(bquote("Drought reshuffles community"), adj=0.5, cex = 0.5,
          line = 0.3)
    mtext(bquote(italic(.("C. penicillatus"))*phantom(.(" increases"))), 
          col=RColorBrewer::brewer.pal(9, 'Paired')[8], adj=0.5, cex = 0.5,
          line = -0.4)
    mtext(bquote(phantom(italic(.("C. penicillatus")))*" increases"),  adj=0.5, cex = 0.5,
          line = -0.4)
  }
  
  if(i == 2014){
    mtext(bquote("Vegetation re-establishes"), adj=0.5, cex = 0.5,
          line = 0.3)
    mtext(bquote(italic(.("D. merriami"))*phantom(.(" increases"))), 
          col=RColorBrewer::brewer.pal(9, 'Paired')[1], adj=0.5, cex = 0.5,
          line = -0.3)
    mtext(bquote(phantom(italic(.("D. merriami")))*" increases"),  adj=0.5, cex = 0.5,
          line = -0.3)
  }
  
  if(i == 2017){
    mtext(bquote("Vegetation stabilises"), adj=0.5, cex = 0.5,
          line = 0.3)
    mtext(bquote(italic(.("D. ordii"))*phantom(.(" and "))*phantom(italic(.("O. torridus")))*phantom(.(" increase"))), 
          col=RColorBrewer::brewer.pal(9, 'Paired')[2], adj=0.5, cex = 0.5,
          line = -0.3)
    mtext(bquote(phantom(italic(.("D. ordii")))*" and "*phantom(italic(.("O. torridus")))*phantom(.(" increase"))), 
          adj=0.5, cex = 0.5,
          line = -0.3)
    mtext(bquote(phantom(italic(.("D. ordii")))*phantom(.(" and "))*italic(.("O. torridus"))*phantom(.(" increase"))), 
          col=RColorBrewer::brewer.pal(9, 'Paired')[4],
          adj=0.5, cex = 0.5,
          line = -0.3)
    mtext(bquote(phantom(italic(.("D. ordii")))*phantom(.(" and "))*phantom(italic(.("O. torridus")))*" increase"), 
          adj=0.5, cex = 0.5,
          line = -0.3)
  }
  
  if(i == 2020){
    mtext(bquote("Dry period reduces abundances"), adj=0.5, cex = 0.5,
          line = 0.3)
    mtext(bquote(italic(.("C. penicillatus"))*phantom(.(" and "))*phantom(italic(.("D. merriami")))*phantom(.(" dominate"))), 
          col=RColorBrewer::brewer.pal(9, 'Paired')[8], adj=0.5, cex = 0.5,
          line = -0.4)
    mtext(bquote(phantom(italic(.("C. penicillatus")))*" and "*phantom(italic(.("D. merriami")))*phantom(.(" dominate"))), 
          adj=0.5, cex = 0.5,
          line = -0.4)
    mtext(bquote(phantom(italic(.("C. penicillatus")))*phantom(.(" and "))*italic(.("D. merriami"))*phantom(.(" dominate"))), 
          col=RColorBrewer::brewer.pal(9, 'Paired')[1], 
          adj=0.5, cex = 0.5,
          line = -0.4)
    mtext(bquote(phantom(italic(.("C. penicillatus")))*phantom(.(" and "))*phantom(italic(.("D. merriami")))*" dominate"), 
          adj=0.5, cex = 0.5,
          line = -0.4)
  }
  
}
dev.off()