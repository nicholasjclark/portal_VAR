#### Analyse candidate models and produce model-based visualisations ####
library(mvgam)
library(scoringRules)
setwd("C:/Users/Nick/Google Drive/Academic Work Folder/Ecological forecasting/mv_portalcasting/rodent_evaluation_ms")
source('Functions/checking_functions.R')

# Load all relevant models
load('Outputs/bench1.rda')
load('Outputs/bench2.rda')
load('Outputs/modvar.rda')

# Load the pre-prepared modelling data
load('data/rodents_data_tsobjects.rda')

# Evaluate forecasts using the multivariate Variogram proper scoring rule;
# first 5 out of sample timepoints were all NA due to the COVID pandemic
bench1_scores <- evaluate_variogram(object = bench1, newdata = data_test)
bench1_scores <- bench1_scores[-c(1:5)]
bench1_lines <- predict(loess(y ~ time, span = 0.6, data = data.frame(y = log(bench1_scores),
                                                        time = seq_len(length(bench1_scores)))),
                        newdata = data.frame(time = seq_len(length(bench1_scores))))
bench2_scores <- evaluate_variogram(object = bench2, newdata = data_test)
bench2_scores <- bench2_scores[-c(1:5)]
bench2_lines <- predict(loess(y ~ time, span = 0.6, data=data.frame(y = log(bench2_scores),
                                                        time = seq_len(length(bench2_scores)))),
                        newdata = data.frame(time = seq_len(length(bench2_scores))))
modvar_scores <- evaluate_variogram(object = modvar, newdata = data_test)
modvar_scores <- modvar_scores[-c(1:5)]
modvar_lines <- predict(loess(y ~ time, span = 0.6, data=data.frame(y = log(modvar_scores),
                                                        time = seq_len(length(modvar_scores)))),
                        newdata = data.frame(time = seq_len(length(modvar_scores))))

jpeg('Figures/variogram_scores.jpg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
plot(1, type = "n", bty = 'L',
     xaxt = 'n',
     ylab = '',
     ylim = c(4, 6.8),
     xlim = c(-1.2, length(bench1_scores)),
     xlab = '')
lines(bench1_lines, lwd = 3.5, col = 'white')
lines(bench1_lines, lwd = 3, col = "#8F2727")
lines(bench2_lines, lwd = 3.5, col = 'white')
lines(bench2_lines, lwd = 3, col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5])
lines(modvar_lines, lwd = 3.5, col = 'white')
lines(modvar_lines, lwd = 3)

points(log(bench1_scores), pch = 16, cex = 1.2, col = 'white')
points(log(bench1_scores), pch = 16, col = "#8F2727")
points(log(bench2_scores), pch = 16, cex = 1.2, col = 'white')
points(log(bench2_scores), pch = 16, col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5])
points(log(modvar_scores), pch = 16, cex = 1.2, col = 'white')
points(log(modvar_scores), pch = 16)

axis(1, at = seq(1, 15, by = 2), labels = seq(1, 15, by = 2) + 5, cex.axis = 1,
     tck= -0.05)
box(bty = 'l', lwd = 2)
text(x = 0.8, y = bench2_lines[1], labels = 'AR',
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1,
     adj = 1)
text(x = 0.8, y = bench1_lines[1] + 0.1, labels = 'GAM-AR',
     col = "#8F2727", cex = 1, adj = 1)
text(x = 0.8, y = modvar_lines[1] - 0.06, labels = 'GAM-VAR',
      cex = 1, adj = 1)
mtext('Forecast horizon (months)', side = 1, outer = TRUE,
      line = 0.25, cex = 0.9)
mtext('log(Variogram score)', side = 2, outer = TRUE,
      line = 0.3, cex = 0.9)
dev.off()

# GAM-VAR is far superior; double check using unconstrained forecasts
sum(evaluate_variogram(object = bench1, newdata = data_test,
                       bound = 100000), na.rm = TRUE)
sum(evaluate_variogram(object = bench2, newdata = data_test,
                       bound = 100000), na.rm = TRUE)
sum(evaluate_variogram(object = modvar, newdata = data_test,
                       bound = 100000), na.rm = TRUE)

# GAM-VAR still superior
# Plot out of sample cumulative distribution functions for competing models
ecdf_plotdat = function(vals, x){
  if(length(which(is.na(vals))) > (length(vals) - 3)){
  } else {
    func <- ecdf(vals)
    func(x)
  }
}

# Use adjusted forecasts to compare predictive CDF
modvar_fcs <- get_fc_constained(modvar)
bench1_fcs <- get_fc_constained(bench1)
bench2_fcs <- get_fc_constained(bench2)

pred_ecdf = function(fc, plot_x, line_col, poly_col){
  pred_cdfs <- do.call(rbind, (lapply(1:NROW(fc), function(x){
    ecdf_plotdat(fc[x,], x = plot_x)
  })))
  
  probs <- c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  cred <- sapply(1:NCOL(pred_cdfs),
                 function(n) quantile(pred_cdfs[,n],
                                      probs = probs, na.rm = TRUE))
  
  polygon(c(plot_x, rev(plot_x)), c(cred[1,], rev(cred[9,])),
          col = poly_col, border = line_col)
}

jpeg('Figures/fc_cdfs.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  s_name <- levels(modvar$obs_data$series)[x]
  truths <- data.frame(y = data_test$y,
                       series = data_test$series, 
                       time = data_test$time) %>%
    dplyr::filter(series == s_name) %>%
    dplyr::select(time, y) %>%
    dplyr::distinct() %>%
    dplyr::arrange(time) %>%
    dplyr::pull(y)
  
  plot_x <- seq(min(truths, na.rm = T),
                max(truths, na.rm = T))
  
  plot(1, type = "n", bty = 'L',
       yaxt = 'n',
       xlab = '',
       ylab = '',
       xlim = c(min(plot_x, na.rm = TRUE), max(plot_x, na.rm = TRUE)),
       ylim = c(0, 1))
  
  # Add empirical quantiles of model CDFs
  pred_ecdf(fc = bench1_fcs[[x]][,tail(1:NCOL(modvar_fcs[[1]]),length(truths))], 
            plot_x = plot_x,
            line_col = "#8F2727",
            poly_col = scales::alpha("#8F2727",
                                     0.35))
  pred_ecdf(fc = bench2_fcs[[x]][,tail(1:NCOL(modvar_fcs[[1]]),length(truths))],
            plot_x = plot_x,
            line_col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5],
            poly_col = scales::alpha(RColorBrewer::brewer.pal(n = 5, 'Blues')[5],
                                     0.35))
  
  pred_ecdf(fc = modvar_fcs[[x]][,tail(1:NCOL(modvar_fcs[[1]]),length(truths))],
            plot_x = plot_x,
            line_col = 'black',
            poly_col = scales::alpha("black",
                                     0.35))
  
  lines(x = plot_x,
        y = ecdf_plotdat(truths,
                         plot_x),
        col = 'white',
        lwd = 2)
  lines(x = plot_x,
        y = ecdf_plotdat(truths,
                         plot_x),
        col = 'black',
        lwd = 1.5)
  box(bty = 'L', lwd = 2)
  
  title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.35,
        xpd = NA)
  if(x %in% c(1,4,7)){
    axis(side = 2)
  } else {
    axis(side = 2, labels = NA)
  }
  if(x == 4){
    title(ylab = 'Cumulative distribution function', xpd = NA, line = 2.25)
  }
  if(x == 8){
    title(xlab = 'Out of sample observed counts',
          xpd = NA, line = 2.25)
  }
}
dev.off()


# Inspect latent trend variance estimates 
# for the various models and compare some trends
jpeg('Figures/trend_sigmas.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 0, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  bench1_sigmas <- MCMCvis::MCMCchains(bench1$model_output, 
                                       paste0('sigma[', x, ']'),
                                       ISB = FALSE)
  bench2_sigmas <- MCMCvis::MCMCchains(bench2$model_output, 
                                       paste0('sigma[', x, ']'),
                                       ISB = FALSE)
  modvar_sigmas <- MCMCvis::MCMCchains(modvar$model_output, 
                                       paste0('sigma[', x, ']'),
                                       ISB = FALSE)
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
                     hist_data$density2,
                     hist_data$density3))
  plot(1, type = "n", bty = 'n',
       ylab = '',
       xlab = '',
       xaxt = 'n',
       yaxt = 'n',
       xlim = c(0, 1),
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
    title(xlab = expression(paste("Trend standard deviation (sqrt(", Sigma['var[i,i]'], "))")),
          xpd = NA, line = 2.25)
  }
}
dev.off()

plot(bench1, 'trend', series = 1, newdata = data_test)
plot(bench2, 'trend', series = 1, newdata = data_test)
plot(modvar, 'trend', series = 1, newdata = data_test)

plot(bench1, 'trend', series = 2, newdata = data_test)
plot(bench2, 'trend', series = 2, newdata = data_test)
plot(modvar, 'trend', series = 2, newdata = data_test)

plot(bench1, 'trend', series = 8, newdata = data_test)
plot(bench2, 'trend', series = 8, newdata = data_test)
plot(modvar, 'trend', series = 8, newdata = data_test)

# Plot residuals for all series to look for any unmodelled 
# systematic variation
jpeg('Figures/resids_time.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
object = modvar
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  plot_resids_time(object = modvar, series = x)
  
  if(x > 6){
    axis(1, at = seq(2, 320,
                     by = 12), labels = seq(1997, 2023), cex.axis = 1,
         tck= -0.05)
  } else {
    axis(1, at = seq(2, 320,
                     by = 12), labels = NA, cex.axis = 1,
         tck= -0.05)
  }
  
  title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.35,
        xpd = NA)
  if(x == 4){
    title(ylab = 'Randomized quantile residuals', xpd = NA, line = 2.25)
  }
  
}
dev.off()

# Plot ACF functions for all series
jpeg('Figures/resid_acfs.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  plot_mod_acfs(object = modvar, series = x)
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
    title(xlab = 'Lag (months)',
          xpd = NA, line = 2.25)
  }
}
dev.off()

# Plot QQ normal functions for all series
jpeg('Figures/resid_qqnorms.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  plot_mod_qq(object = modvar, series = x)
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

# Inspect overdispersion estimates for competing models, on the log scale
jpeg('Figures/phis.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 0, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  bench1_rs <- MCMCvis::MCMCchains(bench1$model_output, 
                                       paste0('r[', x, ']'),
                                       ISB = FALSE)
  bench1_rs <- bench1_rs[bench1_rs < 1000]

  bench2_rs <- MCMCvis::MCMCchains(bench2$model_output, 
                                       paste0('r[', x, ']'),
                                       ISB = FALSE)
  bench2_rs <- bench2_rs[bench2_rs < 1000]

  modvar_rs <- MCMCvis::MCMCchains(modvar$model_output, 
                                       paste0('r[', x, ']'),
                                       ISB = FALSE)
  modvar_rs <- modvar_rs[modvar_rs < 1000]
  
  break_seq <- seq(min(c(bench1_rs, bench2_rs, modvar_rs)), 
                   max(c(bench1_rs, bench2_rs, modvar_rs)), 
                   length.out = floor(1 / log(quantile(modvar_rs, 0.5)) * 300))
  hist_data <- hist(bench1_rs,
                    breaks = break_seq,
                    plot = F)
  hist_data2 <- hist(bench2_rs,
                    breaks = break_seq,
                    plot = F)
  hist_data3 <- hist(modvar_rs,
                    breaks = break_seq,
                    plot = F)
  ylimits <- range(c(hist_data$density,
                     hist_data$density2,
                     hist_data$density3))
  plot(1, type = "n", bty = 'n',
       ylab = '',
       xlab = '',
       xaxt = 'n',
       yaxt = 'n',
       xlim = c(0, max(break_seq)),
       ylim = ylimits * 1.2)

  
  plot_line_hist(s = bench1_rs, 
                 min_val = min(break_seq),
                 max_val = max(break_seq),
                 delta = ((max(break_seq) - min(break_seq)) / 
                            length(break_seq)),
                 line_col = "#8F2727",
                 poly_col = scales::alpha("#8F2727",
                                          0.35))
  plot_line_hist(s = bench2_rs, 
                 min_val = min(break_seq),
                 max_val = max(break_seq),
                 delta = ((max(break_seq) - min(break_seq)) / 
                            length(break_seq)),
                 line_col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5],
                 poly_col = scales::alpha(RColorBrewer::brewer.pal(n = 5, 'Blues')[5],
                                          0.35))
  plot_line_hist(s = modvar_rs, 
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
    title(xlab = expression(paste("Overdispersion (", phi, ")")),
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

plot(bench1, 'forecast', series = 8, newdata = data_test)
plot(bench2, 'forecast', series = 8, newdata = data_test)
plot(modvar, 'forecast', series = 8, newdata = data_test)

#### Analyse and evaluate the GAM-VAR model ####
# First calculate expert-adjusted forecast distributions, constrained by 
# the total number of available traps (196)
cons_fcs <- get_fc_constained(modvar)

# Plot the forecasts, which will also show the Discrete Rank Probability Score
# for the out of sample period
plot_fc_constained(cons_fcs, series = 1, newdata = data_test)
plot_fc_constained(cons_fcs, series = 2, newdata = data_test)
plot_fc_constained(cons_fcs, series = 3, newdata = data_test)
plot_fc_constained(cons_fcs, series = 4, newdata = data_test)
plot_fc_constained(cons_fcs, series = 5, newdata = data_test)
plot_fc_constained(cons_fcs, series = 6, newdata = data_test)
plot_fc_constained(cons_fcs, series = 7, newdata = data_test)
plot_fc_constained(cons_fcs, series = 8, newdata = data_test)
plot_fc_constained(cons_fcs, series = 9, newdata = data_test)

# Plot the NDVI random slope distributions
plot(modvar, 're')

# Plot the distributed lag posterior medians as bivariate heatmaps
jpeg('Figures/dist_lags.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:6, ncol = 3, nrow = 2, byrow = TRUE))
plot_mvgam_smooth(modvar, smooth = 2)
plot_mvgam_smooth(modvar, smooth = 3)
plot_mvgam_smooth(modvar, smooth = 4)
plot_mvgam_smooth(modvar, smooth = 5)
plot_mvgam_smooth(modvar, smooth = 6)
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
  plot_ndvi_contrast(object = modvar, series = x,
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

# Plot trend + residual against NDVI to inspect any support for nonlinear functions
jpeg('Figures/NDVI_leftovers.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(1.5, 1.5, 1, 1),
    oma = c(1.75, 1.75, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  ndvi_retro_residual(object = modvar, 
                      series = x,
                     xlabel = FALSE,
                     show_xlabs = x > 6,
                     n_bins = 20)
  title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.35,
        xpd = NA)
  if(x == 8){
    title(xlab = 'NDVI moving average', xpd = NA, line = 2.25)
  }
  if(x == 4){
    title(ylab = 'Trend + residual', xpd = NA, line = 2.25)
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
  plot_mintemp_conditional(object = modvar, series = x,
                        xlabel = x > 6,
                        ylabel = x %in% c(1,4,7))
  if(x == 4){
    title(ylab = 'Conditional minimum temperature effect (scaled)',
          xpd = NA, line = 2.25)
  }
}
dev.off()


# Plot uncertainty contributions
jpeg('Figures/uncertainty_conts.jpeg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE))
for(x in 1:9){
  plot_mvgam_uncertainty(object = modvar, newdata = data_test,
                         series = x, legend_position = 'none')
  title(main = bquote(italic(.(species_names[x]))), cex.main = 1, line = 0.25,
        xpd = NA)
  if(x == 1){
    text(1, 0.2, label="GAM component", 
         pos = 4, col="white")
    text(1, 0.8, label="Trend component", 
         pos = 4, col="#7C0000")
  }
  if(x == 4){
    title(ylab = 'Contribution to forecast uncertainty', xpd = NA, line = 2.25)
  }
  if(x == 8){
    title(xlab = 'Forecast horizon (months)', xpd = NA, line = 2.25)
  }
}
dev.off()

# Plot AR and VAR coefficients
plot_var_coefs(object = modvar)

# Plot some of the interesting trend comparisons
rstan::stan_hist(modvar$model_output, c('beta_var[2,6]',
                                        'beta_var[6,2]'))
jpeg('Figures/DO_PE_trends.jpeg', width = 6.5, height = 6.5,
     res = 300, units = 'in')
par(mar=c(2.25, 4, 0.5, 1))
layout(matrix(1:3, nrow = 3, byrow = T))
plot_trend_comp(object = modvar, series1 = 2, series2 = 6,
                hide_xlabels = TRUE)
axis(1, at = seq(2, 324,
                 by = 12), labels = NA, cex.axis = 1,
     tck= -0.05)
abline(v = max(data_train$time), col = '#FFFFFF60', lwd = 2.85)
abline(v = max(data_train$time), col = 'black', lwd = 2.5, lty = 'dashed')
text(x = max(data_train$time) + 16, y = 2.4, labels = 'Validation')
arrows(x0 = max(data_train$time) + 6, x1 = max(data_train$time) + 24,
       y0 = 2, y1 = 2, lwd = 1, length = 0.08)
text(x = max(data_train$time) - 15, y = 2.4, labels = 'Training')
arrows(x0 = max(data_train$time) - 6, x1 = max(data_train$time) - 24,
       y0 = 2, y1 = 2, lwd = 1, length = 0.08)
text(x = 70, y = 1.8, labels = expression(italic(Dipodomys~ordii)),
     col = "#8F2727", cex = 1)
text(x = 76, y = -1.4, labels = expression(italic(Peromyscus~eremicus)),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1)
plot_fc_constained(cons_fcs, series = 2, newdata = data_test, 
                   hide_xlabels = TRUE,
                   ylab = 'Posterior predictions')
axis(1, at = seq(2, 324,
                 by = 12), labels = NA, cex.axis = 1,
     tck= -0.05)
plot_fc_constained(cons_fcs, series = 6, newdata = data_test, 
                   hide_xlabels = TRUE,
                   colours = 'blues',
                   ylab = 'Posterior predictions')
axis(1, at = seq(2, 324,
                 by = 12), labels = seq(1997, 2023), cex.axis = 1,
     tck= -0.05)
dev.off()

# Compare the above for the three models
jpeg('Figures/DO_PE_trends_allmods.jpeg', width = 6.5, height = 6.5,
     res = 300, units = 'in')
par(mar=c(2.25, 4, 0.75, 1))
layout(matrix(1:3, nrow = 3, byrow = T))
plot_trend_comp(object = modvar, series1 = 2, series2 = 6,
                hide_xlabels = TRUE)
title('GAM-VAR', adj = 0)
axis(1, at = seq(2, 324,
                 by = 12), labels = NA, cex.axis = 1,
     tck= -0.05)
abline(v = max(data_train$time), col = '#FFFFFF60', lwd = 2.85)
abline(v = max(data_train$time), col = 'black', lwd = 2.5, lty = 'dashed')
text(x = max(data_train$time) + 16, y = 2.4, labels = 'Validation')
arrows(x0 = max(data_train$time) + 6, x1 = max(data_train$time) + 24,
       y0 = 2, y1 = 2, lwd = 1, length = 0.08)
text(x = max(data_train$time) - 15, y = 2.4, labels = 'Training')
arrows(x0 = max(data_train$time) - 6, x1 = max(data_train$time) - 24,
       y0 = 2, y1 = 2, lwd = 1, length = 0.08)
text(x = 70, y = 1.8, labels = expression(italic(Dipodomys~ordii)),
     col = "#8F2727", cex = 1)
text(x = 76, y = -1.4, labels = expression(italic(Peromyscus~eremicus)),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1)
plot_trend_comp(object = bench2, series1 = 2, series2 = 6,
                hide_xlabels = TRUE)
title('GAM-AR', adj = 0)
axis(1, at = seq(2, 324,
                 by = 12), labels = NA, cex.axis = 1,
     tck= -0.05)
abline(v = max(data_train$time), col = '#FFFFFF60', lwd = 2.85)
abline(v = max(data_train$time), col = 'black', lwd = 2.5, lty = 'dashed')
plot_trend_comp(object = bench1, series1 = 2, series2 = 6,
                hide_xlabels = TRUE)
title('AR', adj = 0)
axis(1, at = seq(2, 324,
                 by = 12), labels = seq(1997, 2023), cex.axis = 1,
     tck= -0.05)
abline(v = max(data_train$time), col = '#FFFFFF60', lwd = 2.85)
abline(v = max(data_train$time), col = 'black', lwd = 2.5, lty = 'dashed')
dev.off()


rstan::stan_hist(modvar$model_output, c('beta_var[6,8]',
                                        'beta_var[8,6]'))
jpeg('Figures/PE_PP_trends.jpeg', width = 6.5, height = 6.5,
     res = 300, units = 'in')
par(mar=c(2.25, 4, 0.5, 1))
layout(matrix(1:3, nrow = 3, byrow = T))
plot_trend_comp(object = modvar, series1 = 6, series2 = 8,
                hide_xlabels = TRUE)
axis(1, at = seq(2, 324,
                 by = 12), labels = NA, cex.axis = 1,
     tck= -0.05)
abline(v = max(data_train$time), col = '#FFFFFF60', lwd = 2.85)
abline(v = max(data_train$time), col = 'black', lwd = 2.5, lty = 'dashed')
text(x = max(data_train$time) + 16, y = 2.35, labels = 'Validation')
arrows(x0 = max(data_train$time) + 6, x1 = max(data_train$time) + 24,
       y0 = 1.95, y1 = 1.95, lwd = 1, length = 0.08)
text(x = max(data_train$time) - 15, y = 2.35, labels = 'Training')
arrows(x0 = max(data_train$time) - 6, x1 = max(data_train$time) - 24,
       y0 = 1.95, y1 = 1.95, lwd = 1, length = 0.08)
text(x = 50, y = -2.7, labels = expression(italic(Perognathus~eremicus)),
     col = "#8F2727", cex = 1)
text(x = 44, y = 1.9, labels = expression(italic(Chaetodipus~penicillatus)),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1)
plot_fc_constained(cons_fcs, series = 6, newdata = data_test, 
                   hide_xlabels = TRUE,
                   ylab = 'Posterior predictions')
axis(1, at = seq(2, 324,
                 by = 12), labels = NA, cex.axis = 1,
     tck= -0.05)
plot_fc_constained(cons_fcs, series = 8, newdata = data_test, 
                   hide_xlabels = TRUE,
                   colours = 'blues',
                   ylab = 'Posterior predictions')
axis(1, at = seq(2, 324,
                 by = 12), labels = seq(1997, 2023), cex.axis = 1,
     tck= -0.05)
dev.off()



rstan::stan_hist(modvar$model_output, c('beta_var[4,7]',
                                        'beta_var[7,4]'))
jpeg('Figures/OT_PF_trends.jpeg', width = 6.5, height = 6.5,
     res = 300, units = 'in')
par(mar=c(2.25, 4, 0.5, 1))
layout(matrix(1:3, nrow = 3, byrow = T))
plot_trend_comp(object = modvar, series1 = 4, series2 = 7,
                hide_xlabels = TRUE)
axis(1, at = seq(2, 324,
                 by = 12), labels = NA, cex.axis = 1,
     tck= -0.05)
abline(v = max(data_train$time), col = '#FFFFFF60', lwd = 2.85)
abline(v = max(data_train$time), col = 'black', lwd = 2.5, lty = 'dashed')
text(x = max(data_train$time) + 16, y = 2.35, labels = 'Validation')
arrows(x0 = max(data_train$time) + 6, x1 = max(data_train$time) + 24,
       y0 = 1.95, y1 = 1.95, lwd = 1, length = 0.08)
text(x = max(data_train$time) - 15, y = 2.35, labels = 'Training')
arrows(x0 = max(data_train$time) - 6, x1 = max(data_train$time) - 24,
       y0 = 1.95, y1 = 1.95, lwd = 1, length = 0.08)
text(x = 80, y = 2.5, labels = expression(italic(Onychomys~torridus)),
     col = "#8F2727", cex = 1)
text(x = 44, y = -1.8, labels = expression(italic(Perognathus~flavus)),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1)
plot_fc_constained(cons_fcs, series = 4, newdata = data_test, 
                   hide_xlabels = TRUE,
                   ylab = 'Posterior predictions')
axis(1, at = seq(2, 324,
                 by = 12), labels = NA, cex.axis = 1,
     tck= -0.05)
plot_fc_constained(cons_fcs, series = 7, newdata = data_test, 
                   hide_xlabels = TRUE,
                   colours = 'blues',
                   ylab = 'Posterior predictions')
axis(1, at = seq(2, 324,
                 by = 12), labels = seq(1997, 2023), cex.axis = 1,
     tck= -0.05)
dev.off()


rstan::stan_hist(modvar$model_output, c('beta_var[1,3]',
                                        'beta_var[3,1]'))
jpeg('Figures/DM_OL_trends.jpeg', width = 6.5, height = 6.5,
     res = 300, units = 'in')
par(mar=c(2.25, 4, 0.5, 1))
layout(matrix(1:3, nrow = 3, byrow = T))
plot_trend_comp(object = modvar, series1 = 1, series2 = 3,
                hide_xlabels = TRUE)
axis(1, at = seq(2, 324,
                 by = 12), labels = NA, cex.axis = 1,
     tck= -0.05)
abline(v = max(data_train$time), col = '#FFFFFF60', lwd = 2.85)
abline(v = max(data_train$time), col = 'black', lwd = 2.5, lty = 'dashed')
text(x = max(data_train$time) + 16, y = 2.35, labels = 'Validation')
arrows(x0 = max(data_train$time) + 6, x1 = max(data_train$time) + 24,
       y0 = 1.95, y1 = 1.95, lwd = 1, length = 0.08)
text(x = max(data_train$time) - 15, y = 2.35, labels = 'Training')
arrows(x0 = max(data_train$time) - 6, x1 = max(data_train$time) - 24,
       y0 = 1.95, y1 = 1.95, lwd = 1, length = 0.08)
text(x = 60, y = 1.9, labels = expression(italic(Dipodomys~merriami)),
     col = "#8F2727", cex = 1)
text(x = 44, y = -2.7, labels = expression(italic(Onychomys~leucogaster)),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1)
plot_fc_constained(cons_fcs, series = 1, newdata = data_test, 
                   hide_xlabels = TRUE,
                   ylab = 'Posterior predictions')
axis(1, at = seq(2, 324,
                 by = 12), labels = NA, cex.axis = 1,
     tck= -0.05)
plot_fc_constained(cons_fcs, series = 3, newdata = data_test, 
                   hide_xlabels = TRUE,
                   colours = 'blues',
                   ylab = 'Posterior predictions')
axis(1, at = seq(2, 324,
                 by = 12), labels = seq(1997, 2023), cex.axis = 1,
     tck= -0.05)
dev.off()

# Plot impulse response functions for select species
plot_impulse_responses(object = modvar, 
                       series = 1, 
                       filepath = 'Figures/DM_imp_response.jpg')
plot_impulse_responses(object = modvar, 
                       series = 2, 
                       filepath = 'Figures/DO_imp_response.jpg')
plot_impulse_responses(object = modvar, 
                       series = 5, 
                       filepath = 'Figures/PB_imp_response.jpg')
plot_impulse_responses(object = modvar, 
                       series = 7, 
                       filepath = 'Figures/PF_imp_response.jpg')
plot_impulse_responses(object = modvar, 
                       series = 8, 
                       filepath = 'Figures/PP_imp_response.jpg')

# Plot variance decompositions for select species
var_decomps <- calculate_var_decomps(modvar)

# Plot for species 3, 4, 2 and 8
jpeg('Figures/var_decomps.jpg', width = 6.25, height = 4.25,
     res = 300, units = 'in')
par(mar=c(2, 2, 1, 1),
    oma = c(1.25, 1.25, 0, 0))
layout(matrix(1:4, ncol = 2, nrow = 2, byrow = TRUE))
series = 3
all_creds <- matrix(NA, nrow = 9, ncol = 6)
for(i in 1:9){
  all_creds[i,] <- get_decomp_creds(var_decomps, series, i)[5,]
}
normalise = function(x){
  x / sum(x)
}
all_creds <- apply(all_creds, 2, normalise)

pred_vals <- 1:6
plot(1, type = "n", bty = 'L',
     xlab = 'Forecast horizon (months)',
     xaxt = 'n',
     ylab = '',
     xlim = c(1, 6),
     ylim = c(min(all_creds) - 0.01, max(all_creds) + 0.01))
cols <- rep('grey80', 9)
cols[5] <- "#8F2727"
cols[8] <- RColorBrewer::brewer.pal(n = 5, 'Blues')[5]
for(i in 1:9){
  lines(pred_vals, all_creds[i, ], col = 'white', lwd = 1.5)
  lines(pred_vals, all_creds[i, ], col = cols[i], lwd = 1)
}
lines(pred_vals, all_creds[series, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[series, ], col = 'black', lwd = 3)

lines(pred_vals, all_creds[5, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[5, ], col = cols[5], lwd = 3)
lines(pred_vals, all_creds[8, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[8, ], col = cols[8], lwd = 3)

text(x = 2.1, y = .15, labels = bquote(italic(.(species_names[5]))),
     col = "#8F2727", cex = 1, srt = 350)
text(x = 3.55, y = .072, labels = bquote(italic(.(species_names[8]))),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1)

box(bty = 'L', lwd = 2)
title(main = bquote(italic(.(species_names[series]))), cex.main = 1, line = 0.35,
      xpd = NA)
axis(1, cex.axis = 1,
     labels = NA,
     tck= -0.05)

# Series 4
series = 4
all_creds <- matrix(NA, nrow = 9, ncol = 6)
for(i in 1:9){
  all_creds[i,] <- get_decomp_creds(var_decomps, series, i)[5,]
}
normalise = function(x){
  x / sum(x)
}
all_creds <- apply(all_creds, 2, normalise)

pred_vals <- 1:6
plot(1, type = "n", bty = 'L',
     xlab = 'Forecast horizon (months)',
     xaxt = 'n',
     ylab = '',
     xlim = c(1, 6),
     ylim = c(min(all_creds) - 0.01, max(all_creds) + 0.01))
cols <- rep('grey80', 9)
cols[7] <- "#8F2727"
cols[8] <- RColorBrewer::brewer.pal(n = 5, 'Blues')[5]
for(i in 1:9){
  lines(pred_vals, all_creds[i, ], col = 'white', lwd = 1.5)
  lines(pred_vals, all_creds[i, ], col = cols[i], lwd = 1)
}
lines(pred_vals, all_creds[series, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[series, ], col = 'black', lwd = 3)

lines(pred_vals, all_creds[7, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[7, ], col = cols[7], lwd = 3)
lines(pred_vals, all_creds[8, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[8, ], col = cols[8], lwd = 3)

text(x = 2.4, y = .08, labels = bquote(italic(.(species_names[7]))),
     col = "#8F2727", cex = 1)
text(x = 2.4, y = .135, labels = bquote(italic(.(species_names[8]))),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1,
     srt = 342)

box(bty = 'L', lwd = 2)
title(main = bquote(italic(.(species_names[series]))), cex.main = 1, line = 0.35,
      xpd = NA)
axis(1, cex.axis = 1,
     labels = NA,
     tck= -0.05)

# Series 2
series = 2
all_creds <- matrix(NA, nrow = 9, ncol = 6)
for(i in 1:9){
  all_creds[i,] <- get_decomp_creds(var_decomps, series, i)[5,]
}
normalise = function(x){
  x / sum(x)
}
all_creds <- apply(all_creds, 2, normalise)

pred_vals <- 1:6
plot(1, type = "n", bty = 'L',
     xlab = 'Forecast horizon (months)',
     xaxt = 'n',
     ylab = '',
     xlim = c(1, 6),
     ylim = c(min(all_creds) - 0.01, max(all_creds) + 0.01))
cols <- rep('grey80', 9)
cols[7] <- "#8F2727"
cols[9] <- RColorBrewer::brewer.pal(n = 5, 'Blues')[5]
for(i in 1:9){
  lines(pred_vals, all_creds[i, ], col = 'white', lwd = 1.5)
  lines(pred_vals, all_creds[i, ], col = cols[i], lwd = 1)
}
lines(pred_vals, all_creds[series, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[series, ], col = 'black', lwd = 3)

lines(pred_vals, all_creds[7, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[7, ], col = cols[7], lwd = 3)
lines(pred_vals, all_creds[9, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[9, ], col = cols[9], lwd = 3)

text(x = 2.99, y = .15, labels = bquote(italic(.(species_names[7]))),
     col = "#8F2727", cex = 1, srt = 352)
text(x = 2.75, y = .075, labels = bquote(italic(.(species_names[9]))),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1,
     srt = 7)

box(bty = 'L', lwd = 2)
title(main = bquote(italic(.(species_names[series]))), cex.main = 1, line = 0.35,
      xpd = NA)
axis(1, cex.axis = 1,
     tck= -0.05)

# Series 8
series = 8
all_creds <- matrix(NA, nrow = 9, ncol = 6)
for(i in 1:9){
  all_creds[i,] <- get_decomp_creds(var_decomps, series, i)[5,]
}
normalise = function(x){
  x / sum(x)
}
all_creds <- apply(all_creds, 2, normalise)

pred_vals <- 1:6
plot(1, type = "n", bty = 'L',
     xlab = 'Forecast horizon (months)',
     xaxt = 'n',
     ylab = '',
     xlim = c(1, 6),
     ylim = c(min(all_creds) - 0.01, max(all_creds) + 0.01))
cols <- rep('grey80', 9)
cols[9] <- "#8F2727"
cols[3] <- RColorBrewer::brewer.pal(n = 5, 'Blues')[5]
for(i in 1:9){
  lines(pred_vals, all_creds[i, ], col = 'white', lwd = 1.5)
  lines(pred_vals, all_creds[i, ], col = cols[i], lwd = 1)
}
lines(pred_vals, all_creds[series, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[series, ], col = 'black', lwd = 3)

lines(pred_vals, all_creds[9, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[9, ], col = cols[9], lwd = 3)
lines(pred_vals, all_creds[3, ], col = 'white', lwd = 3.5)
lines(pred_vals, all_creds[3, ], col = cols[3], lwd = 3)

text(x = 3, y = .159, labels = bquote(italic(.(species_names[9]))),
     col = "#8F2727", cex = 1)
text(x = 3, y = .082, labels = bquote(italic(.(species_names[3]))),
     col = RColorBrewer::brewer.pal(n = 5, 'Blues')[5], cex = 1,
     srt = 5)

box(bty = 'L', lwd = 2)
title(main = bquote(italic(.(species_names[series]))), cex.main = 1, line = 0.35,
      xpd = NA)
axis(1, cex.axis = 1,
     tck= -0.05)
mtext('Forecast horizon (months)', side = 1, outer = TRUE,
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
modvar_mus <- MCMCvis::MCMCchains(modvar$model_output, 'mus')
for(i in seq(1999, 2022, by = 3)){
  index <- unique(data_train$time[which(data_train$year == i & data_train$month == 7)])
  med_counts = data.frame(series = levels(data_train$series),
                          med_count = unlist(lapply(seq_len(length(levels(data_train$series))), function(x){
                            exp(quantile(modvar_mus[,seq(x,dim(modvar_mus)[2],
                                                         by = length(levels(data_train$series)))][,index], probs = 0.5))
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