# Beyond single-species models: leveraging multispecies forecasts to navigate the dynamics of ecological predictability

<img src="Figures/model_definition.png" width = 620 alt="Dynamic Generalized Additive Model for forecasting rodent capture time series"/>

## Summary
This repository contains R code to extract data and replicate analyses in the manuscript titled *Beyond single-species models: leveraging multispecies forecasts to navigate the dynamics of ecological predictability* (currently in In-Press; a preprint of a previous version is hosted on *biorxiv* at the following DOI: [https://doi.org/10.32942/X2TS34](https://doi.org/10.32942/X2TS34)). This work shows how to build and interrogate multivariate Dynamic Generalized Additive Models (DGAMs) that can simultaneously learn useful multispecies dependencies and shared environmental effects, while also producing reliable probabilistic forecasts (the core model of the paper is presented in its full mathematical form above).

## Required R libraries
`portalcasting`  
`tidyverse`   
`forecast`  
`cmdstanr` (and `Cmdstan`)  
[`mvgam`](https://github.com/nicholasjclark/mvgam)  
`scoringRules`  
`extraDistr`  
`nleqslv`

## Workflow
Raw rodent capture and covariate data have already been downloaded from the latest version of the `portalr` database and extracted to the `data` directory. Data can be prepared for analysis / modeling following instructions in the `1.prep_data.R` script. Models are built using the [`mvgam`](https://github.com/nicholasjclark/mvgam) R package in the `2.models.R` script. Note that a working version of `Stan` is required to condition models on observed data, along with either the `cmdstanr` or `rstan` interface. This script will produce large model objects (stored as class `mvgam`) that unfortunately cannot be uploaded to `Github` due to their size. Analysis of these models is completed in the `3.analysis.R` script. Figures are produced throughout the workflow. All of these figures are stored in the `Figures` directory. 

## Other `mvgam` resources

If you are interested in using the `mvgam` package to build similar analyses, you may be interested in looking over some of the other resources that are publically available. A series of <a href="https://nicholasjclark.github.io/mvgam/" target="_blank">vignettes cover data formatting, forecasting and several extended case studies of DGAMs</a>. A number of other examples have also been compiled:

- <a href="https://www.youtube.com/watch?v=0zZopLlomsQ"
  target="_blank">Ecological Forecasting with Dynamic Generalized Additive
  Models</a>
- <a href="https://ecogambler.netlify.app/blog/distributed-lags-mgcv/"
  target="_blank">Distributed lags (and hierarchical distributed lags)
  using <code>mgcv</code> and <code>mvgam</code></a>
- <a href="https://ecogambler.netlify.app/blog/vector-autoregressions/"
  target="_blank">State-Space Vector Autoregressions in
  <code>mvgam</code></a>
- <a href="https://www.youtube.com/watch?v=RwllLjgPUmM"
  target="_blank">Ecological Forecasting with Dynamic GAMs; a tutorial and
  detailed case study</a>
- <a href="https://ecogambler.netlify.app/blog/interpreting-gams/"
  target="_blank">How to interpret and report nonlinear effects from
  Generalized Additive Models</a>
- <a href="https://www.youtube.com/watch?v=_fnDz2Bz3h8"
  target="_blank">Introduction to Stan and Hamiltonian Monte Carlo</a>
- <a href="https://ecogambler.netlify.app/blog/phylogenetic-smooths-mgcv/"
  target="_blank">Phylogenetic smoothing using <code>mgcv</code></a>
- <a href="https://ecogambler.netlify.app/blog/time-varying-seasonality/"
  target="_blank">Incorporating time-varying seasonality in forecast
  models</a>

