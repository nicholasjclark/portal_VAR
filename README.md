# Forecasting rodent community regime transitions with Dynamic Generalized Additive Models

## Summary
This repository contains R code to extract data and replicate analyses in the manuscript titled *Forecasting rodent population dynamics and community transitions with dynamic nonlinear models* (Clark et al 2023, currently on hosted on *biorxiv* at the following DOI: XXXX)

## Required R libraries
`portalcasting`  
`tidyverse`   
`forecast`  
`cmdstanr` (and `Cmdstan`)  
[`mvgam`](https://github.com/nicholasjclark/mvgam)  
`scoringRules`

## Workflow
Raw rodent capture and covariate data have already been downloaded from the latest version of the `portalr` database and extracted to the `data` directory. Data can be prepared for analysis / modeling following instructions in the `1.prep_data.R` script. Models are built using `mvgam` functionality in the `2.models.R` script. A working version of `cmdstanr` is required to condition models on observed data. This script will produce large model objects (stored as class `mvgam`) that unfortunately cannot be uploaded to `Github` due to their size. Analysis of these models is completed in the `3.analysis.R` script. Figures are produced throughout the workflow. All of these figures are stored in the `Figures` directory. 
