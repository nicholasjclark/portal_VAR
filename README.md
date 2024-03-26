# Multi-species dependencies improve forecasts of population dynamics in a long-term monitoring study

## Summary
This repository contains R code to extract data and replicate analyses in the manuscript titled *Multi-species dependencies improve forecasts of population dynamics in a long-term monitoring study* (currently in review; a preprint of a previous version is hosted on *biorxiv* at the following DOI: [https://doi.org/10.32942/X2TS34](https://doi.org/10.32942/X2TS34))

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
Raw rodent capture and covariate data have already been downloaded from the latest version of the `portalr` database and extracted to the `data` directory. Data can be prepared for analysis / modeling following instructions in the `1.prep_data.R` script. Models are built using `mvgam` functionality in the `2.models.R` script. A working version of `cmdstanr` is required to condition models on observed data. This script will produce large model objects (stored as class `mvgam`) that unfortunately cannot be uploaded to `Github` due to their size. Analysis of these models is completed in the `3.analysis.R` script. Figures are produced throughout the workflow. All of these figures are stored in the `Figures` directory. 
