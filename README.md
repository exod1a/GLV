
* Remember to set the working directory when running these analyses to be ecoevoanalysis.

This directory is for the running a Bayesian model in Stan to describe multi-species dynamics. The model is based on a generalised Lotka-Volterra model. The base for the code comes from the paper Gibson, T. and Gerber, G. (2018). *Robust and Scalable Models of Microbiome Dynamics* [Paper presentation]. Proceedings of the 35th International Conference on Machine Learning, Stockholm, Sweden, PMLR 80. http://proceedings.mlr.press/v80/gibson18a/gibson18a.pdf

### main.R: 
Runs entire analysis pipeline

### setup.R: 
Load packages and functions

### results_summary.Rmd: 
An overview of what has been done and what remains

### results_summary.html: 
The output of results_summary.Rmd

### gLVModel.stan 
The Stan code for the model

### input_data.R

### run_gLV.R

### interaction_matrices.R

### trajectories.R

### get_metadata.R

### posterior_checks.R
