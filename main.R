
# get packages
source("R/model/setup.R")

# set seed for pipeline
seed <- 101
set.seed(seed)
# for saving files
run_num <- 26     # don't forget to increment 
# run posterior checks
run_pos_checks <- F

# set parameters to run the pipeline
source("R/model/input_data.R")
# make cluster plot and save to file
source("R/model/cluster_plot.R")
# run the simulation and set up data
source("R/model/run_gLV.R")
# get interaction matrices
source("R/model/ints_growth_rates.R")
# plot predicted trajectories with CI and forecasting
source("R/model/trajectories.R")
# create metadata sheet
source("R/model/get_metadata.R")
# run posterior checks
if (run_pos_checks)
  source("R/model/posterior_checks.R")
