
# get necessary functions
source("R/model/functions.R")

setwd("/Users/chandlerross/Desktop/ResearchTurku/ecoevoanalysis")

# relevant packages
R_packages <- c("reshape2", "dplyr", "devtools", "seqtime", "wesanderson", "magrittr", 
                "rstan", "viridisLite", "lattice", "cowplot", "igraph")

# load packages
for (p in R_packages) library(p, character.only = TRUE)
