
#-------------------------------------------------------------------------------------------
# setup model
# simulation time
time <- 3; dt <- 0.05
# species clusters
n_clust <- sample(2:4, 1)      # 2-4 clusters
# number of species
n_spec  <- 5*n_clust - 2
# cluster assignments 
cl      <- split(1:n_spec, sample(n_clust, n_spec, repl=TRUE))
cl_assignments <- rep(NA, n_spec)
for (i in 1:n_clust)
  cl_assignments[cl[[i]]] <- i
# cluster interaction matrix
cl_int  <- matrix(floor(rnorm(n_clust^2, mean=0, sd=1.5)), nrow = n_clust, ncol = n_clust)
diag(cl_int) <- 0
# self interactions
sf_int  <- rep(-5, n_spec)
# create matrix
A       <- build_matrix(cl, cl_int, sf_int)
# growth rates
b       <- runif(n_spec)
# inital abundances
y       <- runif(n_spec)

# Stan parameters
nu_a     <- 0.1
nu_b     <- 0.1
nu_w     <- 0.1
n_chains <- 4
if (n_clust < 4) {
  n_iter <- 2000
} else {
  n_iter <- 3000 
}
# model version
version  <- 1.5
# path for files
pathway  <- paste("R/model/runs/run", run_num, sep = "")
dir.create(pathway)
pathway  <- paste(pathway, "/", sep = "")
