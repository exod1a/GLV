# functions

# visualise interaction matrix
#
# x  = matrix, to be plotted
# xl = string, x-axis label
# yl = string, y-axis label
# ti = string, plot title
plotmat <- function(x, xl, yl, ti) {
  cols <- viridis(20)
  x <- as.matrix(t(apply(x, 2, rev)))
  levelplot(x, xlab = xl, ylab = yl, main = ti, 
            scales=list(x=list(at=seq(1,dim(x)[1],1)), 
                        y=list(at=seq(1,dim(x)[1],1))), 
            col.regions = cols, par.settings=list(panel.background=list(col="white")))
}

# help to set the limits for the posterior check graphs 
#
# a = numeric vector
# b = numeric vector
lims <- function(a, b) {
  lims <- c()
  lims[1] <- min(min(a$x), min(b$x))  # lower bound x
  lims[2] <- max(max(a$x), max(b$x))  # upper bound x
  lims[3] <- min(min(a$y), min(b$y))  # lower bound y
  lims[4] <- max(max(a$y), max(b$y))  # upperbound y
  
  return(lims)
}

# Run generalised Lotka-Volterra simulation 
#
# n_spec = integer, number of species
# A      = n_spec x n_spec matrix, interaction coefficients
# b      = n_spec x 1 vector, growth rates
# y      = n_spec x 1 vector, initial abundances
# time   = real, total simulation time
# dt     = positive real, time step 
gLV_sim <- function(n_spec, A = matrix(sample(c(-2:2), n_spec^2, replace = T), nrow = n_spec), 
                    b = runif(n_spec), y = runif(n_spec), time = 5, dt = 0.05) {
  
  # run simulation
  res <- glv(N = n_spec, A, b, y, tend = time, tstep = dt)
  
  # prepare data
  res %<>% t() %>% as.data.frame()
  names(res)      <- paste0("s", seq_len(n_spec))
  
  return(res)
}

# solve quadratic formula
# a, b, and c are the coefficients
# returns both complex roots. If repeated root, it only returns one number
quad <- function(a, b, c)
{
  a <- as.complex(a)
  answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
              (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
  if(all(Im(answer) == 0)) answer <- Re(answer)
  if(answer[1] == answer[2]) return(answer[1])
  answer
}

# build learned matrix from Stan model
# summary = Stan summary
# n_clust = integer, the number of clusters
# row     = integer, row of summary to use
get_matrix <- function(summary, n_clust, n_spec, rownum = 1) {
  # find proper indices of interaction coefficients in summary
  elems <- c((2*n_spec+1):(2*n_spec+n_clust^2-n_clust))
  # appropriately sized matrix
  int_matrix <- diag(0, nrow = n_clust, ncol = n_clust)
  
  t_left  <- 1
  t_right <- 0                   # upper triangular section, right and down
  
  b_left  <- n_clust^2 - n_clust
  b_right <- b_left + 1          # lower triangular section, up and right
  
  for (i in 1:(n_clust-1)) {
    t_right <- t_right + n_clust - i
    b_right <- b_right - (n_clust - i)
    # fill matrix with proper indices
    int_matrix[i, (i+1):n_clust] <- summary[elems[t_left:t_right], rownum]
    int_matrix[(i+1):n_clust, i] <- summary[elems[b_left:b_right], rownum]
    
    t_left <- t_right + 1
    b_left <- b_right - 1
  }
  return(int_matrix) 
}

# build learned matrix from Stan model
# mat     = Stan summary mean of b_raw
# n_clust = integer, the number of clusters
# row     = integer, row of summary to use
get_matrix2 <- function(mat, n_clust, n_spec) {
  # find proper indices of interaction coefficients in summary
  elems <- c((2*n_spec+1):(2*n_spec+n_clust^2-n_clust))
  # appropriately sized matrix
  int_matrix <- diag(0, nrow = n_clust, ncol = n_clust)
  
  t_left  <- 1
  t_right <- 0                   # upper triangular section, right and down
  
  b_left  <- n_clust^2 - n_clust
  b_right <- b_left + 1          # lower triangular section, up and right
  
  for (i in 1:(n_clust-1)) {
    t_right <- t_right + n_clust - i
    b_right <- b_right - (n_clust - i)
    # fill matrix with proper indices
    int_matrix[i, (i+1):n_clust] <- mat[t_left:t_right]
    int_matrix[(i+1):n_clust, i] <- mat[b_left:b_right]
    
    t_left <- t_right + 1
    b_left <- b_right - 1
  }
  return(int_matrix) 
}

# find the confidence interval of the gLV trajectories
#
# samples = results from the Stan model
# y       = vector of initial abundances
# b       = vector of growth rates
# time    = total simulation time
# dt      = time step
# cl      = list of vectors of varying lengths, the vectors hold the cluster assignments
abundance_CI <- function(samples, y, b, time, dt, cl) {
  # determine the number of species, samples, clusters, indices
  n_spec    <- dim(extract(samples, "a1")[["a1"]])[2]
  n_samples <- dim(extract(samples, "a1")[["a1"]])[1]
  c         <- dim(extract(samples, "b_raw")[["b_raw"]])[2]
  n_clust   <- quad(1, -1, -c)[1]
  summary   <- summary(samples)$summary
  counter   <- 0 # for simulations that don't converge
  
  # create a set of rows in data frame for each species
  species <- c()
  for (i in 1:n_spec) {
    species <- c(species, rep(paste("s", i, sep = ""), time/dt+1))
  }
  
  # for storing simulation values
  df <- data.frame("time" = rep(seq(0, time, dt), n_spec),
                   "species" = species)
  
  # loop over every sample
  for (i in 1:n_samples) {
    mat       <- extract(samples, "b_raw")[["b_raw"]][i, ]
    A         <- get_matrix2(mat, n_clust, n_spec)         # cluster interaction matrix sample
    self_int  <- extract(samples, "a2")[["a2"]][i, ]       # self-interactions
    A         <- build_matrix(cl, int_matrix, self_int)    # full interaction matrix
    b         <- extract(samples, "a1")[["a1"]][i, ]       # growth rates sample
    sim       <- gLV_sim(n_spec, A, b, y, time, dt)        # run simulation with samples
    
    if (dim(sim)[1] == (time/dt + 1)) {
      df[, i+2-counter] <- melt(sim, id.vars = 0)$value    # add simulation result to df column
    } else {
      counter <- counter + 1
    }
  }
  
  # data frame to hold the means, and confidence intervals
  stats <- data.frame("time" = rep(seq(0, time, dt), n_spec),
                      "species" = species, "mean" = NA, 
                      "lower_2.5" = NA, "upper_2.5" = NA)
  
  # calculate their values
  for (i in 1:nrow(df)) {
    stats[i, 3] <- rowMeans(df[i, 3:ncol(df)])
    err         <- 2 * sd(df[i, 3:ncol(df)])
    stats[i, 4] <- stats[i, 3] - err
    stats[i, 5] <- stats[i, 3] + err
  }
  
  return(stats)
}


# hierarchical clustering
# df        = time series data frame where the rows are times
# n_clusts  = integer, number of clusters for cutting the tree
# metric    = string, distance metric for dissimilarity matrix
hr_clust <- function(df, n_clusts, metric = "euclidean") {
  # get dissimilarity matrix
  dis <- dist(t(df), method = metric)
  # run clustering algorithm
  tree <- hclust(dis)
  # find cluster assignments
  cl_nums <- cutree(tree, k = n_clusts)
  
  return(cl_nums)
}

# interaction matrix
#A <- matrix(sample(c(-2:2), n_spec^2, replace = T), nrow = n_spec)
#diag(A) <- -3

# create interaction matrix for clusters
# cl     = list of vectors of varying lengths, the vectors hold the cluster assignments
# cl_int = matrix containing cluster interactions (diagonal is always zero), cl_int_ij = effect of j on i
# sf_int = vector, the self-interactions
build_matrix <- function(cl, cl_int, sf_int) {
  
  A <- diag(sf_int)
  
  for (i in 1:(length(cl)-1)) {
    for (j in (i+1):length(cl)) {
      # effect of i on j
      A[cl[[j]], cl[[i]]] <- cl_int[j, i]
      # effect of j on i
      A[cl[[i]], cl[[j]]] <- cl_int[i, j]
    }
  }
  
  return(A)
}

# set non-zero matrix elements to one
normalise_matrix <- function(mat) {
  
  n       <- dim(mat)[1]          # get number of rows
  new_mat <- melt(mat)[[3]]       # turn into vector
  
  for (i in 1:length(new_mat)) {
    if (new_mat[i] != 0)
      new_mat[i] <- 1             # if non-zero, set to 1
  }  
  
  return(matrix(new_mat, nrow = n))
}

# rearrange columns so that clusters are together
# mat = matrix
# cl  = list of vectors of varying lengths, cluster assignments
rearrange_matrix <- function(mat, cl) {
  
  
  
}



