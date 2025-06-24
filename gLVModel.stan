// Predict Growth Rates and Interaction Coefficients of a generalised 
// Lotka-Volterra (gLV) system.

data {

  int<lower=1> N;               // number of species
  int<lower=2> K;               // number of time points
  matrix<lower=0>[K, N] x;      // Data: abundance values for T time points and N species
  vector[K] t;                  // Data: time values
  int<lower=1> n_cl;            // number of clusters
  int<lower=1> cl_gr[N];        // cluster assigments

  // Hyperparameters
  real<lower=0> nu_a;           // for the variance of the growth rates and self-interactions
  real<lower=0> nu_b;           // for the variance of (interspecies) interaction coefficients
  real<lower=0> nu_w;           // variance of random noise
  
}

transformed data {
  
  // find off diagonal indices of n_cl x n_cl matrix
  int<lower=1, upper=n_cl> idxs[n_cl*n_cl-n_cl, 2];
  int count = 1;
  
  for (i in 1:(n_cl-1)) {
    for (j in (i+1):n_cl) {
      idxs[count] = {i, j};
      idxs[n_cl*n_cl-n_cl-count+1] = {j, i};
      
      count += 1;
    }
  }
  
}

parameters {

  vector[N] a1;                 // growth rates
  vector[N] a2;                 // self-interactions
  vector[n_cl*n_cl-n_cl] b_raw; // indexes for interaction matrix (diagonals removed)
  
  real<lower=0> var_a;          // variance of growth rates and self-interactions
  real<lower=0> var_b;          // variance of interspecies interactions
  real<lower=0> var_w;          // variance of stochasticity

}

model {
  
  matrix[n_cl, n_cl] b;         // interaction matrix
  
  for (i in 1:(n_cl*n_cl-n_cl)) {
    b[idxs[i, 1], idxs[i, 2]] = b_raw[i];
    
    if (i <= n_cl) {
      b[i, i] = 0;              // bii's are not parameters 
    }
  }
  
  // Priors
  b_raw ~ normal(0, var_b);
  
  var_a ~ inv_chi_square(nu_a);
  var_b ~ inv_chi_square(nu_b); // nu in (0, inf)
  var_w ~ inv_chi_square(nu_w);

  a1 ~ normal(0, var_a);
  a2 ~ normal(0, var_a);

  // Likelihood
  for (i in 1:N) {
    for (k in 1:(K-1)) {
      real dt = t[k+1] - t[k];
      real ints = 0;
      
      for (j in 1:N) {
        ints += b[cl_gr[i], cl_gr[j]] * x[k, j];  // make sure they are the proper clusters in b 
      }
      
      x[k+1, i] ~ normal(x[k, i] + x[k, i]*(a1[i] + a2[i]*x[k, i] + ints)*dt, 0.01);  // put dt*var_w back
    }
  }
  
}

generated quantities {
  
  matrix[K, N] x_sim;           // simulated data from the posterior
  matrix[n_cl, n_cl] b;         // interaction matrix
  
  for (i in 1:(n_cl*n_cl-n_cl)) {
    b[idxs[i, 1], idxs[i, 2]] = b_raw[i];
    
    if (i <= n_cl) {
      b[i, i] = 0;              // bii's are not parameters 
    }
  }
  
  for (i in 1:N) {
    x_sim[1, i] = x[1, i];
    for (k in 1:(K-1)) {
      real dt = t[k+1] - t[k];
      real ints = 0;
      
      for (j in 1:N) {
        ints += b[cl_gr[i], cl_gr[j]] * x[k, j];
      }      
      
      x_sim[k+1, i] = normal_rng(x[k, i] + x[k, i]*(a1[i] + a2[i]*x[k, i] + ints)*dt, dt*var_w);
    }
  }
  
}
