// Predict Growth Rates and Interaction Coefficients of a generalised 
// Lotka-Volterra (gLV) system.

data {

  int<lower=1> N;              // number of species
  int<lower=1> N_clust;        // number of clusters
  int<lower=2> K;              // number of time points
  matrix<lower=0>[K, N] x;     // Data: abundance values for T time points and N species
  vector[K] t;                 // Data: time values

  // Hyperparameters
  real<lower=0> nu_a;          // for the variance of the growth rates and self-interactions
  real<lower=0> nu_b;          // for the variance of (interspecies) interaction coefficients
  real<lower=0> nu_w;          // variance of random noise
  
  real<lower=0> theta;        
  real<lower=0> rho;           // for stick-breaking process
  
}

parameters {

  vector[N] a1;                             // growth rates
  matrix[N_clust, N_clust] b;               // interaction matrix
  
  real<lower=0> var_a;                      // variance of growth rates and self-interactions
  real<lower=0> var_b;                      // variance of interspecies interactions
  real<lower=0> var_w;                      // variance of stochasticity
  
  real<lower=0> alpha;                      // beta function parameter
  vector<lower=0, upper=1>[N_clust-1] brks; // breaks for stick-breaking process
  
}

transformed parameters {
  
  simplex[N_clust] pi_c;
  {
    real plus = 0;
    pi_c[1] = brks[1];
    plus = pi_c[1];
    for (i in 2:(N_clust-1)) {
      pi_c[i] = (1 - plus) * brks[i];
      plus += pi_c[i];
    }
    pi_c[N_clust] = 1 - plus;
  }
  
}

model {
  
  // Priors
  var_a ~ inv_chi_square(nu_a);
  var_b ~ inv_chi_square(nu_b);       
  var_w ~ inv_chi_square(nu_w);

  alpha ~ gamma(theta, rho);
  brks ~ beta(1, alpha);              // for stick-breaking process

  a1 ~ normal(0, var_a);
  
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        b[i, j] ~ normal(0, var_b);
      } else {
        b[i, j] ~ normal(0, var_a);
      }
    }
  }

  // Likelihood
  for (i in 1:N) {
    for (k in 1:(K-1)) {
      real dt = t[k+1] - t[k];
      
      x[k+1, i] ~ normal(x[k, i] + x[k, i] * (a1[i] + dot_product(b[i, ], x[k, ])) * dt, dt * var_w);
    }
  }
  
}
