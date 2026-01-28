data {
  int n, n_prey, n_predator;
  array[n] int prey, predator;
  vector[n] min_years_since, max_years_since;
}

parameters {
  vector[n_predator] alpha;
  vector<lower=0>[n_predator] sigma;
  vector<lower=0>[n_predator] tau;
  matrix[n_prey, n_predator] z;
  vector<lower=min_years_since, upper=max_years_since>[n] u;
}

transformed parameters {
  matrix[n_prey, n_predator] epsilon = diag_post_multiply(z, tau);
  vector[n] log_lik;
  {
    vector[n] mu = alpha[predator];
    for (i in 1:n) {
      mu[i] += epsilon[prey[i], predator[i]];
      log_lik[i] = normal_lpdf(u[i] | mu[i], sigma[predator[i]]);
    }
  }
  real lprior = std_normal_lpdf(alpha) + exponential_lpdf(sigma | 1)
                + exponential_lpdf(tau | 1) + std_normal_lpdf(to_vector(z));
}

model {
  target += sum(log_lik) + lprior;
}
