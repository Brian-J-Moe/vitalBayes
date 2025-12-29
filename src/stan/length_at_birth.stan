// birth_probit.stan
// Binomial regression with probit link for length-at-birth (b50)
// Probit is preferred for birth models based on threshold-crossing interpretation
// where latent developmental readiness follows a normal distribution

data {
  int<lower=1> N;                    // Number of observations
  vector<lower=0>[N] length;         // Observed lengths
  array[N] int<lower=0, upper=1> status; // Birth status (0=embryo, 1=free-swimming)
  
  // Priors for b50 (log scale)
  real prior_b50_mu;
  real<lower=0> prior_b50_sigma;
  
  // Priors for slope (log scale for positivity)
  real prior_slope_mu;
  real<lower=0> prior_slope_sigma;
}

parameters {
  real log_b50;                      // Log length at 50% birth probability
  real log_slope;                    // Log slope (transition steepness)
}

transformed parameters {
  real<lower=0> b50 = exp(log_b50);
  real<lower=0> slope = exp(log_slope);
  
  // Linear predictor (probit scale)
  vector[N] eta = slope * (length - b50);
  
  // Probability of being free-swimming (probit link)
  vector<lower=0, upper=1>[N] p = Phi(eta);
}

model {
  // Priors
  log_b50 ~ normal(prior_b50_mu, prior_b50_sigma);
  log_slope ~ normal(prior_slope_mu, prior_slope_sigma);
  
  // Likelihood
  status ~ bernoulli(p);
}

generated quantities {
  // Log-likelihood for LOO
  vector[N] log_lik;
  
  // Posterior predictive
  vector[N] p_pred = p;
  array[N] int status_rep;
  
  // Derived quantities
  // For probit: inverse of Phi at p = 0.05 is approx -1.645, at p = 0.95 is approx 1.645
  real<lower=0> b05 = b50 - 1.645 / slope;  // Length at 5% birth prob
  real<lower=0> b95 = b50 + 1.645 / slope;  // Length at 95% birth prob
  real<lower=0> transition_width = b95 - b05;
  
  // Summary statistics
  real mean_p_embryo;
  real mean_p_freeswim;
  int n_embryo = 0;
  int n_freeswim = 0;
  real sum_p_embryo = 0;
  real sum_p_freeswim = 0;
  int n_correct = 0;
  real prop_correct_rep;
  
  for (i in 1:N) {
    log_lik[i] = bernoulli_lpmf(status[i] | p[i]);
    status_rep[i] = bernoulli_rng(p[i]);
    
    if (status_rep[i] == status[i]) {
      n_correct += 1;
    }
    
    if (status[i] == 0) {
      // Embryo
      n_embryo += 1;
      sum_p_embryo += p[i];
    } else {
      // Free-swimming
      n_freeswim += 1;
      sum_p_freeswim += p[i];
    }
  }
  
  // Mean predicted probability for each group
  // For embryos, expected p should be low; for free-swimming, expected p should be high
  mean_p_embryo = n_embryo > 0 ? sum_p_embryo / n_embryo : 0;
  mean_p_freeswim = n_freeswim > 0 ? sum_p_freeswim / n_freeswim : 0;
  prop_correct_rep = n_correct * 1.0 / N;
}
