// length_at_maturity_single.stan
// Single-sex binomial regression for length-at-50%-maturity (L50)
// Uses probit link function with lognormal priors on L50
// Probit chosen for threshold-crossing interpretation: latent developmental
// readiness is normally distributed, individual matures when readiness > 0

data {
  int<lower=1> N;                    // Number of observations
  vector<lower=0>[N] length;         // Observed lengths
  array[N] int<lower=0, upper=1> mature; // Maturity status (0=immature, 1=mature)
  
  // Priors for L50 (log scale)
  real prior_L50_mu;
  real<lower=0> prior_L50_sigma;
  
  // Priors for slope (log scale for positivity)
  real prior_slope_mu;
  real<lower=0> prior_slope_sigma;
}

parameters {
  real log_L50;                      // Log length at 50% maturity
  real log_slope;                    // Log slope (transition steepness)
}

transformed parameters {
  real<lower=0> L50 = exp(log_L50);
  real<lower=0> slope = exp(log_slope);
  
  // Linear predictor
  vector[N] eta = slope * (length - L50);
  
  // Probability of maturity (probit link)
  vector<lower=0, upper=1>[N] p = Phi(eta);
}

model {
  // Priors
  log_L50 ~ normal(prior_L50_mu, prior_L50_sigma);
  log_slope ~ normal(prior_slope_mu, prior_slope_sigma);
  
  // Likelihood
  mature ~ bernoulli(p);
}

generated quantities {
  // Log-likelihood for LOO
  vector[N] log_lik;
  
  // Posterior predictive
  vector[N] p_pred = p;
  array[N] int mature_rep;
  
  // Derived quantities (probit scale: use inv_Phi aka qnorm)
  // Phi^{-1}(0.05) = -1.645, Phi^{-1}(0.95) = 1.645
  real<lower=0> L05 = L50 - 1.6449 / slope;   // Length at 5% maturity
  real<lower=0> L95 = L50 + 1.6449 / slope;   // Length at 95% maturity
  real<lower=0> transition_width = L95 - L05;  // Width of transition zone
  
  // Summary statistics
  real mean_p_mature;
  real mean_p_immature;
  int n_mature = 0;
  int n_immature = 0;
  real sum_p_mature = 0;
  real sum_p_immature = 0;
  int n_correct = 0;
  real prop_correct_rep;
  
  for (i in 1:N) {
    log_lik[i] = bernoulli_lpmf(mature[i] | p[i]);
    mature_rep[i] = bernoulli_rng(p[i]);
    
    // Track classification accuracy
    if (mature_rep[i] == mature[i]) {
      n_correct += 1;
    }
    
    // Accumulate for calibration checks
    if (mature[i] == 1) {
      n_mature += 1;
      sum_p_mature += p[i];
    } else {
      n_immature += 1;
      sum_p_immature += p[i];
    }
  }
  
  mean_p_mature = n_mature > 0 ? sum_p_mature / n_mature : 0;
  mean_p_immature = n_immature > 0 ? sum_p_immature / n_immature : 0;
  prop_correct_rep = n_correct * 1.0 / N;
}
