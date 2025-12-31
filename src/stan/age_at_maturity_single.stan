// age_at_maturity_single.stan
// Single-sex binomial regression for age-at-50%-maturity (t50)
// Uses probit link function with lognormal priors on t50
// Probit chosen for threshold-crossing interpretation: latent developmental
// readiness is normally distributed, individual matures when readiness > 0

data {
  int<lower=1> N;                    // Number of observations
  vector<lower=0>[N] age;            // Observed ages
  array[N] int<lower=0, upper=1> mature; // Maturity status (0=immature, 1=mature)
  
  // Priors for t50 (log scale)
  real prior_t50_mu;
  real<lower=0> prior_t50_sigma;
  
  // Priors for slope (log scale for positivity)
  real prior_slope_mu;
  real<lower=0> prior_slope_sigma;
}

parameters {
  real log_t50;                      // Log age at 50% maturity
  real log_slope;                    // Log slope (transition steepness)
}

transformed parameters {
  real<lower=0> t50 = exp(log_t50);
  real<lower=0> slope = exp(log_slope);
  
  // Linear predictor with clamping for numerical stability
  vector[N] eta_raw = slope * (age - t50);
  vector[N] eta;
  for (i in 1:N) {
    eta[i] = fmin(fmax(eta_raw[i], -20.0), 20.0);
  }
  
  // Probability of maturity (probit link)
  vector<lower=0, upper=1>[N] p = Phi(eta);
}

model {
  // Priors
  log_t50 ~ normal(prior_t50_mu, prior_t50_sigma);
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
  
  // Derived quantities (probit scale: Phi^{-1}(0.05) = -1.645, Phi^{-1}(0.95) = 1.645)
  real<lower=0> t05 = t50 - 1.6449 / slope;   // Age at 5% maturity
  real<lower=0> t95 = t50 + 1.6449 / slope;   // Age at 95% maturity
  real<lower=0> transition_width = t95 - t05;  // Width of transition zone
  
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
    
    if (mature_rep[i] == mature[i]) {
      n_correct += 1;
    }
    
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
