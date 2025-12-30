// length_at_maturity_twosex.stan
// Two-sex binomial regression for length-at-50%-maturity (L50)
// Uses probit link with optional partial pooling between sexes
// Non-centered parameterization with half-normal priors on tau
// Probit chosen for threshold-crossing interpretation: latent developmental
// readiness is normally distributed, individual matures when readiness > 0

data {
  int<lower=1> N;                    // Number of observations
  vector<lower=0>[N] length;         // Observed lengths
  array[N] int<lower=0, upper=1> mature; // Maturity status
  array[N] int<lower=1, upper=2> sex;    // Sex indicator (1=female, 2=male)
  
  int<lower=0, upper=1> use_pooling; // 0=independent, 1=partial pooling
  
  // Priors for L50 (log scale) - vectors of length 2 for [female, male]
  vector[2] prior_L50_mu;
  vector<lower=0>[2] prior_L50_sigma;
  
  // Priors for slope (log scale)
  real prior_slope_mu;
  real<lower=0> prior_slope_sigma;
  
  // Hyperprior scale for partial pooling (half-normal scale)
  real<lower=0> prior_tau_L50;
}

transformed data {
  // Count observations per sex
  int n_female = 0;
  int n_male = 0;
  for (i in 1:N) {
    if (sex[i] == 1) n_female += 1;
    else n_male += 1;
  }
}

parameters {
  // For pooling: population mean and between-sex SD
  real mu_L50;
  real<lower=0> tau_L50;
  
  // Raw deviates for non-centered parameterization (always present)
  vector[2] raw_L50;
  
  // Slope parameters (not pooled - often genuinely differ between sexes)
  vector[2] log_slope;
}

transformed parameters {
  // Construct log-scale L50
  vector[2] log_L50;
  
  if (use_pooling == 1) {
    // Non-centered: log_L50[s] = mu + tau * raw[s]
    for (s in 1:2) {
      log_L50[s] = mu_L50 + tau_L50 * raw_L50[s];
    }
  } else {
    // No pooling: use raw values directly as log(L50)
    // raw_L50 serves as log_L50_direct in this case
    log_L50 = raw_L50;
  }
  
  // Transform to natural scale
  vector<lower=0>[2] L50 = exp(log_L50);
  vector<lower=0>[2] slope = exp(log_slope);
  
  // Probabilities (probit link)
  vector<lower=0, upper=1>[N] p;
  
  for (i in 1:N) {
    int s = sex[i];
    real eta = slope[s] * (length[i] - L50[s]);
    p[i] = Phi(eta);
  }
}

model {
  // Priors
  if (use_pooling == 1) {
    // Hyperprior for population mean
    mu_L50 ~ normal(mean(prior_L50_mu), mean(prior_L50_sigma));
    
    // Half-normal prior on between-sex SD
    tau_L50 ~ normal(0, prior_tau_L50);
    
    // Standard normal for non-centered deviates
    raw_L50 ~ std_normal();
  } else {
    // Independent priors per sex (raw_L50 is log_L50 directly)
    for (s in 1:2) {
      raw_L50[s] ~ normal(prior_L50_mu[s], prior_L50_sigma[s]);
    }
    
    // Weakly constrain unused params to avoid sampling issues
    mu_L50 ~ normal(0, 1);
    tau_L50 ~ normal(0, 1);
  }
  
  // Slope priors (same for both sexes but estimated separately)
  log_slope ~ normal(prior_slope_mu, prior_slope_sigma);
  
  // Likelihood
  mature ~ bernoulli(p);
}

generated quantities {
  // Log-likelihood
  vector[N] log_lik;
  
  // Posterior predictive
  vector[N] p_pred = p;
  array[N] int mature_rep;
  
  // Derived quantities per sex (probit scale: Phi^{-1}(0.05) = -1.645)
  vector[2] L05;
  vector[2] L95;
  vector[2] transition_width;
  
  for (s in 1:2) {
    L05[s] = L50[s] - 1.6449 / slope[s];
    L95[s] = L50[s] + 1.6449 / slope[s];
    transition_width[s] = L95[s] - L05[s];
  }
  
  // Sex difference
  real L50_diff = L50[1] - L50[2];  // Female - Male
  
  // Summary statistics by sex
  real mean_p_mature_f = 0;
  real mean_p_mature_m = 0;
  real mean_p_immature_f = 0;
  real mean_p_immature_m = 0;
  
  int n_mature_f = 0;
  int n_mature_m = 0;
  int n_immature_f = 0;
  int n_immature_m = 0;
  
  real sum_p_mature_f = 0;
  real sum_p_mature_m = 0;
  real sum_p_immature_f = 0;
  real sum_p_immature_m = 0;
  
  int n_correct = 0;
  real prop_correct_rep;
  
  for (i in 1:N) {
    int s = sex[i];
    
    log_lik[i] = bernoulli_lpmf(mature[i] | p[i]);
    mature_rep[i] = bernoulli_rng(p[i]);
    
    if (mature_rep[i] == mature[i]) {
      n_correct += 1;
    }
    
    // Accumulate by sex
    if (s == 1) {
      if (mature[i] == 1) {
        n_mature_f += 1;
        sum_p_mature_f += p[i];
      } else {
        n_immature_f += 1;
        sum_p_immature_f += p[i];
      }
    } else {
      if (mature[i] == 1) {
        n_mature_m += 1;
        sum_p_mature_m += p[i];
      } else {
        n_immature_m += 1;
        sum_p_immature_m += p[i];
      }
    }
  }
  
  mean_p_mature_f = n_mature_f > 0 ? sum_p_mature_f / n_mature_f : 0;
  mean_p_mature_m = n_mature_m > 0 ? sum_p_mature_m / n_mature_m : 0;
  mean_p_immature_f = n_immature_f > 0 ? sum_p_immature_f / n_immature_f : 0;
  mean_p_immature_m = n_immature_m > 0 ? sum_p_immature_m / n_immature_m : 0;
  prop_correct_rep = n_correct * 1.0 / N;
}
