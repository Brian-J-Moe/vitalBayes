// growth_twosex_k.stan
// Two-sex growth model with traditional k parameterization
// Supports optional partial pooling via hyperparameters (use_pooling flag)
// Non-centered parameterization for numerical stability
// Supports von Bertalanffy (1), Gompertz (2), and Logistic (3) growth curves

data {
  int<lower=1> N;                    // Number of observations
  vector<lower=0>[N] length;         // Observed lengths
  vector<lower=0>[N] age;            // Observed ages
  array[N] int<lower=1, upper=2> sex;// Sex indicator (1=female, 2=male)
  array[N] int<lower=1> row_id;      // Original row indices
  
  int<lower=1, upper=3> which_model; // 1=VBGM, 2=Gompertz, 3=Logistic
  int<lower=0, upper=1> robust;      // 0=lognormal, 1=Student-t errors
  int<lower=0, upper=1> use_pooling; // 0=independent, 1=partial pooling
  
  // Priors for Linf (log scale) - vectors of length 2 for [female, male]
  vector[2] prior_Linf_mu;
  vector<lower=0>[2] prior_Linf_sigma;
  vector<lower=0>[2] Linf_lower;
  
  // Priors for L0 (log scale)
  vector[2] prior_L0_mu;
  vector<lower=0>[2] prior_L0_sigma;
  
  // Priors for k (log scale)
  vector[2] prior_k_mu;
  vector<lower=0>[2] prior_k_sigma;
  
  // Hyperprior scales for partial pooling (half-normal scale)
  real<lower=0> prior_tau_Linf;
  real<lower=0> prior_tau_L0;
  real<lower=0> prior_tau_k;
  
  // Priors for sigma
  real loc_sig;
  real<lower=0> scale_sig;
}

transformed data {
  vector[N] log_length = log(length);
  
  // Count observations per sex
  int n_female = 0;
  int n_male = 0;
  for (i in 1:N) {
    if (sex[i] == 1) n_female += 1;
    else n_male += 1;
  }
}

parameters {
  // Population means (for pooling)
  real mu_Linf;
  real mu_L0;
  real mu_k;
  
  // Between-sex standard deviations (for pooling)
  real<lower=0> tau_Linf;
  real<lower=0> tau_L0;
  real<lower=0> tau_k;
  
  // Raw deviates (always present)
  vector[2] raw_Linf;
  vector[2] raw_L0;
  vector[2] raw_k;
  
  // Sex-specific observation error
  vector<lower=0>[2] sigma;
  
  // Degrees of freedom (if robust)
  array[robust] real<lower=2> nu_raw;
}

transformed parameters {
  // Construct log-scale parameters
  vector[2] log_Linf;
  vector[2] log_L0;
  vector[2] log_k;
  
  if (use_pooling == 1) {
    // Non-centered: log_param = mu + tau * raw
    for (s in 1:2) {
      log_Linf[s] = mu_Linf + tau_Linf * raw_Linf[s];
      log_L0[s] = mu_L0 + tau_L0 * raw_L0[s];
      log_k[s] = mu_k + tau_k * raw_k[s];
    }
  } else {
    // No pooling: raw values are direct log parameters
    log_Linf = raw_Linf;
    log_L0 = raw_L0;
    log_k = raw_k;
  }
  
  // Transform to natural scale with constraints
  vector<lower=0>[2] Linf;
  vector<lower=0>[2] L0;
  vector<lower=0>[2] k;
  
  for (s in 1:2) {
    Linf[s] = fmax(exp(log_Linf[s]), Linf_lower[s]);
    L0[s] = exp(log_L0[s]);
    k[s] = exp(log_k[s]);
  }
  
  // Constrain L0 < Linf
  vector<lower=0>[2] L0_constrained;
  for (s in 1:2) {
    L0_constrained[s] = fmin(L0[s], Linf[s] * 0.99);
  }
  
  // Degrees of freedom
  real nu = robust ? nu_raw[1] : 100.0;
  
  // Predicted lengths
  vector<lower=0>[N] mu;
  
  for (i in 1:N) {
    int s = sex[i];
    real a = age[i];
    
    if (which_model == 1) {
      // von Bertalanffy
      mu[i] = Linf[s] - (Linf[s] - L0_constrained[s]) * exp(-k[s] * a);
    } else if (which_model == 2) {
      // Gompertz
      real r0 = log(Linf[s] / L0_constrained[s]);
      mu[i] = Linf[s] * exp(-r0 * exp(-k[s] * a));
    } else {
      // Logistic
      real A = exp(-k[s] * a);
      mu[i] = Linf[s] / (1 + A * (Linf[s] / L0_constrained[s] - 1));
    }
  }
}

model {
  // Priors
  if (use_pooling == 1) {
    // Hyperpriors
    mu_Linf ~ normal(mean(prior_Linf_mu), 1);
    mu_L0 ~ normal(mean(prior_L0_mu), 1);
    mu_k ~ normal(mean(prior_k_mu), 1);
    
    // Half-normal priors on between-sex SD
    tau_Linf ~ normal(0, prior_tau_Linf);
    tau_L0 ~ normal(0, prior_tau_L0);
    tau_k ~ normal(0, prior_tau_k);
    
    // Standard normal for non-centered raw deviates
    raw_Linf ~ std_normal();
    raw_L0 ~ std_normal();
    raw_k ~ std_normal();
    
  } else {
    // Independent priors per sex (raw values are log parameters)
    for (s in 1:2) {
      raw_Linf[s] ~ normal(prior_Linf_mu[s], prior_Linf_sigma[s]);
      raw_L0[s] ~ normal(prior_L0_mu[s], prior_L0_sigma[s]);
      raw_k[s] ~ normal(prior_k_mu[s], prior_k_sigma[s]);
    }
    
    // Weakly constrain unused params
    mu_Linf ~ normal(0, 1);
    mu_L0 ~ normal(0, 1);
    mu_k ~ normal(0, 1);
    tau_Linf ~ normal(0, 1);
    tau_L0 ~ normal(0, 1);
    tau_k ~ normal(0, 1);
  }
  
  // Observation error priors
  sigma ~ cauchy(loc_sig, scale_sig);
  
  if (robust) {
    nu_raw[1] ~ gamma(2, 0.1);
  }
  
  // Likelihood
  for (i in 1:N) {
    int s = sex[i];
    if (robust) {
      log_length[i] ~ student_t(nu, log(mu[i]), sigma[s]);
    } else {
      log_length[i] ~ normal(log(mu[i]), sigma[s]);
    }
  }
}

generated quantities {
  // Log-likelihood
  vector[N] log_lik;
  
  // Posterior predictive
  vector[N] y_pred;
  vector[N] y_rep;
  vector[N] residual;
  vector[N] std_residual;
  
  // Sex differences (always computed)
  real Linf_diff = Linf[1] - Linf[2];  // Female - Male
  real L0_diff = L0[1] - L0[2];
  real k_diff = k[1] - k[2];
  
  // Sex-specific summary stats
  real mean_residual_f = 0;
  real mean_residual_m = 0;
  real rmse_f = 0;
  real rmse_m = 0;
  int n_in_CI_f = 0;
  int n_in_CI_m = 0;
  
  // Temporary accumulators
  real sum_resid_f = 0;
  real sum_resid_m = 0;
  real sum_sq_f = 0;
  real sum_sq_m = 0;
  
  for (i in 1:N) {
    int s = sex[i];
    
    // Log-likelihood
    if (robust) {
      log_lik[i] = student_t_lpdf(log_length[i] | nu, log(mu[i]), sigma[s]);
    } else {
      log_lik[i] = normal_lpdf(log_length[i] | log(mu[i]), sigma[s]);
    }
    
    // Predictions
    y_pred[i] = mu[i];
    
    if (robust) {
      real log_rep = student_t_rng(nu, log(mu[i]), sigma[s]);
      y_rep[i] = exp(log_rep);
    } else {
      y_rep[i] = lognormal_rng(log(mu[i]), sigma[s]);
    }
    
    // Residuals
    residual[i] = length[i] - mu[i];
    std_residual[i] = (log_length[i] - log(mu[i])) / sigma[s];
    
    // Accumulate by sex
    if (s == 1) {
      sum_resid_f += residual[i];
      sum_sq_f += square(residual[i]);
      
      real lower_ci = exp(log(mu[i]) - 1.96 * sigma[s]);
      real upper_ci = exp(log(mu[i]) + 1.96 * sigma[s]);
      if (length[i] >= lower_ci && length[i] <= upper_ci) {
        n_in_CI_f += 1;
      }
    } else {
      sum_resid_m += residual[i];
      sum_sq_m += square(residual[i]);
      
      real lower_ci = exp(log(mu[i]) - 1.96 * sigma[s]);
      real upper_ci = exp(log(mu[i]) + 1.96 * sigma[s]);
      if (length[i] >= lower_ci && length[i] <= upper_ci) {
        n_in_CI_m += 1;
      }
    }
  }
  
  // Finalize sex-specific stats
  if (n_female > 0) {
    mean_residual_f = sum_resid_f / n_female;
    rmse_f = sqrt(sum_sq_f / n_female);
  }
  if (n_male > 0) {
    mean_residual_m = sum_resid_m / n_male;
    rmse_m = sqrt(sum_sq_m / n_male);
  }
}
