// growth_single_k.stan
// Single-sex growth model with traditional k parameterization
// Supports von Bertalanffy (1), Gompertz (2), and Logistic (3) growth curves
// Observation model: lognormal (robust=0) or Student-t (robust=1)

data {
  int<lower=1> N;                    // Number of observations
  vector<lower=0>[N] length;         // Observed lengths
  vector<lower=0>[N] age;            // Observed ages
  array[N] int<lower=1> row_id;      // Original row indices for tracking
  
  int<lower=1, upper=3> which_model; // 1=VBGM, 2=Gompertz, 3=Logistic
  int<lower=0, upper=1> robust;      // 0=lognormal, 1=Student-t errors
  
  // Priors for Linf (log scale)
  real prior_Linf_mu;
  real<lower=0> prior_Linf_sigma;
  real<lower=0> Linf_lower;          // Hard lower bound (observed max)
  
  // Priors for L0 (log scale)
  real prior_L0_mu;
  real<lower=0> prior_L0_sigma;
  
  // Priors for k (log scale)
  real prior_k_mu;
  real<lower=0> prior_k_sigma;
  
  // Priors for sigma (observation error)
  real loc_sig;
  real<lower=0> scale_sig;
}

transformed data {
  vector[N] log_length = log(length);
}

parameters {
  real<lower=log(Linf_lower)> log_Linf;  // Log asymptotic length
  real log_L0;                            // Log length at birth
  real log_k;                             // Log growth rate
  real<lower=0> sigma;                    // Observation error (log scale)
  array[robust] real<lower=2> nu_raw;    // Degrees of freedom (Student-t only)
}

transformed parameters {
  real<lower=Linf_lower> Linf = exp(log_Linf);
  real<lower=0> L0 = exp(log_L0);
  real<lower=0> k = exp(log_k);
  
  // Constrain L0 < Linf
  real<lower=0> L0_constrained = fmin(L0, Linf * 0.99);
  
  // Degrees of freedom for robust model
  real nu = robust ? nu_raw[1] : 100.0;
  
  // Predicted mean lengths
  vector<lower=0>[N] mu;
  
  for (i in 1:N) {
    real a = age[i];
    
    if (which_model == 1) {
      // von Bertalanffy
      mu[i] = Linf - (Linf - L0_constrained) * exp(-k * a);
    } else if (which_model == 2) {
      // Gompertz
      real r0 = log(Linf / L0_constrained);
      mu[i] = Linf * exp(-r0 * exp(-k * a));
    } else {
      // Logistic
      real A = exp(-k * a);
      mu[i] = Linf / (1 + A * (Linf / L0_constrained - 1));
    }
  }
}

model {
  // Priors
  log_Linf ~ normal(prior_Linf_mu, prior_Linf_sigma);
  log_L0 ~ normal(prior_L0_mu, prior_L0_sigma);
  log_k ~ normal(prior_k_mu, prior_k_sigma);
  sigma ~ cauchy(loc_sig, scale_sig);
  
  if (robust) {
    nu_raw[1] ~ gamma(2, 0.1);  // Prior for df, peaked around 20
  }
  
  // Likelihood
  if (robust) {
    // Student-t on log scale
    for (i in 1:N) {
      log_length[i] ~ student_t(nu, log(mu[i]), sigma);
    }
  } else {
    // Lognormal
    log_length ~ normal(log(mu), sigma);
  }
}

generated quantities {
  // Log-likelihood for LOO
  vector[N] log_lik;
  
  // Posterior predictive
  vector[N] y_pred;
  vector[N] y_rep;
  vector[N] residual;
  vector[N] std_residual;
  
  // Summary statistics
  real mean_residual;
  real sd_residual;
  real rmse;
  real mae;
  int<lower=0, upper=N> n_in_CI;
  
  // Derived quantity: k already available as parameter
  
  for (i in 1:N) {
    // Log-likelihood
    if (robust) {
      log_lik[i] = student_t_lpdf(log_length[i] | nu, log(mu[i]), sigma);
    } else {
      log_lik[i] = normal_lpdf(log_length[i] | log(mu[i]), sigma);
    }
    
    // Posterior predictive (mean prediction)
    y_pred[i] = mu[i];
    
    // Posterior predictive replications
    if (robust) {
      real log_rep = student_t_rng(nu, log(mu[i]), sigma);
      y_rep[i] = exp(log_rep);
    } else {
      y_rep[i] = lognormal_rng(log(mu[i]), sigma);
    }
    
    // Residuals
    residual[i] = length[i] - mu[i];
    std_residual[i] = (log_length[i] - log(mu[i])) / sigma;
  }
  
  // Summary statistics
  mean_residual = mean(residual);
  sd_residual = sd(residual);
  rmse = sqrt(mean(square(residual)));
  mae = mean(abs(residual));
  
  // Count observations within 95% CI
  {
    int count = 0;
    for (i in 1:N) {
      real lower_ci = exp(log(mu[i]) - 1.96 * sigma);
      real upper_ci = exp(log(mu[i]) + 1.96 * sigma);
      if (length[i] >= lower_ci && length[i] <= upper_ci) {
        count += 1;
      }
    }
    n_in_CI = count;
  }
}
