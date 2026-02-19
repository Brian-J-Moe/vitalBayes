// growth_single_maturity.stan
// Single-sex growth model with maturity-based parameterization
// k is derived from Lmat and tmat, providing biological interpretability
// Supports von Bertalanffy (1), Gompertz (2), and Logistic (3) growth curves

data {
  int<lower=1> N;                    // Number of observations
  vector<lower=0>[N] length;         // Observed lengths
  vector<lower=0>[N] age;            // Observed ages
  array[N] int<lower=1> row_id;      // Original row indices for tracking

  int<lower=1, upper=3> which_model; // 1=VBGM, 2=Gompertz, 3=Logistic
  int<lower=0, upper=1> robust;      // 0=lognormal, 1=Student-t errors

  // Priors for Linf
  real Lmax;
  real<lower=0> alpha_delta;
  real<lower=0> beta_delta;

  // Priors for L0 (log scale)
  real prior_L0_mu;
  real<lower=0> prior_L0_sigma;

  // Priors for Lmat (log scale)
  real prior_Lmat_mu;
  real<lower=0> prior_Lmat_sigma;

  // Priors for tmat (log scale)
  real prior_tmat_mu;
  real<lower=0> prior_tmat_sigma;

  // Priors for sigma (observation error)
  real loc_sig;
  real<lower=0> scale_sig;
}

transformed data {
  vector[N] log_length = log(length);
}

parameters {
  real<lower=0> delta_Linf;               // Length difference above Lmax
  real log_L0;                            // Log length at birth
  real log_Lmat;                          // Log length at maturity
  real log_tmat;                          // Log age at maturity
  real<lower=0> sigma;                    // Observation error (log scale)
  array[robust] real<lower=2> nu_raw;    // Degrees of freedom (Student-t only)
}

transformed parameters {
  real Linf = Lmax + delta_Linf;
  real<lower=0> L0 = exp(log_L0);
  real<lower=0> Lmat = exp(log_Lmat);
  real<lower=0> tmat = exp(log_tmat);

  // Constrain L0 < Lmat < Linf for biological plausibility
  real L0_constrained = fmin(L0, Lmat * 0.95);
  real Lmat_constrained = fmin(fmax(Lmat, L0_constrained * 1.05), Linf * 0.99);

  // Derive k from maturity parameters with safeguards
  real<lower=0> k;

  // Degrees of freedom for robust model
  real nu = robust ? nu_raw[1] : 100.0;

  // Predicted mean lengths
  vector<lower=0>[N] mu;

  // Ensure safe denominators and ratios for log()
  real Linf_minus_L0 = fmax(Linf - L0_constrained, 1e-6);
  real Linf_minus_Lmat = fmax(Linf - Lmat_constrained, 1e-6);

  // Compute k based on growth model
  if (which_model == 1) {
    // von Bertalanffy: k = (1/tmat) * ln((Linf - L0)/(Linf - Lmat))
    real ratio = fmax(Linf_minus_L0 / Linf_minus_Lmat, 1.001);
    k = (1.0 / tmat) * log(ratio);
  } else if (which_model == 2) {
    // Gompertz: k = (1/tmat) * ln[ln(Linf/L0) / ln(Linf/Lmat)]
    real r0 = log(fmax(Linf / L0_constrained, 1.001));
    real rmt = log(fmax(Lmat_constrained / Linf, 1e-6));
    // rmt is negative, so -r0/rmt is positive
    real ratio = fmax(-r0 / rmt, 1.001);
    k = (1.0 / tmat) * log(ratio);
  } else {
    // Logistic: k = (1/tmat) * ln((Lmat*(Linf-L0)) / (L0*(Linf-Lmat)))
    real num = Lmat_constrained * Linf_minus_L0;
    real den = fmax(L0_constrained * Linf_minus_Lmat, 1e-6);
    real ratio = fmax(num / den, 1.001);
    k = (1.0 / tmat) * log(ratio);
  }

  // Compute predicted lengths
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
  delta_Linf ~ gamma(alpha_delta, beta_delta);
  log_L0 ~ normal(prior_L0_mu, prior_L0_sigma);
  log_Lmat ~ normal(prior_Lmat_mu, prior_Lmat_sigma);
  log_tmat ~ normal(prior_tmat_mu, prior_tmat_sigma);
  sigma ~ cauchy(loc_sig, scale_sig);

  if (robust) {
    nu_raw[1] ~ gamma(2, 0.1);
  }

  // Likelihood
  if (robust) {
    for (i in 1:N) {
      log_length[i] ~ student_t(nu, log(mu[i]), sigma);
    }
  } else {
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

  for (i in 1:N) {
    // Log-likelihood
    if (robust) {
      log_lik[i] = student_t_lpdf(log_length[i] | nu, log(mu[i]), sigma);
    } else {
      log_lik[i] = normal_lpdf(log_length[i] | log(mu[i]), sigma);
    }

    // Posterior predictive
    y_pred[i] = mu[i];

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
