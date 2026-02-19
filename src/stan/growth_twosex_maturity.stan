// growth_twosex_maturity.stan
// Two-sex growth model with maturity-based parameterization
// k is derived from Lmat and tmat for biological interpretability
// Linf uses delta-gamma parameterization: Linf = Lmax + delta, delta ~ Gamma
//   This eliminates boundary pile-up when the posterior favors Linf near Lmax
// Pooling on Linf uses soft pairwise shrinkage on the reconstructed log scale,
//   preserving proportional interpretation of tau_Linf across species sizes
// Supports selective pooling: pool_maturity controls whether Lmat/tmat are pooled
// Non-centered parameterization for L0, Lmat, tmat (when pooled)
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
  int<lower=0, upper=1> pool_maturity; // 0=direct priors for Lmat/tmat, 1=also pool Lmat/tmat

  // Linf: delta-gamma parameterization
  // delta_Linf = Linf - Lmax ~ Gamma(alpha_delta, beta_delta)
  // alpha_delta = 1/CV_Linf^2, beta_delta = 1/(mu_delta * CV_Linf^2)
  // where mu_delta = Lmax * (Linf_multiplier - 1)
  vector<lower=0>[2] Lmax;           // Observed max length per sex
  vector<lower=0>[2] alpha_delta;    // Gamma shape per sex
  vector<lower=0>[2] beta_delta;     // Gamma rate per sex

  // Priors for L0 (log scale)
  vector[2] prior_L0_mu;
  vector<lower=0>[2] prior_L0_sigma;

  // Priors for Lmat (log scale)
  vector[2] prior_Lmat_mu;
  vector<lower=0>[2] prior_Lmat_sigma;

  // Priors for tmat (log scale)
  vector[2] prior_tmat_mu;
  vector<lower=0>[2] prior_tmat_sigma;

  // Hyperprior scales for partial pooling (half-normal scale)
  real<lower=0> prior_tau_Linf;      // On log(Linf) scale — proportional interpretation
  real<lower=0> prior_tau_L0;
  real<lower=0> prior_tau_Lmat;
  real<lower=0> prior_tau_tmat;

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
  // Linf: delta parameterization (no hierarchical mean needed)
  vector<lower=0>[2] delta_Linf;     // Excess above Lmax per sex
  real<lower=0> tau_Linf;            // Between-sex SD on log(Linf) for soft pooling

  // L0: standard non-centered parameterization
  real mu_L0;
  real<lower=0> tau_L0;
  vector[2] raw_L0;

  // Lmat/tmat: non-centered when pooled, direct when not
  real mu_Lmat;                      // Only used when pool_maturity == 1
  real<lower=0> tau_Lmat;            // Only used when pool_maturity == 1
  vector[2] raw_Lmat;

  real mu_tmat;                      // Only used when pool_maturity == 1
  real<lower=0> tau_tmat;            // Only used when pool_maturity == 1
  vector[2] raw_tmat;

  // Sex-specific observation error
  vector<lower=0>[2] sigma;

  // Degrees of freedom (if robust)
  array[robust] real<lower=2> nu_raw;
}

transformed parameters {
  // --- Linf: reconstruct from delta ---
  vector<lower=0>[2] Linf;
  for (s in 1:2)
    Linf[s] = Lmax[s] + delta_Linf[s];

  // --- L0, Lmat, tmat: construct log-scale parameters ---
  vector[2] log_L0;
  vector[2] log_Lmat;
  vector[2] log_tmat;

  if (use_pooling == 1) {
    // L0 always uses hierarchical structure when pooling is on
    for (s in 1:2) {
      log_L0[s] = mu_L0 + tau_L0 * raw_L0[s];
    }

    // Lmat and tmat: hierarchical only if pool_maturity == 1
    if (pool_maturity == 1) {
      for (s in 1:2) {
        log_Lmat[s] = mu_Lmat + tau_Lmat * raw_Lmat[s];
        log_tmat[s] = mu_tmat + tau_tmat * raw_tmat[s];
      }
    } else {
      // Direct parameterization for Lmat/tmat (no pooling)
      log_Lmat = raw_Lmat;
      log_tmat = raw_tmat;
    }
  } else {
    // No pooling: raw values are direct log parameters
    log_L0 = raw_L0;
    log_Lmat = raw_Lmat;
    log_tmat = raw_tmat;
  }

  // Transform to natural scale
  vector<lower=0>[2] L0;
  vector<lower=0>[2] Lmat;
  vector<lower=0>[2] tmat;

  for (s in 1:2) {
    L0[s] = exp(log_L0[s]);
    Lmat[s] = exp(log_Lmat[s]);
    tmat[s] = exp(log_tmat[s]);
  }

  // Constrain L0 < Lmat < Linf for biological plausibility
  vector<lower=0>[2] L0_constrained;
  vector<lower=0>[2] Lmat_constrained;

  for (s in 1:2) {
    L0_constrained[s] = fmin(L0[s], Lmat[s] * 0.95);
    Lmat_constrained[s] = fmin(fmax(Lmat[s], L0_constrained[s] * 1.05), Linf[s] * 0.99);
  }

  // Derive k from maturity parameters
  vector<lower=0>[2] k;

  for (s in 1:2) {
    real Linf_minus_L0 = fmax(Linf[s] - L0_constrained[s], 1e-6);
    real Linf_minus_Lmat = fmax(Linf[s] - Lmat_constrained[s], 1e-6);

    if (which_model == 1) {
      // von Bertalanffy: k = (1/tmat) * ln[(Linf - L0) / (Linf - Lmat)]
      real ratio = fmax(Linf_minus_L0 / Linf_minus_Lmat, 1.001);
      k[s] = (1.0 / tmat[s]) * log(ratio);
    } else if (which_model == 2) {
      // Gompertz: k = (1/tmat) * ln[ln(Linf/L0) / ln(Linf/Lmat)]
      real r0 = log(fmax(Linf[s] / L0_constrained[s], 1.001));
      real rmt = log(fmax(Lmat_constrained[s] / Linf[s], 1e-6));
      real ratio = fmax(-r0 / rmt, 1.001);
      k[s] = (1.0 / tmat[s]) * log(ratio);
    } else {
      // Logistic: k = (1/tmat) * ln[Lmat*(Linf-L0) / (L0*(Linf-Lmat))]
      real num = Lmat_constrained[s] * Linf_minus_L0;
      real den = fmax(L0_constrained[s] * Linf_minus_Lmat, 1e-6);
      real ratio = fmax(num / den, 1.001);
      k[s] = (1.0 / tmat[s]) * log(ratio);
    }
  }

  // Degrees of freedom
  real nu;
  if (robust == 1) nu = nu_raw[1];
  else nu = 100.0;

  // Predicted lengths
  vector<lower=0>[N] mu;

  for (i in 1:N) {
    int s = sex[i];
    real a = age[i];

    if (which_model == 1) {
      mu[i] = Linf[s] - (Linf[s] - L0_constrained[s]) * exp(-k[s] * a);
    } else if (which_model == 2) {
      real r0 = log(Linf[s] / L0_constrained[s]);
      mu[i] = Linf[s] * exp(-r0 * exp(-k[s] * a));
    } else {
      real A = exp(-k[s] * a);
      mu[i] = Linf[s] / (1 + A * (Linf[s] / L0_constrained[s] - 1));
    }
  }
}

model {
  // === Linf priors: always independent gamma on delta ===
  delta_Linf ~ gamma(alpha_delta, beta_delta);

  // === Pooling structure ===
  if (use_pooling == 1) {
    // Soft pairwise shrinkage on log(Linf) — proportional interpretation
    // log(Linf[1]) - log(Linf[2]) ~ Normal(0, tau_Linf)
    // Equivalent to: Linf[1]/Linf[2] ~ LogNormal(0, tau_Linf)
    log(Linf[1]) - log(Linf[2]) ~ normal(0, tau_Linf);
    tau_Linf ~ normal(0, prior_tau_Linf);

    // L0: standard non-centered hierarchical
    mu_L0 ~ normal(mean(prior_L0_mu), mean(prior_L0_sigma));
    tau_L0 ~ normal(0, prior_tau_L0);
    raw_L0 ~ std_normal();

    if (pool_maturity == 1) {
      // Full pooling: Lmat/tmat also get hierarchical structure
      // Widened hyperpriors (3x SD) to prevent over-constraint
      mu_Lmat ~ normal(mean(prior_Lmat_mu), 3 * mean(prior_Lmat_sigma));
      mu_tmat ~ normal(mean(prior_tmat_mu), 3 * mean(prior_tmat_sigma));

      tau_Lmat ~ normal(0, prior_tau_Lmat);
      tau_tmat ~ normal(0, prior_tau_tmat);

      raw_Lmat ~ std_normal();
      raw_tmat ~ std_normal();

      // Anchoring priors (3x widened) to prevent collapse
      for (s in 1:2) {
        target += normal_lpdf(log_Lmat[s] | prior_Lmat_mu[s], 3 * prior_Lmat_sigma[s]);
        target += normal_lpdf(log_tmat[s] | prior_tmat_mu[s], 3 * prior_tmat_sigma[s]);
      }
    } else {
      // Selective pooling: Lmat/tmat use direct sex-specific priors only
      for (s in 1:2) {
        raw_Lmat[s] ~ normal(prior_Lmat_mu[s], prior_Lmat_sigma[s]);
        raw_tmat[s] ~ normal(prior_tmat_mu[s], prior_tmat_sigma[s]);
      }

      // Weakly constrain unused hyperparams
      mu_Lmat ~ normal(0, 1);
      mu_tmat ~ normal(0, 1);
      tau_Lmat ~ normal(0, 0.1);
      tau_tmat ~ normal(0, 0.1);
    }

  } else {
    // === No pooling: independent priors per sex ===
    for (s in 1:2) {
      raw_L0[s] ~ normal(prior_L0_mu[s], prior_L0_sigma[s]);
      raw_Lmat[s] ~ normal(prior_Lmat_mu[s], prior_Lmat_sigma[s]);
      raw_tmat[s] ~ normal(prior_tmat_mu[s], prior_tmat_sigma[s]);
    }

    // Weakly constrain unused hierarchical params
    mu_L0 ~ normal(0, 1);
    mu_Lmat ~ normal(0, 1);
    mu_tmat ~ normal(0, 1);
    tau_Linf ~ normal(0, 0.1);
    tau_L0 ~ normal(0, 0.1);
    tau_Lmat ~ normal(0, 0.1);
    tau_tmat ~ normal(0, 0.1);
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

  // Sex differences (Female - Male)
  real Linf_diff = Linf[1] - Linf[2];
  real L0_diff = L0[1] - L0[2];
  real Lmat_diff = Lmat[1] - Lmat[2];
  real tmat_diff = tmat[1] - tmat[2];
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
