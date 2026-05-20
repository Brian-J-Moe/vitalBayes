#' Simulate artificial von Bertalanffy growth data with truncated-exponential ages
#' and sex-specific logistic maturity ogives
#'
#' Simulates individual-level VB parameters (L0, Lmat, tmat, Linf) by sex, draws ages
#' from a truncated exponential (older ages increasingly rare), computes VB mean length,
#' adds observation error to obtain fl, then samples maturity from a logistic ogive on
#' observed fl with sex-specific steepness (w95).
#'
#' VB formulation:
#'   k  = (1/tmat) * log((Linf - L0) / (Linf - Lmat))
#'   mu = Linf - (Linf - L0) * exp(-k * age)
#'
#' @param n_female Integer. Number of female (non-embryo) observations.
#' @param n_male Integer. Number of male (non-embryo) observations.
#' @param n_embryo Integer. Number of embryo observations.
#' @param age_max_female Numeric. Maximum age (years) for females.
#' @param age_max_male Numeric. Maximum age (years) for males.
#' @param age_tail_frac Numeric in (0,1). Reference age = age_tail_frac * age_max for tuning the tail.
#' @param age_tail_prob Numeric in (0,1). Approx target P(age > age_tail_frac*age_max) for the underlying exponential.
#' @param sigma_fl Numeric. SD of observation error (cm) added to VB mean length.
#' @param mat_w95_female Numeric. Female ogive width (cm) between 5% and 95% maturity.
#' @param mat_w95_male Numeric. Male ogive width (cm) between 5% and 95% maturity.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#' @param round_digits Integer or NULL. If not NULL, round fl and age to this many digits.
#'
#' @return data.table with columns: sex, mat, fl, age, embryo
#' @export
simulate_vb_growth_data <- function(n_female = 150L,
                                    n_male   = 120L,
                                    n_embryo = 26L,
                                    age_max_female = 35,
                                    age_max_male   = 31,
                                    age_tail_frac  = 0.70,
                                    age_tail_prob  = 0.30,
                                    sigma_fl = 12,
                                    mat_w95_female = 20,
                                    mat_w95_male   = 18,
                                    seed = 123,
                                    round_digits = 1L) {

  if (!is.null(seed)) set.seed(seed)

  # ---- internal helpers ----

  #' @keywords internal
  .rtnorm <- function(n, mean, sd, lower = -Inf, upper = Inf, max_iter = 200L) {
    x <- stats::rnorm(n, mean, sd)
    bad <- (!is.finite(x)) | (x < lower) | (x > upper)
    it <- 0L
    while (any(bad) && it < max_iter) {
      x[bad] <- stats::rnorm(sum(bad), mean, sd)
      bad <- (!is.finite(x)) | (x < lower) | (x > upper)
      it <- it + 1L
    }
    x <- pmax(x, lower)
    x <- pmin(x, upper)
    x
  }

  #' Draw from Exp(rate) truncated to [0, max_age] via inverse CDF
  #' @keywords internal
  .rtruncexp <- function(n, rate, max_age) {
    umax <- 1 - exp(-rate * max_age)
    u <- stats::runif(n, min = 0, max = umax)
    -log1p(-u) / rate
  }

  #' Convert desired rarity of old ages into exponential rate
  #' @keywords internal
  .rate_from_tail <- function(age_max, tail_frac, tail_prob) {
    a <- tail_frac * age_max
    -log(tail_prob) / a
  }

  #' Enforce Linf > Lmat > L0 with small gaps; resample offenders from sex-specific priors
  #' @keywords internal
  .enforce_constraints <- function(L0, Lmat, tmat, Linf,
                                   Lmat_mean, Lmat_sd,
                                   Linf_mean, Linf_sd,
                                   gap_L0_Lmat = 20,
                                   gap_Lmat_Linf = 50,
                                   max_iter = 500L) {

    it <- 0L
    bad <- (Lmat <= (L0 + gap_L0_Lmat)) | (Linf <= (Lmat + gap_Lmat_Linf)) |
      (!is.finite(L0 + Lmat + tmat + Linf))

    while (any(bad) && it < max_iter) {
      idx <- which(bad)
      Lmat[idx] <- .rtnorm(length(idx), mean = Lmat_mean, sd = Lmat_sd, lower = 1)
      Linf[idx] <- .rtnorm(length(idx), mean = Linf_mean, sd = Linf_sd, lower = 1)

      bad <- (Lmat <= (L0 + gap_L0_Lmat)) | (Linf <= (Lmat + gap_Lmat_Linf)) |
        (!is.finite(L0 + Lmat + tmat + Linf))
      it <- it + 1L
    }

    bad <- (Lmat <= (L0 + gap_L0_Lmat)) | (Linf <= (Lmat + gap_Lmat_Linf))
    if (any(bad)) {
      Lmat[bad] <- L0[bad] + gap_L0_Lmat + abs(Lmat[bad] - (L0[bad] + gap_L0_Lmat))
      Linf[bad] <- Lmat[bad] + gap_Lmat_Linf + abs(Linf[bad] - (Lmat[bad] + gap_Lmat_Linf))
    }

    list(L0 = L0, Lmat = Lmat, tmat = tmat, Linf = Linf)
  }

  #' @keywords internal
  .vb_mu <- function(age, L0, Lmat, tmat, Linf) {
    k <- (1.0 / tmat) * log((Linf - L0) / (Linf - Lmat))
    Linf - (Linf - L0) * exp(-k * age)
  }

  #' Sample maturity from logistic ogive on observed length
  #' w95 controls steepness: smaller => steeper.
  #' @keywords internal
  .sample_mat_logistic <- function(fl, Lmat, w95) {
    if (!is.finite(w95) || w95 <= 0) stop("w95 must be > 0")
    beta <- (2 * log(19)) / w95
    p <- stats::plogis(beta * (fl - Lmat))
    stats::rbinom(length(p), size = 1L, prob = p)
  }

  #' Simulate one sex group
  #' @keywords internal
  .simulate_sex <- function(sex, n,
                            L0_mean,   L0_sd,
                            Lmat_mean, Lmat_sd,
                            tmat_mean, tmat_sd,
                            Linf_mean, Linf_sd,
                            age_max,   mat_w95) {

    L0   <- .rtnorm(n, mean = L0_mean,   sd = L0_sd, lower = 1)
    Lmat <- .rtnorm(n, mean = Lmat_mean, sd = Lmat_sd, lower = 1)
    tmat <- .rtnorm(n, mean = tmat_mean, sd = tmat_sd, lower = 0.1)
    Linf <- .rtnorm(n, mean = Linf_mean, sd = Linf_sd, lower = 1)

    pars <- .enforce_constraints(
      L0 = L0, Lmat = Lmat, tmat = tmat, Linf = Linf,
      Lmat_mean = Lmat_mean, Lmat_sd = Lmat_sd,
      Linf_mean = Linf_mean, Linf_sd = Linf_sd
    )

    rate <- .rate_from_tail(age_max = age_max, tail_frac = age_tail_frac, tail_prob = age_tail_prob)
    age  <- .rtruncexp(n = n, rate = rate, max_age = age_max)

    mu <- .vb_mu(age = age, L0 = pars$L0, Lmat = pars$Lmat, tmat = pars$tmat, Linf = pars$Linf)
    fl <- stats::rnorm(n, mean = mu, sd = sigma_fl)
    fl <- pmax(fl, 0.1)

    mat <- .sample_mat_logistic(fl = fl, Lmat = pars$Lmat, w95 = mat_w95)

    data.table::data.table(
      sex    = sex,
      mat    = as.integer(mat),
      fl     = fl,
      age    = age,
      embryo = FALSE
    )
  }

  # ---- simulate sexes ----
  female_dt <- .simulate_sex(
    sex = "female", n = as.integer(n_female),
    L0_mean   = 47,    L0_sd   = 16.8,
    Lmat_mean = 167,   Lmat_sd = 12.4,
    tmat_mean = 14.6,  tmat_sd = 1.7,
    Linf_mean = 236.3, Linf_sd = 20.2,
    age_max = age_max_female,
    mat_w95 = mat_w95_female
  )

  male_dt <- .simulate_sex(
    sex = "male", n = as.integer(n_male),
    L0_mean   = 45,    L0_sd   = 11.3,
    Lmat_mean = 142,   Lmat_sd = 8.7,
    tmat_mean = 12.1,  tmat_sd = 1.4,
    Linf_mean = 202.9, Linf_sd = 15.8,
    age_max = age_max_male,
    mat_w95 = mat_w95_male
  )

  # ---- embryo records ----
  embryo_dt <- data.table::data.table(
    sex    = NA_character_,
    mat    = NA_real_,
    fl     = .rtnorm(as.integer(n_embryo), mean = 46.2, sd = 30.5, lower = 22, upper = 52),
    age    = NA_real_,
    embryo = TRUE
  )

  growth_data <- data.table::rbindlist(list(female_dt, male_dt, embryo_dt),
                                       use.names = TRUE, fill = TRUE)

  if (!is.null(round_digits)) {
    growth_data[, `:=`(
      fl  = round(fl,  digits = round_digits),
      age = round(age, digits = round_digits)
    )]
  }

  growth_data[]
}


