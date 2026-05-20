# Simulate artificial von Bertalanffy growth data with truncated-exponential ages and sex-specific logistic maturity ogives

Simulates individual-level VB parameters (L0, Lmat, tmat, Linf) by sex,
draws ages from a truncated exponential (older ages increasingly rare),
computes VB mean length, adds observation error to obtain fl, then
samples maturity from a logistic ogive on observed fl with sex-specific
steepness (w95).

## Usage

``` r
simulate_vb_growth_data(
  n_female = 150L,
  n_male = 120L,
  n_embryo = 26L,
  age_max_female = 35,
  age_max_male = 31,
  age_tail_frac = 0.7,
  age_tail_prob = 0.3,
  sigma_fl = 12,
  mat_w95_female = 20,
  mat_w95_male = 18,
  seed = 123,
  round_digits = 1L
)
```

## Arguments

- n_female:

  Integer. Number of female (non-embryo) observations.

- n_male:

  Integer. Number of male (non-embryo) observations.

- n_embryo:

  Integer. Number of embryo observations.

- age_max_female:

  Numeric. Maximum age (years) for females.

- age_max_male:

  Numeric. Maximum age (years) for males.

- age_tail_frac:

  Numeric in (0,1). Reference age = age_tail_frac \* age_max for tuning
  the tail.

- age_tail_prob:

  Numeric in (0,1). Approx target P(age \> age_tail_frac\*age_max) for
  the underlying exponential.

- sigma_fl:

  Numeric. SD of observation error (cm) added to VB mean length.

- mat_w95_female:

  Numeric. Female ogive width (cm) between 5% and 95% maturity.

- mat_w95_male:

  Numeric. Male ogive width (cm) between 5% and 95% maturity.

- seed:

  Integer or NULL. RNG seed for reproducibility.

- round_digits:

  Integer or NULL. If not NULL, round fl and age to this many digits.

## Value

data.table with columns: sex, mat, fl, age, embryo

## Details

VB formulation: k = (1/tmat) \* log((Linf - L0) / (Linf - Lmat)) mu =
Linf - (Linf - L0) \* exp(-k \* age)
