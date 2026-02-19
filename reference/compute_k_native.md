# Compute Native Growth Coefficient for a Specific Model

Derives the growth coefficient \\k\\ for a specified growth model from
biological milestones \\(L\_\infty, L_0, L\_{mat}, t\_{mat})\\.

## Usage

``` r
compute_k_native(
  Linf,
  L0,
  Lmat,
  tmat,
  growth_model = c("vb", "gompertz", "logistic"),
  warn = TRUE
)
```

## Arguments

- Linf:

  Numeric vector. Asymptotic length.

- L0:

  Numeric vector. Length at birth.

- Lmat:

  Numeric vector. Length at maturity.

- tmat:

  Numeric vector. Age at maturity.

- growth_model:

  Character. Growth model: `"vb"`, `"gompertz"`, or `"logistic"`.

- warn:

  Logical. Warn on invalid values?

## Value

Numeric vector of native k values. Invalid values returned as NA.

## Details

The native k formulas for each model are:

**Von Bertalanffy:** \$\$k = \frac{1}{t\_{mat}}
\ln\left(\frac{L\_\infty - L_0}{L\_\infty - L\_{mat}}\right)\$\$

**Gompertz:** \$\$k = -\frac{1}{t\_{mat}} \ln\left(\frac{\ln(L\_\infty /
L\_{mat})}{\ln(L\_\infty / L_0)}\right)\$\$

**Logistic:** \$\$k = -\frac{1}{t\_{mat}}
\ln\left(\frac{L\_\infty/L\_{mat} - 1}{L\_\infty/L_0 - 1}\right)\$\$

## Examples

``` r
k_vb <- compute_k_native(Linf = 126, L0 = 35, Lmat = 83, tmat = 47,
                         growth_model = "vb")
k_gomp <- compute_k_native(Linf = 108, L0 = 35, Lmat = 83, tmat = 47,
                           growth_model = "gompertz")
```
