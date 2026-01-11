# Compute Von Bertalanffy-Equivalent Growth Coefficient

Derives the von Bertalanffy growth coefficient \\k\\ from biological
milestones \\(L\_\infty, L_0, L\_{mat}, t\_{mat})\\. This enables
Chen-Watanabe mortality estimation using posteriors from *any* growth
model (von Bertalanffy, Gompertz, or Logistic).

## Usage

``` r
compute_k_vb_equivalent(Linf, L0, Lmat, tmat, warn = TRUE)
```

## Arguments

- Linf:

  Numeric vector. Asymptotic length posterior draws.

- L0:

  Numeric vector. Length at birth posterior draws.

- Lmat:

  Numeric vector. Length at maturity posterior draws.

- tmat:

  Numeric vector. Age at maturity posterior draws.

- warn:

  Logical. If `TRUE` (default) warns when draws produce invalid \\k\\
  values.

## Value

Numeric vector of VB-equivalent \\k\\ values. Invalid values (from
non-positive arguments to log) are returned as `NA`.

## Details

The key insight is that while the three growth models use different
functional forms and produce different numerical \\k\\ values, they all
estimate the same underlying biological quantities: asymptotic length,
birth size, and maturity milestones. The VB-equivalent \\k\\ is computed
as:

\$\$k\_{VB}^{equiv} = \frac{1}{t\_{mat}} \ln\left(\frac{L\_\infty -
L_0}{L\_\infty - L\_{mat}}\right)\$\$

This is the growth coefficient that would produce a von Bertalanffy
curve passing through the biological milestones \\(0, L_0)\\ and
\\(t\_{mat}, L\_{mat})\\ with asymptote \\L\_\infty\\.

When the input growth fit is von Bertalanffy with maturity-based
parameterization, this computation exactly reproduces the fitted \\k\\.
When the input is Gompertz or Logistic, it produces the VB-equivalent
\\k\\ that encodes the same biological growth information.

## Why This Matters

The von Bertalanffy model tends to produce unstable \\L\_\infty\\
estimates when data are sparse at older ages — a common situation in
elasmobranch research. Gompertz and Logistic models often provide more
reliable fits in these cases. By deriving VB-equivalent \\k\\ from
biological milestones, users can:

- Select the growth model that best fits their data

- Still use Chen-Watanabe mortality estimation

- Maintain theoretical coherence (CW was derived under VB assumptions)

## See also

[`M_chen_watanabe`](https://brian-j-moe.github.io/vitalBayes/reference/M_chen_watanabe.md)
for the L0-parameterized CW model,
[`extract_growth_parameters`](https://brian-j-moe.github.io/vitalBayes/reference/extract_growth_parameters.md)
for posterior extraction.

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract from any growth model posterior
draws <- extract_growth_parameters(growth_fit, sex = 1)
k_vb <- compute_k_vb_equivalent(
  Linf = draws$Linf,
  L0   = draws$L0,
  Lmat = draws$Lmat,
  tmat = draws$tmat
)

# Direct specification for sensitivity analysis
k_vb <- compute_k_vb_equivalent(
  Linf = rnorm(1000, 100, 5),
  L0   = rnorm(1000, 25, 2),
  Lmat = rnorm(1000, 70, 3),
  tmat = rnorm(1000, 10, 1)
)
} # }
```
