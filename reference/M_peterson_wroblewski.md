# Peterson-Wroblewski Natural Mortality Model

Computes weight-based natural mortality following Peterson & Wroblewski
(1984). Mortality scales allometrically with body weight.

## Usage

``` r
M_peterson_wroblewski(
  age,
  Linf,
  L0,
  k,
  lw_fun,
  growth_model = c("vb", "gompertz", "logistic")
)
```

## Arguments

- age:

  Numeric vector of ages at which to compute mortality.

- Linf:

  Asymptotic length.

- L0:

  Length at birth.

- k:

  Growth coefficient (native to the fitted model).

- lw_fun:

  Function mapping length to weight in grams: `lw_fun(L)`.

- growth_model:

  Character. Growth model for length prediction: `"vb"`, `"gompertz"`,
  or `"logistic"`. Default `"vb"`.

## Value

Numeric vector of instantaneous mortality rates.

## Details

The model expresses mortality as a power function of body weight:
\$\$M(W) = 1.92 \cdot W^{-0.25}\$\$ where \\W\\ is body weight in grams.

This model is growth-model-agnostic: it only requires predicted body
weight at age, which can be derived from any growth model via a
length-weight relationship. The \\-0.25\\ exponent reflects metabolic
scaling theory: metabolic rate scales approximately as \\W^{0.75}\\
(Kleiber's law), and mortality is assumed proportional to mass-specific
metabolic rate, yielding a \\W^{-0.25}\\ dependence.

Because \\k\\ does not appear in the mortality equation itself, the
native growth coefficient and native growth model should always be used
for \\L(t)\\ prediction — no VB-equivalent conversion is needed.

## References

Peterson, I., & Wroblewski, J. S. (1984). Mortality rate of fishes in
the pelagic ecosystem. *Canadian Journal of Fisheries and Aquatic
Sciences*, 41(7), 1117-1120.

## Examples

``` r
if (FALSE) { # \dontrun{
lw_fun <- function(L) 0.0001 * L^3.1  # Length in cm, weight in g
ages <- seq(0.5, 30, by = 0.5)

M <- M_peterson_wroblewski(
  age = ages, Linf = 100, L0 = 25, k = 0.1,
  lw_fun = lw_fun, growth_model = "gompertz"
)
} # }
```
