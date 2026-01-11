# Peterson-Wroblewski Natural Mortality Model

Computes weight-based natural mortality following Peterson & Wroblewski
(1984).

## Usage

``` r
M_peterson_wroblewski(age, Linf, k, t0, lw_fun)
```

## Arguments

- age:

  Numeric vector of ages at which to compute mortality.

- Linf:

  Asymptotic length.

- k:

  von Bertalanffy growth coefficient.

- t0:

  Theoretical age at length zero.

- lw_fun:

  Function mapping length to weight in grams: `lw_fun(length)`.

## Value

Numeric vector of instantaneous mortality rates.

## Details

The model expresses mortality as an allometric function of body weight:
\$\$M(W) = 1.92 W^{-0.25}\$\$ where \\W\\ is body weight in grams.

## References

Peterson, I., & Wroblewski, J. S. (1984). Mortality rate of fishes in
the pelagic ecosystem. Canadian Journal of Fisheries and Aquatic
Sciences, 41(7), 1117-1120.

## Examples

``` r
if (FALSE) { # \dontrun{
lw <- function(lt) 1e-5 * lt^3  # Weight in grams
ages <- seq(0.1, 30, length.out = 100)
M <- M_peterson_wroblewski(ages, Linf = 200, k = 0.15, t0 = -1, lw_fun = lw)
plot(ages, M, type = "l")
} # }
```
