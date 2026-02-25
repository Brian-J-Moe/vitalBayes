# Understanding vitalBayes: The Statistics Behind Elasmobranch Life History Modeling

## Introduction

Welcome to the statistical foundations of vitalBayes! This guide is
designed for researchers who want to understand not just *how* to use
the package, but *why* the models are built the way they are. Whether
you’re preparing a manuscript, explaining your methods to reviewers, or
simply curious about the statistical machinery under the hood, this
document will walk you through the key concepts step by step.

### What Makes vitalBayes Different?

Most life history analyses treat birth size estimation, maturity, and
growth as separate problems solved with separate tools. vitalBayes takes
a different approach: it recognizes that these biological processes are
fundamentally connected, and the statistical framework should reflect
that connection.

**The Core Philosophy**: An individual that matures at length
\\L\_{mat}\\ and age \\t\_{mat}\\ *must* pass through that point on its
growth curve. By building this constraint directly into our growth
models, we can use well-estimated maturity parameters to stabilize
poorly-estimated growth parameters. This is especially valuable when age
data are sparse or uncertain.

### The Three-Stage Workflow

The vitalBayes framework proceeds through three interconnected
estimation stages, with a fourth stage connecting the biological
parameters to population-level processes:

1.  **Birth Size Estimation**: Distinguishing embryos from free-swimming
    individuals to estimate size at birth
2.  **Maturity Estimation**: Modeling the probability of reproductive
    maturity as a function of length and/or age
3.  **Growth Modeling**: Fitting growth curves that incorporate
    information from the previous stages
4.  **Mortality Estimation**: Deriving age-specific natural mortality
    schedules from growth posteriors

Each stage builds on the previous one, with posterior distributions from
earlier models informing priors for later models. This creates a
coherent analytical pipeline that properly propagates uncertainty
throughout.

## Thinking About Thresholds: Why We Use the Probit Link

Before diving into specific models, let’s explore a fundamental choice
that shapes the entire vitalBayes framework: the probit link function.
Understanding this choice requires thinking carefully about what we’re
actually modeling when we analyze birth timing or maturity.

### The Threshold-Crossing Concept

Imagine you’re observing a population of developing embryos. At some
point, each individual transitions from embryo to free-swimming neonate.
But when exactly does this happen?

In reality, the timing isn’t perfectly predictable. Even among siblings
from the same litter, some might be ready to be born at slightly smaller
sizes than others due to subtle differences in development rate,
maternal condition, or environmental factors. If we could measure
“developmental readiness” on some continuous scale, we’d expect it to
vary across individuals in a roughly bell-shaped (normal) distribution.

**Think about it**: If an individual is born when its “readiness”
crosses some threshold, and readiness varies normally across
individuals, what would the relationship between size and birth
probability look like?

This intuition leads us to the **latent variable formulation**. We
imagine each individual has a latent (unobserved) readiness score:

\\z_i^\* = \beta \cdot (L_i - b\_{50}) + \epsilon_i, \quad \text{where }
\epsilon_i \sim \mathcal{N}(0, 1)\\

The individual transitions (is born, becomes mature, etc.) when this
latent variable exceeds zero:

\\y_i = \begin{cases} 1 & \text{if } z_i^\* \> 0 \\ 0 & \text{otherwise}
\end{cases}\\

Working through the probability:

\\P(y_i = 1) = P(z_i^\* \> 0) = P(\epsilon_i \> -\beta(L_i - b\_{50})) =
\Phi(\beta(L_i - b\_{50}))\\

where \\\Phi(\cdot)\\ is the cumulative distribution function of the
standard normal distribution. This is exactly the probit model!

### Why Not Logit?

You might be wondering: “But I’ve always used logistic regression for
binary outcomes. Why switch to probit?”

Both probit and logit links produce S-shaped curves that look nearly
identical in most datasets. The practical difference in fit is usually
negligible. However, the *interpretation* differs substantially:

**The Logit Interpretation**: The logit link models log-odds ratios. A
one-unit increase in length multiplies the odds of maturity by
\\e^\beta\\. While mathematically elegant, odds ratios aren’t a natural
way to think about developmental biology. When was the last time you
read a fisheries paper reporting “the odds of maturity increased by a
factor of 1.3 per centimeter”?

**The Probit Interpretation**: The probit link directly models the
threshold-crossing process described above. The slope \\\beta\\ tells
you how many standard deviations of developmental readiness correspond
to one unit of length. This maps cleanly onto biological intuition about
individual variation in development.

By using probit for both birth and maturity models, vitalBayes maintains
a consistent interpretation across the entire life history: each
transition represents crossing a normally-distributed developmental
threshold.

### Derived Quantities

One practical consequence of the probit link is how we calculate derived
quantities like “length at 5% maturity” or “length at 95% maturity.”

For a probit model, if we want the length at which probability equals
\\p\\:

\\L_p = L\_{50} + \frac{\Phi^{-1}(p)}{\beta}\\

For the commonly-reported 5% and 95% points:

- \\\Phi^{-1}(0.05) \approx -1.645\\
- \\\Phi^{-1}(0.95) \approx 1.645\\

So the **transition width** (the length range over which 90% of the
transition occurs) is:

\\\Delta\_{trans} = L\_{95} - L\_{05} = \frac{2 \times 1.645}{\beta}
\approx \frac{3.29}{\beta}\\

This is slightly narrower than the equivalent logit calculation (which
uses \\\log(19) \approx 2.944\\ instead of 1.645), meaning probit slopes
are numerically larger than logit slopes for the same curve steepness.

## Birth Size Estimation

The first stage of the vitalBayes workflow estimates birth size by
comparing lengths of embryos and free-swimming individuals.

### The Model

We observe lengths \\L_i\\ for individuals with known status \\y_i\\ (0
= embryo, 1 = free-swimming). The model is:

\\y_i \sim \text{Bernoulli}(p_i)\\ \\p_i = \Phi\[\beta \cdot (L_i -
b\_{50})\]\\

where:

- \\b\_{50}\\ is the length at which 50% of individuals have
  transitioned to free-swimming status (our primary parameter of
  interest)
- \\\beta \> 0\\ controls how sharply the transition occurs

**What \\b\_{50}\\ represents**: This isn’t exactly “birth length” in
the sense of the length at which parturition occurs. Rather, it’s the
length at which we’d expect a 50-50 chance that a randomly-selected
individual of that size is free-swimming vs. still an embryo. For
practical purposes in population modeling, this is often more useful
than trying to pinpoint an exact birth moment.

### Prior Specification

Both parameters must be positive, so we use lognormal priors:

\\\log(b\_{50}) \sim \mathcal{N}(\mu\_{b\_{50}}, \sigma\_{b\_{50}}^2)\\
\\\log(\beta) \sim \mathcal{N}(\mu\_\beta, \sigma\_\beta^2)\\

The prior mean for \\b\_{50}\\ is computed automatically from the data
as the midpoint between the largest embryo and smallest free-swimming
individual:

\\\hat{b}\_{50}^{prior} = \frac{\max\\L_i : y_i = 0\\ + \min\\L_i : y_i
= 1\\}{2}\\

This data-derived prior center helps chains start in reasonable
parameter space while remaining appropriately vague.

### Practical Considerations

#### What if there’s no overlap?

In some datasets, the largest embryo is smaller than the smallest
free-swimmer. This actually makes estimation easier—the transition is
clearly constrained to fall in that gap. The model handles this
gracefully.

#### What if there’s substantial overlap?

This is biologically interesting! It suggests high individual variation
in birth timing, which manifests as a shallow slope \\\beta\\ and wide
transition zone. The model will estimate this appropriately, though
you’ll want to think carefully about whether the “overlap” might reflect
misclassification in your data.

## Maturity Estimation

Maturity ogives are fundamental to elasmobranch population dynamics. The
vitalBayes approach extends standard methods in two key ways: using the
probit link and implementing partial pooling for two-sex models.

### Basic Maturity Model

For a single-sex model, we observe maturity status \\y_i \in \\0, 1\\\\
as a function of length or age:

\\y_i \sim \text{Bernoulli}(p_i)\\ \\p_i = \Phi\[\beta \cdot (x_i -
x\_{50})\]\\

where \\x\\ is either length \\L\\ or age \\t\\, and \\x\_{50}\\ is the
corresponding 50% maturity point.

### The Challenge of Imbalanced Sex Ratios

Here’s where things get interesting. Elasmobranch sampling often
produces highly imbalanced sex ratios due to sexual segregation by
depth, season, or habitat, gear selectivity, and differential
catchability.

**Consider this scenario**: You’ve sampled 150 female gulper sharks with
good maturity data, but only 34 males. You want sex-specific maturity
estimates. What are your options?

1.  **Pool the sexes**: Ignore sexual dimorphism entirely. Simple, but
    potentially biologically wrong.
2.  **Fit separately**: Independent estimates per sex. Honest, but the
    male estimate will have huge uncertainty.
3.  **Something in between?**: What if we could let the sexes inform
    each other while still allowing for genuine differences?

Option 3 is **partial pooling**, and it’s one of the most powerful
features of hierarchical Bayesian modeling.

### Partial Pooling: Borrowing Strength Without Losing Flexibility

The core idea is elegant: we model sex-specific parameters as draws from
a common population distribution. Mathematically:

\\\log(L\_{50,s}) = \mu\_{L\_{50}} + \tau\_{L\_{50}} \cdot \eta_s, \quad
\eta_s \sim \mathcal{N}(0, 1)\\

where:

- \\s \in \\1, 2\\\\ indexes sex (1 = female, 2 = male)
- \\\mu\_{L\_{50}}\\ is the population mean (log-scale)
- \\\tau\_{L\_{50}}\\ is the between-sex standard deviation
- \\\eta_s\\ are standardized deviates

**How partial pooling works**: The key parameter is \\\tau\_{L\_{50}}\\,
which the model estimates from the data.

- When \\\tau\\ is estimated to be **large**: The sexes are genuinely
  different, and estimates stay close to what you’d get from independent
  fitting.
- When \\\tau\\ is estimated to be **small**: The data suggest the sexes
  are similar, and estimates shrink toward each other.
- When **one sex has sparse data**: That sex’s estimate shrinks toward
  the better-estimated sex, reducing uncertainty without forcing
  equality.

The model learns the appropriate degree of pooling from the data itself!

#### Why This Matters

For our imbalanced sex example with 150 females and 34 males:

- The female \\L\_{50}\\ estimate is well-informed by data
- The male \\L\_{50}\\ estimate, if fit independently, would have a wide
  credible interval
- With partial pooling, the male estimate “borrows” information from the
  female estimate, resulting in a more precise (and usually more
  accurate) estimate
- If the males truly are different from females, the data will push
  \\\tau\\ larger and the estimates will separate

This isn’t “making up” data—it’s appropriately using biological
knowledge (that conspecific sexes are related) to improve inference.

### The Non-Centered Parameterization: A Technical But Important Detail

You might have noticed something strange in the model specification: why
write \\\mu + \tau \cdot \eta\\ instead of just sampling
\\\log(L\_{50,s}) \sim \mathcal{N}(\mu, \tau^2)\\ directly?

This is called **non-centered parameterization**, and it solves a
computational problem called the “funnel.”

**The Funnel Problem**: In hierarchical models, when \\\tau\\ is small,
the group-level parameters (here, \\L\_{50,s}\\) must all be very close
to \\\mu\\. This creates a “funnel” shape in the posterior: at small
\\\tau\\, there’s a narrow region where parameters can live, but at
larger \\\tau\\, they can spread out.

Hamiltonian Monte Carlo (HMC) samplers like Stan’s have trouble
navigating this geometry. The step size that works well in the wide part
of the funnel is too large for the narrow part, causing divergent
transitions.

**The Solution**: By parameterizing in terms of \\\eta_s \sim
\mathcal{N}(0, 1)\\ and computing \\L\_{50,s} = \exp(\mu + \tau \cdot
\eta_s)\\, we “unfold” the funnel. Now \\\eta_s\\ are always standard
normal regardless of \\\tau\\, and the sampler can explore freely.

This is why vitalBayes uses the non-centered form—it’s essential for
reliable inference, especially with the small sample sizes typical in
elasmobranch studies.

### Prior on Between-Sex Variation

The between-sex standard deviation \\\tau\\ receives a half-normal
prior:

\\\tau \sim \mathcal{N}^+(0, \sigma\_\tau^2)\\

**Why half-normal instead of half-Cauchy?**

Half-Cauchy priors are popular in hierarchical models (Gelman 2006
recommended them), but they have very heavy tails. For small samples,
this can lead to overseparation (the prior puts substantial probability
on large \\\tau\\ values, potentially inducing separation between groups
even when the data provide little evidence for it) and computational
issues (extreme \\\tau\\ values can cause numerical problems).

The half-normal prior is more regularizing—it gently pulls \\\tau\\
toward zero while still allowing the data to estimate larger values if
warranted. For the moderate sample sizes typical in elasmobranch work
(tens to low hundreds of individuals), this extra regularization
improves both estimation and computation.

## Growth Models

Growth modeling is where vitalBayes brings its integrative design to
full effect. We implement three classic growth functions, each with two
parameterization options, and connect them to the upstream birth and
maturity stages.

### The Three Growth Functions

All three models are implemented in \\L_0\\-based form, where \\L_0\\ is
length at birth—a directly observable quantity, unlike the VB parameter
\\t_0\\ (the theoretical age at length zero). This means priors can be
set from embryo or neonate measurements rather than from a mathematical
abstraction.

#### von Bertalanffy Growth Model (VBGM)

\\L(t) = L\_\infty - (L\_\infty - L_0) e^{-kt}\\

The von Bertalanffy model describes growth as a constant exponential
approach to \\L\_\infty\\. The absolute growth rate \\dL/dt\\ is highest
at birth and declines monotonically, with an inflection point at \\0.632
\\ L\_\infty\\. This pattern suits species where somatic growth is
fastest in early life and steadily decelerates as the organism
approaches asymptotic size. The VB model is the most widely used growth
function in fisheries science and provides the baseline parameterization
for the Chen-Watanabe natural mortality model (see [Section 7:
Connecting Growth to Mortality](#mortality)).

#### Gompertz Growth Model

\\L(t) = L\_\infty \exp\left\[-\ln\left(\frac{L\_\infty}{L_0}\right)
e^{-kt}\right\]\\

The Gompertz model describes growth where the *specific* growth rate
(proportional growth rate, \\\frac{1}{L}\frac{dL}{dt}\\) decreases
exponentially with time. Growth is initially rapid but slows earlier and
more abruptly than under VB dynamics, with the inflection point
occurring at \\L\_\infty / e \approx 0.368 \\ L\_\infty\\—much earlier
in ontogeny than the VB inflection. This makes the Gompertz particularly
suitable for species exhibiting rapid juvenile growth followed by
pronounced deceleration. Many small coastal elasmobranchs (bonnetheads,
Atlantic sharpnose sharks, some skate species) show this pattern, with
growth effectively ceasing near reproductive maturity.

#### Logistic Growth Model

\\L(t) = \frac{L\_\infty}{1 + \left(\frac{L\_\infty}{L_0} - 1\right)
e^{-kt}}\\

The Logistic model is symmetric around its inflection point at
\\L\_\infty / 2\\, with growth accelerating before the midpoint and
decelerating after. This sigmoidal trajectory can be appropriate for
species where early juvenile growth is initially slow (due to
nutritional constraints, habitat transitions, or ontogenetic diet
shifts) before accelerating during a rapid growth phase and then
decelerating toward asymptotic size. In practice, the Logistic model
produces the tightest approach to \\L\_\infty\\—it reaches asymptotic
size faster than either VB or Gompertz for equivalent parameter values.

**When might you prefer one model over another?**

The choice among growth models encodes assumptions about the organism’s
growth dynamics. The VB model, with maximum growth rate at birth, suits
species with continuously decelerating growth. The Gompertz, with its
earlier inflection, works well for species that show rapid juvenile
growth but slow dramatically before maturity. The Logistic, with its
symmetric inflection at \\L\_\infty / 2\\, accommodates species with an
initial slow phase followed by accelerating growth.

In many cases, all three produce similar-looking curves and LOO-CV
differences are modest. When this occurs, the VB is a reasonable default
both for consistency with the existing literature and because it
produces the most stable mortality estimates downstream. When LOO-CV
clearly favors an alternative, that model should be preferred. See
[`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md)
for a worked comparison example.

### Why \\L\_\infty\\ Must Exceed the Maximum Observed Length

This is one of those issues that seems technical at first but has deep
biological implications.

#### The Problem

When we fit a growth curve without constraints, the estimation algorithm
finds parameter values that minimize prediction error across the
dataset. For \\L\_\infty\\, this typically means fitting the asymptote
to go through the *center* of the length-at-age scatter at old ages.

But think about what this implies biologically. At any given age,
there’s natural variation in length—some individuals are larger than
average, some smaller. The largest individuals at each age represent the
upper tail of this distribution.

If the fitted \\L\_\infty\\ falls *below* these largest individuals,
we’re saying they’ve already exceeded their theoretical maximum size.
That’s biologically nonsensical!

#### The Correct Interpretation

\\L\_\infty\\ should represent the expected length of a very old
individual, or the length at which growth rate approaches zero. It’s an
*upper bound*, not a central tendency.

#### Why We Might Underestimate True \\L\_{max}\\

The observed maximum in your sample is almost certainly *not* the true
maximum in the population, for several reasons. First, finite samples
won’t capture the rarest, largest individuals. Second, the largest
individuals are also the oldest and have accumulated years of mortality
risk (mortality truncation). Third, size-selective fishing often removes
the largest individuals before they can be sampled.

#### The Delta-Gamma Parameterization

Rather than estimating \\L\_\infty\\ directly with a hard lower bound at
\\L\_{max}\\, vitalBayes estimates the **excess** above \\L\_{max}\\:

\\\delta_L = L\_\infty - L\_{max}, \quad \delta_L \> 0\\

and then reconstructs \\L\_\infty = L\_{max} + \delta_L\\ in the Stan
model.

**Why not a truncated lognormal?**

The original approach placed a lognormal prior on \\L\_\infty\\ with a
hard lower bound at \\L\_{max}\\. When the prior mean is close to the
truncation point (e.g., \\1.05 \times L\_{max}\\), the truncated
lognormal piles up density right at the boundary. This creates a flat
region in the posterior surface that HMC’s leapfrog integrator struggles
to navigate — the sampler bounces off the wall, producing poor
exploration and slow mixing. This pathology is especially pronounced for
the logistic growth model, which approaches the asymptote faster than VB
and can genuinely prefer \\L\_\infty\\ very close to \\L\_{max}\\.

The excess \\\delta_L\\ receives a **gamma prior**:

\\\delta_L \sim \text{Gamma}(\alpha\_\delta, \beta\_\delta)\\

with shape and rate parameters derived from the user’s `Linf_multiplier`
and `CV_Linf`:

\\\mu\_\delta = L\_{max} \times (\text{Linf\\multiplier} - 1)\\

\\\alpha\_\delta = \frac{1}{CV\_{L\_\infty}^2}, \quad \beta\_\delta =
\frac{1}{\mu\_\delta \times CV\_{L\_\infty}^2}\\

At the default settings (`Linf_multiplier = 1.05`, `CV_Linf = 0.20`),
this produces a gamma with mean \\0.05 \times L\_{max}\\ and shape
parameter 25 — a moderately peaked distribution centered at “5% above
observed maximum” that allows the data to pull \\\delta_L\\ toward zero
if the logistic model warrants it, or push it higher if the VB model
needs more headroom.

This approach eliminates boundary pile-ups (because \\\delta_L\\ lives
on \\(0, \infty)\\ with no truncation), provides natural positive skew
(small excesses are more probable than large ones, but the gamma’s right
tail accommodates uncertainty), and produces a lighter right tail than
lognormal (preventing biologically absurd asymptotic lengths while still
allowing flexibility). The user interface is unchanged — you still
specify `Linf_multiplier` and `CV_Linf` in the same way.

**Why gamma rather than lognormal for \\\delta_L\\?** The gamma’s
lighter right tail is appropriate here: \\\delta_L\\ represents a modest
biological excess, and placing too much prior mass on very large values
of \\\delta_L\\ would allow \\L\_\infty\\ to drift far above
observations without data support. The gamma concentrates more mass near
its mode while still being flexible enough to accommodate data-driven
deviations. The lognormal’s heavier tail would be more appropriate if we
had genuine uncertainty about the *order of magnitude* of \\\delta_L\\,
which we typically don’t — we’re uncertain about whether it’s 3% or 8%
of \\L\_{max}\\, not whether it’s 5% or 500%.

### The Maturity-Based Parameterization

Traditional growth models estimate the Brody growth coefficient \\k\\
directly from length-at-age data. This works reasonably well when you
have abundant data spanning the full range of ages. But for many
elasmobranch datasets, especially for rare or data-limited species,
\\k\\ is poorly constrained because age estimation is uncertain, sample
sizes are small, and old individuals (near the asymptote) are rare.

vitalBayes offers an alternative: **deriving \\k\\ from maturity
parameters**.

#### The Key Insight

Here’s the biological fact we exploit: an individual that matures at
length \\L\_{mat}\\ and age \\t\_{mat}\\ must, by definition, lie on the
growth curve at the point \\(t\_{mat}, L\_{mat})\\.

If we know (or can estimate) the maturity milestone, we can substitute
it into the growth equation and solve for \\k\\. The derivation is
model-specific, because each growth function has a different functional
form.

**Derivation for the von Bertalanffy model**:

Starting from: \\L\_{mat} = L\_\infty - (L\_\infty - L_0) e^{-k \cdot
t\_{mat}}\\

Rearranging: \\\frac{L\_\infty - L\_{mat}}{L\_\infty - L_0} = e^{-k
\cdot t\_{mat}}\\

Taking logs: \\-k \cdot t\_{mat} = \ln\left(\frac{L\_\infty -
L\_{mat}}{L\_\infty - L_0}\right)\\

Solving for \\k\\: \\k\_{VB} = \frac{1}{t\_{mat}}
\ln\left(\frac{L\_\infty - L_0}{L\_\infty - L\_{mat}}\right)\\

#### Gompertz and Logistic Derivations

The same substitution logic applied to the Gompertz and Logistic
equations yields model-specific \\k\\ expressions. For the Gompertz,
substituting \\L(t\_{mat}) = L\_{mat}\\ and solving gives:

\\k_g = \frac{1}{t\_{mat}} \ln\\\left(\frac{\ln(L\_\infty /
L_0)}{\ln(L\_\infty / L\_{mat})}\right)\\

And for the Logistic:

\\k_l = \frac{1}{t\_{mat}} \ln\\\left(\frac{L\_{mat}(L\_\infty -
L_0)}{L_0(L\_\infty - L\_{mat})}\right)\\

In each case, \\k\\ is a deterministic function of \\(L\_\infty, L_0,
L\_{mat}, t\_{mat})\\—parameters that are either directly estimated in
the Stan model or informed by upstream birth and maturity fits. This
eliminates the need to specify a prior on \\k\\ directly and breaks the
\\L\_\infty\\-\\k\\ posterior correlation.

An important subtlety: \\k\\ has different meanings across the three
models and should not be directly compared between them. For example,
Gompertz \\k_g\\ is generally larger than VB \\k\_{VB}\\ for the same
species, because the Gompertz’s exponentially declining growth rate
requires a faster initial rate constant to reach the same asymptote.
What matters for biological inference is the resulting growth
trajectory, not the numerical value of \\k\\.

#### Why This Helps

The maturity-based parameterization offers several advantages:

1.  **Reduced correlation**: The traditional \\(L\_\infty, k)\\
    parameter pair is notoriously correlated—when \\L\_\infty\\ is high,
    \\k\\ tends to be low, and vice versa. This ridge-like posterior
    geometry slows MCMC sampling. By replacing \\k\\ with \\(L\_{mat},
    t\_{mat})\\, we break this correlation.

2.  **Informative anchoring**: The maturity milestone typically falls
    *within* the observed data range, unlike \\L\_\infty\\ which
    extrapolates beyond observations. This anchoring stabilizes the
    curve.

3.  **Prior information propagation**: When maturity models are fit
    first (Stage 2), their posteriors provide informative priors for the
    growth model. Uncertainty from maturity estimation flows naturally
    into growth uncertainty.

4.  **Biological interpretability**: \\(L\_{mat}, t\_{mat})\\ have
    direct biological meaning that researchers can validate against
    independent studies.

**When to use maturity-based vs. k-based parameterization?**

Use **maturity-based** when you have maturity data and have already fit
maturity ogives, when age data are sparse or uncertain, when you want to
create a biologically integrated analysis, or when you’re working with a
data-limited species.

Use **k-based** when you don’t have maturity data, when you have
abundant high-quality age data spanning the full growth trajectory, when
you want to compare with historical studies that report \\k\\, or when
you’re doing sensitivity analyses to compare parameterization
approaches.

### Selective Pooling in Growth Models

When partial pooling is combined with the maturity-based
parameterization, a subtle problem arises. If the upstream maturity
models
([`fit_bayesian_maturity()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_maturity.md))
were themselves fit with `use_pooling = TRUE`, the maturity parameters
(\\L\_{mat}\\, \\t\_{mat}\\) already contain pooled estimates. Pooling
them *again* in the growth model creates a **double-pooling** problem
that can over-shrink sex differences toward the population mean, in
extreme cases reversing genuine biological dimorphism.

**The Double-Pooling Problem**: Consider a species where females mature
at 72 cm and males at 65 cm. If the maturity model already pulled these
estimates slightly toward each other (say to 71 and 66), and then the
growth model applies hierarchical pooling again, the estimates might
compress further (say to 70 and 67). Each stage of pooling is
individually appropriate, but the cumulative effect can wash out a real
biological signal.

vitalBayes addresses this with **selective pooling**: when
maturity-based parameterization is used with partial pooling enabled,
the default behavior pools only \\L\_\infty\\ and \\L_0\\
hierarchically, while \\L\_{mat}\\ and \\t\_{mat}\\ receive direct
sex-specific priors from the upstream maturity fits. This is controlled
by the `pool_maturity` argument in
[`fit_bayesian_growth()`](https://brian-j-moe.github.io/vitalBayes/reference/fit_bayesian_growth.md),
which auto-detects whether the maturity inputs are vitalBayes model fits
and adjusts accordingly.

The logic is straightforward: parameters that have already been
regularized by a hierarchical model shouldn’t be regularized again.
Parameters that haven’t (\\L\_\infty\\, \\L_0\\) still benefit from
pooling. For the full decision guide and code examples, see
[`vignette("fit_bayesian_growth")`](https://brian-j-moe.github.io/vitalBayes/articles/fit_bayesian_growth.md).

### Observation Model

Growth observations are modeled with lognormal errors:

\\\log(L_i) \sim \mathcal{N}(\log(\mu_i), \sigma^2)\\

where \\\mu_i = L(t_i; \theta)\\ is the predicted length from the growth
model.

**Why lognormal instead of normal errors?**

Lengths can’t be negative, and normal errors could produce negative
predictions while lognormal can’t. Measurement error is often
*proportional* to size—a 1 cm error matters more for a 30 cm juvenile
than a 100 cm adult, and lognormal naturally captures this
multiplicative structure. On the log scale, we’re essentially fitting a
linear model with constant variance, which is statistically
well-behaved.

For datasets with outliers or heavy-tailed residuals, vitalBayes also
supports a **Student-t observation model** via the `robust = TRUE`
argument. The Student-t distribution has heavier tails than the normal
(lognormal), which downweights extreme observations rather than letting
them dominate the posterior. This is particularly useful for species
where occasional extreme measurements (gear-related injuries,
measurement errors, or genuinely unusual individuals) would otherwise
exert undue influence on parameter estimates.

## Prior Elicitation from CV: Making Priors Intuitive

One of the package’s practical design choices is allowing users to
specify priors using **coefficients of variation (CV)** on the natural
scale, rather than standard deviations on the log scale.

### The Problem with Log-Scale Priors

Traditional Bayesian software asks for priors like: “The log of
\\L\_\infty\\ has a normal prior with mean 4.5 and SD 0.3.”

But what does that *mean* in terms of actual lengths? It’s not
immediately obvious. And when parameters span different scales (birth
length vs. asymptotic length vs. growth rate), using consistent SDs on
the log scale produces inconsistent prior beliefs on the natural scale.

### The CV Solution

The coefficient of variation is the ratio of standard deviation to mean:

\\CV = \frac{\sigma}{\mu}\\

This is **scale-invariant**: a CV of 0.2 always means “20% relative
uncertainty,” whether you’re talking about a parameter that’s 30 cm or
300 cm.

vitalBayes accepts priors specified as a natural-scale mean (your best
guess in interpretable units) and a CV (your uncertainty relative to
that mean).

**Default CVs in vitalBayes**:

| Parameter                             | Default CV | Context                                                         |
|:--------------------------------------|:-----------|:----------------------------------------------------------------|
| \\L\_\infty\\                         | 0.20       | Growth model; informed by \\L\_{max}\\                          |
| \\L_0\\                               | 0.30       | Growth model; fewer neonatal observations                       |
| \\k\\                                 | 0.50       | Growth model (k-based only); often poorly constrained           |
| \\L\_{mat}\\                          | 0.20       | Growth model (maturity-based); well-estimated from maturity fit |
| \\t\_{mat}\\                          | 0.30       | Growth model (maturity-based); more uncertain than \\L\_{mat}\\ |
| \\b\_{50}\\, \\L\_{50}\\, \\t\_{50}\\ | 0.30       | Birth and maturity models                                       |

These defaults encode reasonable beliefs for data-limited species but
can be adjusted based on prior knowledge. When upstream model fits are
provided (e.g., `birth_stanfit`, `length.mature_stanfit`), their
posterior summaries replace the CV-based defaults.

### What Are the Actual Prior Distributions?

Most positive parameters in vitalBayes receive **lognormal priors**:

\\\log(\theta) \sim \mathcal{N}(\mu\_{log}, \sigma\_{log}^2)\\

This applies to \\b\_{50}\\, \\\beta\\ (birth model); \\L\_{50}\\,
\\t\_{50}\\, slope (maturity models); and \\L_0\\, \\k\\, \\L\_{mat}\\,
\\t\_{mat}\\ (growth models). Lognormal priors are a natural choice for
strictly positive biological quantities because they enforce positivity,
they produce right-skewed distributions on the natural scale (consistent
with the idea that a parameter is more likely to be somewhat larger than
your best guess than dramatically smaller), and they are conjugate-ish
with the log-scale observation models used throughout the package.

The CV specification is converted to log-scale hyperparameters
\\(\mu\_{log}, \sigma\_{log})\\ before being passed to Stan (see below).

\\L\_\infty\\ is the exception. Rather than a truncated lognormal, it
receives a **gamma prior on the excess above \\L\_{max}\\** — the
delta-gamma parameterization described in the [\\L\_\infty\\ section
above](#linf).

The observation error \\\sigma\\ receives a half-normal prior: \\\sigma
\sim \text{Half-Normal}(0, s)\\, where \\s\\ is a scale parameter
(default 1). Hierarchical between-sex standard deviations \\\tau\\ also
receive half-normal priors (see the [partial pooling
section](#partial-pooling)).

### Converting CV to Log-Scale Priors

The user specifies a natural-scale mean \\\mu\_\theta\\ and
\\CV\_\theta\\, producing \\\sigma\_\theta = CV\_\theta \times
\mu\_\theta\\. These must be converted to log-scale hyperparameters for
the lognormal prior in Stan. The analytical approximation is:

\\\mu\_{log} \approx \log(\mu\_\theta) - \frac{\sigma\_{log}^2}{2},
\quad \sigma\_{log} \approx CV\_\theta\\

This works well for small CVs (\\\lesssim 0.3\\), but breaks down at
larger values because the lognormal is increasingly skewed and the
relationship between natural-scale moments and log-scale parameters
becomes nonlinear.

vitalBayes instead uses **simulation-based moment matching**: it draws
samples from a truncated normal on the natural scale (lower bound at
zero, ensuring positivity), log-transforms them, and takes the empirical
mean and SD of the log-transformed samples as \\\mu\_{log}\\ and
\\\sigma\_{log}\\. This properly handles the positivity constraint and
produces accurate log-scale priors even when the CV is large enough
(e.g., \\CV_k = 0.5\\) that the analytical approximation would introduce
meaningful bias.

## Model Assessment and Comparison

Fitting a model is only half the battle. We also need to assess how well
it fits the data and compare alternative model formulations.

### Posterior Predictive Checks

The fundamental question: **does the model produce data that look like
the real data?**

For each posterior draw \\\theta^{(m)}\\, we can generate a replicated
dataset:

\\L_i^{rep(m)} \sim \text{LogNormal}(\log(\mu_i^{(m)}), \sigma^{(m)})\\

Comparing these replications to the observed data reveals model
problems: systematic bias (do residuals trend with age or fitted
values?), underdispersion (is the model overconfident, with observations
falling outside predicted intervals too often?), and overdispersion (is
the model too uncertain?).

vitalBayes computes several summary metrics including RMSE (root mean
squared prediction error), MAE (mean absolute error), coverage
(proportion of observations within 95% predictive intervals), and
calibration (for binary models, whether predicted probabilities match
observed frequencies).

### Model Comparison with LOO-CV

When comparing multiple models (VB vs. Gompertz vs. Logistic, or
maturity-based vs. k-based), we need a principled way to assess which
explains the data best while accounting for model complexity.

**Leave-one-out cross-validation (LOO-CV)** answers the question: “How
well does the model predict each observation when that observation is
left out of fitting?”

**The expected log pointwise predictive density (ELPD)**:

\\\text{elpd}\_{loo} = \sum\_{i=1}^{n} \log p(y_i \| y\_{-i})\\

This measures how surprised the model would be by each observation if it
hadn’t seen it. Higher (less negative) is better.

vitalBayes uses Pareto-smoothed importance sampling (PSIS-LOO) for
efficient computation without actually refitting the model \\n\\ times.

The `loo` package provides elpd_diff (difference in ELPD from the best
model), se_diff (standard error of that difference), and LOOIC (LOO
information criterion, \\-2 \times \text{elpd}\_{loo}\\, analogous to
AIC).

A rule of thumb: if `elpd_diff` is more than 2–3 times its standard
error, the model difference is meaningful.

## Connecting Growth to Mortality

The vitalBayes analytical pipeline extends beyond growth estimation.
Growth posteriors—specifically \\L\_\infty\\, \\L_0\\, and \\k\\—provide
the raw material for deriving age-specific natural mortality schedules,
which are essential inputs for population models and stock assessments.

### Why Growth Parameters Matter for Mortality

The fundamental biological link is that growth and mortality are both
governed by metabolic processes. Larger, slower-growing organisms tend
to experience lower natural mortality rates. This relationship has been
formalized in several empirical and semi-theoretical mortality
estimators, each of which takes growth curve parameters as inputs.

**The Chen-Watanabe framework** is the primary mortality model in
vitalBayes. It derives age-specific mortality \\M(t)\\ from the growth
trajectory:

\\M(t) = \frac{M\_\infty \cdot L\_\infty}{L(t)}\\

where \\M\_\infty\\ is the asymptotic mortality rate (the mortality rate
at maximum body size) and \\L(t)\\ is the growth function. As an
individual grows, its mortality declines in inverse proportion to its
length—large fish die at lower rates than small fish.

The critical insight is that \\M\_\infty = k\_{VB}\\, where \\k\_{VB}\\
is the von Bertalanffy growth coefficient. Chen and Watanabe derived
this relationship under VB assumptions: as \\t \to \infty\\, \\L(t) \to
L\_\infty\\ and \\M(t) \to k\\. This means the growth coefficient
governs both the rate of growth and the floor of mortality.

#### Growth-Model-Agnostic Mortality

A design choice in vitalBayes is that mortality estimation works with
any of the three supported growth models, not just von Bertalanffy. The
package introduces a normalized growth coefficient \\G(t) = L(t) /
L\_\infty\\ that allows the CW mortality model to be written as:

\\M(t) = \frac{M\_\infty}{G(t)}\\

where \\G(t)\\ takes a different functional form for each growth model
but produces the same mortality behavior: high mortality at young ages
that declines as the organism grows. Because \\M\_\infty = k\_{VB}\\
regardless of which growth model fits the data, non-VB models compute a
VB-equivalent \\k\\ from their own biological milestones. For the full
mathematical framework, see
`vignette("chen_watanabe_reparameterization")`.

#### Stochastic Mortality Estimation

In practice, uncertainty in growth parameters propagates directly into
mortality estimates via Monte Carlo sampling. The
[`get_stochastic_mortality()`](https://brian-j-moe.github.io/vitalBayes/reference/get_stochastic_mortality.md)
function draws posterior samples from the growth model, computes
age-specific mortality for each draw, and returns the full posterior
distribution of \\M(t)\\. When a fitted growth model is available, draws
come from the joint posterior preserving all parameter correlations. For
literature-based analyses, manual parameters can be specified as mean-SD
pairs, with optional bivariate sampling of \\L\_{mat}\\ and \\t\_{mat}\\
via Cholesky factorization to preserve their positive biological
correlation. This is covered in detail in
[`vignette("mortality_estimation")`](https://brian-j-moe.github.io/vitalBayes/articles/mortality_estimation.md),
which also describes two additional mortality estimators
(Peterson-Wroblewski and Lorenzen) and scaling procedures for
calibrating model-derived schedules to empirical estimates.

## Putting It All Together: The Integrated Workflow

Let’s trace how information flows through a complete vitalBayes
analysis.

### Stage 1 → Stage 2: Birth Informs Growth

The birth model produces a posterior for \\b\_{50}\\ (size at 50%
free-swimming probability). This becomes the **prior for \\L_0\\**
(length at birth) in growth models. The posterior uncertainty from birth
estimation carries forward, so that poorly resolved birth size produces
appropriately wider growth posteriors.

### Stage 2 → Stage 3: Maturity Informs Growth

The maturity models produce posteriors for \\L\_{50}\\ and \\t\_{50}\\.
For maturity-based growth parameterization, these become priors for
\\L\_{mat}\\ and \\t\_{mat}\\, with the growth coefficient \\k\\ derived
rather than estimated.

This integration means that uncertainty propagates (poor birth or
maturity data leads to wider growth posteriors), that stages reinforce
each other (good maturity data stabilizes growth estimation), and that
biological consistency emerges naturally (a single coherent story about
the species’ life history).

### Stage 3 → Stage 4: Growth Informs Mortality

Growth posteriors provide the parameters needed for natural mortality
estimation. The growth coefficient \\k\\ determines the asymptotic
mortality rate, \\L\_\infty\\ sets the scale, and the growth trajectory
\\L(t)\\ governs how mortality declines with age. The full posterior
from the growth model propagates through, so that mortality schedules
carry appropriate uncertainty from all upstream stages.

**Why not fit everything simultaneously?**

A fully joint model fitting birth, maturity, growth, and mortality
together would be conceptually elegant, but computationally prohibitive.
The joint posterior would be very high-dimensional, it would be harder
to identify which component is misbehaving during diagnostics, and
researchers often have different data for different stages. The
sequential approach is a principled approximation that maintains
biological coherence while keeping computation tractable and each stage
individually diagnosable.

## Summary: Key Takeaways

The vitalBayes framework rests on several foundational design choices:

1.  **Probit link functions** for birth and maturity models, reflecting
    a threshold-crossing interpretation of developmental transitions

2.  **Partial pooling** for two-sex models, allowing adaptive shrinkage
    that stabilizes sparse-sex estimates without forcing equality, with
    **selective pooling** in growth models to avoid double-pooling
    maturity parameters

3.  **Non-centered parameterization** for hierarchical models, enabling
    reliable MCMC sampling even with small samples

4.  **Maturity-based growth parameterization** that derives \\k\\ from
    observable maturity metrics, with model-specific derivations for von
    Bertalanffy, Gompertz, and Logistic growth functions

5.  **CV-based prior elicitation** for intuitive specification of prior
    uncertainty

6.  **Proper \\L\_\infty\\ constraints** ensuring the asymptote is
    interpretable as an upper bound, not a central tendency

7.  **Lognormal observation models** with optional Student-t robustness,
    capturing the multiplicative nature of biological measurement error

8.  **Growth-to-mortality integration** connecting individual growth
    trajectories to age-specific natural mortality schedules via the
    Chen-Watanabe and other empirical estimators

Each choice was made deliberately to address specific challenges in
elasmobranch life history estimation: sparse samples, imbalanced sex
ratios, uncertain aging, and the need for biologically coherent
analysis.

## Further Reading

For those wanting to dive deeper:

**Bayesian Foundations**:

- McElreath, R. (2020). *Statistical Rethinking* (2nd ed.). An excellent
  introduction to Bayesian thinking.
- Gelman et al. (2013). *Bayesian Data Analysis* (3rd ed.). The
  comprehensive reference.

**Hierarchical Models**:

- Betancourt & Girolami (2015). “Hamiltonian Monte Carlo for
  Hierarchical Models.” Technical but invaluable for understanding
  non-centered parameterization.

**Growth Models**:

- Smart et al. (2016). “Multimodel approaches in shark and ray growth
  studies.” *Fish and Fisheries*. Reviews growth model selection.
- Cailliet et al. (2006). “Age and growth studies of chondrichthyan
  fishes.” *Environmental Biology of Fishes*. Classic methodology
  reference.

**Natural Mortality**:

- Chen & Watanabe (1989). “Age dependence of natural mortality
  coefficient in fish population dynamics.” *Nippon Suisan Gakkaishi*.
  Foundation for the CW reparameterization.
- Lorenzen (1996). “The relationship between body weight and natural
  mortality in juvenile and adult fish: a comparison of natural
  ecosystems and aquaculture.” *Journal of Fish Biology*. Weight-based
  mortality estimator.
- Then et al. (2015). “Evaluating the predictive performance of
  estimators of natural mortality rate.” *Methods in Ecology and
  Evolution*. Empirical comparison of mortality estimators.

**Stan and cmdstanr**:

- [Stan User’s Guide](https://mc-stan.org/users/documentation/).
  Comprehensive documentation.
- [cmdstanr documentation](https://mc-stan.org/cmdstanr/). R interface
  guide.

------------------------------------------------------------------------

*This document is part of the vitalBayes R package. For bug reports,
feature requests, or questions, please visit the GitHub repository.*
