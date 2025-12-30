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

The vitalBayes framework proceeds through three interconnected stages:

1.  **Birth Size Estimation**: Distinguishing embryos from free-swimming
    individuals to estimate size at birth
2.  **Maturity Estimation**: Modeling the probability of reproductive
    maturity as a function of length and/or age
3.  **Growth Modeling**: Fitting growth curves that incorporate
    information from the previous stages

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
read a shark paper reporting “the odds of maturity increased by a factor
of 1.3 per centimeter”?

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
probit link (for consistency with the birth model) and implementing
partial pooling for two-sex models.

### Basic Maturity Model

For a single-sex model, we observe maturity status \\y_i \in \\0, 1\\\\
as a function of length or age:

\\y_i \sim \text{Bernoulli}(p_i)\\ \\p_i = \Phi\[\beta \cdot (x_i -
x\_{50})\]\\

where \\x\\ is either length \\L\\ or age \\t\\, and \\x\_{50}\\ is the
corresponding 50% maturity point.

### The Challenge of Imbalanced Sex Ratios

Here’s where things get interesting. Elasmobranch sampling often
produces highly imbalanced sex ratios due to:

- Sexual segregation by depth, season, or habitat
- Gear selectivity
- Differential catchability

**Consider this scenario**: You’ve sampled 150 female gulper sharks with
good maturity data, but only 23 males. You want sex-specific maturity
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

For our gulper shark example with 150 females and 23 males:

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

This is called the **non-centered parameterization**, and it solves a
nasty computational problem called the “funnel.”

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
this can lead to:

1.  **Overseparation**: The prior puts substantial probability on large
    \\\tau\\ values, potentially inducing separation between groups even
    when the data provide little evidence for it
2.  **Computational issues**: Extreme \\\tau\\ values can cause
    numerical problems

The half-normal prior is more regularizing—it gently pulls \\\tau\\
toward zero while still allowing the data to estimate larger values if
warranted. For the moderate sample sizes typical in elasmobranch work
(tens to low hundreds of individuals), this extra regularization
improves both estimation and computation.

## Growth Models

Growth modeling is where vitalBayes really shows its integrative
approach. We implement three classic growth functions, each with two
parameterization options.

### The Three Growth Functions

#### von Bertalanffy Growth Model (VBGM)

The workhorse of fisheries biology:

\\L(t) = L\_\infty - (L\_\infty - L_0) e^{-kt}\\

**Biological interpretation**: Growth rate is proportional to the
remaining “growth potential” (\\L\_\infty - L\\). As an individual
approaches its asymptotic size, growth slows down. This reflects
metabolic constraints where anabolic processes (building tissue) scale
with surface area while catabolic processes (maintaining tissue) scale
with volume.

#### Gompertz Growth Model

\\L(t) = L\_\infty \exp\left\[-\ln\left(\frac{L\_\infty}{L_0}\right)
e^{-kt}\right\]\\

**Biological interpretation**: The *specific* growth rate (proportional
growth rate, \\\frac{1}{L}\frac{dL}{dt}\\) decreases exponentially with
time. This represents a gradual “exhaustion” of growth capacity.

#### Logistic Growth Model

\\L(t) = \frac{L\_\infty}{1 + \left(\frac{L\_\infty}{L_0} - 1\right)
e^{-kt}}\\

**Biological interpretation**: Growth rate is maximized at intermediate
sizes, representing some form of density-dependent or resource-dependent
feedback.

**When might you prefer one model over another?**

All three models produce similar-looking curves for many datasets, but
they differ in their curvature and inflection points. The Gompertz and
Logistic have an inflection point (where growth rate is maximized) at
intermediate ages, while the VBGM has maximum growth rate at birth. For
species with rapid early growth that decelerates continuously, the VBGM
often fits best. For species with a more pronounced growth spurt at
intermediate ages, Gompertz or Logistic may be preferred.

vitalBayes makes it easy to fit all three and compare them using
leave-one-out cross-validation (LOO-CV).

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

**The vitalBayes Solution**: We constrain \\L\_\infty \> L\_{max}\\ (the
maximum observed length) in the Stan model:

``` stan
real<lower=Linf_lower> log_Linf;  // where Linf_lower = log(Lmax)
```

This ensures that:

1.  The asymptote always exceeds observed data
2.  \\L\_\infty\\ is interpretable as a true upper bound
3.  The model doesn’t produce the absurdity of individuals “exceeding”
    their maximum possible size

#### Why We Might Underestimate True \\L\_{max}\\

The observed maximum in your sample is almost certainly *not* the true
maximum in the population, for several reasons:

1.  **Sampling limitation**: Finite samples won’t capture the rarest,
    largest individuals
2.  **Mortality truncation**: The largest individuals are also the
    oldest, and have accumulated years of mortality risk
3.  **Fishing pressure**: Size-selective fishing often removes the
    largest individuals before they can be sampled

For these reasons, vitalBayes sets the prior mean for \\L\_\infty\\ at
\\1.05 \times L\_{max}\\ (5% above observed maximum), reflecting the
expectation that true asymptotic size modestly exceeds what we’ve
observed.

### The Maturity-Based Parameterization: Our Key Innovation

Traditional growth models estimate the Brody growth coefficient \\k\\
directly from length-at-age data. This works reasonably well when you
have abundant data spanning the full range of ages. But for many
elasmobranch datasets, especially for rare or data-limited species,
\\k\\ is poorly constrained because:

1.  Age estimation is uncertain
2.  Sample sizes are small
3.  Old individuals (near the asymptote) are rare

vitalBayes offers an alternative: **deriving \\k\\ from maturity
parameters**.

#### The Key Insight

Here’s the biological fact we exploit: an individual that matures at
length \\L\_{mat}\\ and age \\t\_{mat}\\ must, by definition, lie on the
growth curve at the point \\(t\_{mat}, L\_{mat})\\.

If we know (or can estimate) the maturity milestone, we can substitute
it into the growth equation and solve for \\k\\.

**Derivation for the von Bertalanffy model**:

Starting from: \\L\_{mat} = L\_\infty - (L\_\infty - L_0) e^{-k \cdot
t\_{mat}}\\

Rearranging: \\\frac{L\_\infty - L\_{mat}}{L\_\infty - L_0} = e^{-k
\cdot t\_{mat}}\\

Taking logs: \\-k \cdot t\_{mat} = \ln\left(\frac{L\_\infty -
L\_{mat}}{L\_\infty - L_0}\right)\\

Solving for \\k\\: \\k = \frac{1}{t\_{mat}} \ln\left(\frac{L\_\infty -
L_0}{L\_\infty - L\_{mat}}\right)\\

Similar derivations exist for Gompertz and Logistic models.

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

Use **maturity-based** when:

- You have maturity data and have already fit maturity ogives
- Age data are sparse or uncertain
- You want to create a biologically integrated analysis
- You’re working with a data-limited species

Use **k-based** when:

- You don’t have maturity data
- You have abundant, high-quality age data spanning the full growth
  trajectory
- You want to compare with historical studies that report k
- You’re doing sensitivity analyses to compare parameterization
  approaches

### Observation Model: Why Lognormal?

Growth observations are modeled with lognormal errors:

\\\log(L_i) \sim \mathcal{N}(\log(\mu_i), \sigma^2)\\

where \\\mu_i = L(t_i; \theta)\\ is the predicted length from the growth
model.

**Why lognormal instead of normal errors?**

1.  **Positive support**: Lengths can’t be negative. Normal errors could
    produce negative predictions; lognormal can’t.

2.  **Multiplicative errors**: Measurement error is often *proportional*
    to size. A 1 cm error matters more for a 30 cm juvenile than a 100
    cm adult. Lognormal naturally captures this.

3.  **Heteroscedasticity**: With normal errors, variance is constant
    across all sizes. With lognormal, variance increases with predicted
    length, matching empirical patterns.

On the log scale, we’re essentially fitting a linear model with constant
variance, which is statistically well-behaved.

## Prior Elicitation from CV: Making Priors Intuitive

One of vitalBayes’s practical innovations is allowing users to specify
priors using **coefficients of variation (CV)** on the natural scale,
rather than standard deviations on the log scale.

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

vitalBayes accepts priors specified as:

1.  **Natural-scale mean**: Your best guess in interpretable units (cm,
    years, etc.)
2.  **CV**: Your uncertainty relative to that mean

**Default CVs in vitalBayes**:

| Parameter | Default CV | Rationale |
|:---|:---|:---|
| \\L\_\infty\\ | 0.20 | Usually informed by \\L\_{max}\\; moderate uncertainty |
| \\L_0\\ | 0.30 | Fewer neonatal observations; higher uncertainty |
| \\k\\ | 0.50 | Often poorly constrained; high uncertainty |
| \\L\_{50}\\, \\t\_{50}\\ | 0.30 | Transition zone variability |
| \\b\_{50}\\ | 0.30 | Embryo sample rarity |

These defaults encode reasonable beliefs for data-limited species but
can be adjusted based on prior knowledge.

### Converting to Log-Scale

Under the hood, vitalBayes converts these CV-based specifications to
appropriate log-scale parameters for the lognormal priors used in Stan.
The conversion uses simulation-based moment matching:

1.  Draw samples from a truncated normal on the natural scale
2.  Take logarithms
3.  Compute empirical mean and SD of the log-samples

This approach properly handles the positivity constraint and produces
priors that accurately reflect the intended uncertainty.

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
problems:

- **Systematic bias**: Do residuals trend with age or fitted values?
- **Underdispersion**: Is the model overconfident? Do observations fall
  outside predicted intervals too often?
- **Overdispersion**: Is the model too uncertain?

vitalBayes computes several summary metrics:

- **RMSE**: Root mean squared prediction error
- **MAE**: Mean absolute error
- **Coverage**: Proportion of observations within 95% predictive
  intervals
- **Calibration**: For binary models, do predicted probabilities match
  observed frequencies?

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

The `loo` package provides:

- **elpd_diff**: Difference in ELPD from the best model
- **se_diff**: Standard error of that difference
- **LOOIC**: LOO information criterion (\\-2 \times
  \text{elpd}\_{loo}\\), analogous to AIC

A rule of thumb: if `elpd_diff` is more than 2-3 times its standard
error, the model difference is meaningful.

## Putting It All Together: The Integrated Workflow

Let’s trace how information flows through a complete vitalBayes
analysis.

### Stage 1 → Stage 2: Birth Informs Maturity

The birth model produces a posterior for \\b\_{50}\\ (size at 50%
free-swimming probability). This becomes the **prior for \\L_0\\**
(length at birth) in growth models. The posterior uncertainty from birth
estimation carries forward.

### Stage 2 → Stage 3: Maturity Informs Growth

The maturity models produce posteriors for \\L\_{50}\\ and \\t\_{50}\\.
For maturity-based growth parameterization, these become priors for
\\L\_{mat}\\ and \\t\_{mat}\\, with the growth coefficient \\k\\ derived
rather than estimated.

This integration means that:

1.  **Uncertainty propagates**: Poor birth or maturity data leads to
    wider growth posteriors
2.  **Stages reinforce each other**: Good maturity data stabilizes
    growth estimation
3.  **Biological consistency**: A single coherent story emerges about
    the species’ life history

**Why not fit everything simultaneously?**

A fully joint model fitting birth, maturity, and growth together would
be conceptually elegant, but:

1.  **Computational cost**: The joint posterior would be very
    high-dimensional
2.  **Diagnostics**: It’s harder to identify which component is
    misbehaving
3.  **Modularity**: Researchers often have different data for different
    stages
4.  **Sequential priors**: Weakly informative priors at each stage
    compound appropriately

The three-stage approach is a principled approximation that maintains
biological coherence while keeping computation tractable.

## Summary: Key Takeaways

The vitalBayes framework rests on several foundational choices:

1.  **Probit link functions** for birth and maturity models, reflecting
    a threshold-crossing interpretation of developmental transitions

2.  **Partial pooling** for two-sex models, allowing adaptive shrinkage
    that stabilizes sparse-sex estimates without forcing equality

3.  **Non-centered parameterization** for hierarchical models, enabling
    reliable MCMC sampling even with small samples

4.  **Maturity-based growth parameterization** that derives \\k\\ from
    observable maturity metrics, reducing parameter correlation and
    enabling biological integration

5.  **CV-based prior elicitation** for intuitive specification of prior
    uncertainty

6.  **Proper \\L\_\infty\\ constraints** ensuring the asymptote is
    interpretable as an upper bound, not a central tendency

7.  **Lognormal observation models** capturing the multiplicative nature
    of biological measurement error

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

**Stan and cmdstanr**:

- [Stan User’s Guide](https://mc-stan.org/users/documentation/).
  Comprehensive documentation.
- [cmdstanr documentation](https://mc-stan.org/cmdstanr/). R interface
  guide.

------------------------------------------------------------------------

*This document is part of the vitalBayes R package. For bug reports,
feature requests, or questions, please visit the GitHub repository.*
