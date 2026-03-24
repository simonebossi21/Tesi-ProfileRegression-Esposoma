# Bayesian Profile Regression on HELIX Exposome Data

Identification of childhood exposure profiles associated with overweight/obesity risk using Bayesian Profile Regression (Molitor et al., 2010) applied to the HELIX European birth cohort.

## Overview

This project applies a Dirichlet Process Mixture model to 1,301 children from six European birth cohorts (UK, France, Spain, Lithuania, Norway, Greece) to identify latent subgroups defined by combinations of 20 environmental exposures and estimate subgroup-specific obesity risk.

Implementation via the R package [PReMiuM](https://cran.r-project.org/package=PReMiuM) (Liverani et al., 2015).

## Covariates

20 categorical variables from 5 exposome families (4 per family):

| Family | Variables |
|--------|-----------|
| Air Pollution | NO₂, PM₂.₅, PM₁₀, PM absorbance |
| Metals | Pb, Hg, Cd, As |
| POPs | DDE, HCB, PCB153, PFOS |
| Lifestyle | Alcohol (pregnancy), Fish, Fruit, Physical activity |
| Built Environment | Population density, NDVI, Green space (300m), Traffic noise |

Continuous variables discretized into tertiles (Low/Medium/High). Outcome: overweight/obesity (binary).

## Pipeline

1. **Data preparation** — Variable selection, discretization, outcome coding
2. **MCMC sampling** — PReMiuM with nSweeps=20,000, nBurn=30,000, nFilter=2
3. **Convergence diagnostics** — Geweke test, ESS, trace plots, ACF
4. **Post-processing** — Dissimilarity matrix, optimal clustering, risk profiles
5. **Risk analysis** — Cluster-specific probabilities, Odds Ratios, PAF
6. **Profile characterization** — Dirichlet φ estimation, heatmaps, cluster comparisons
7. **Assumption verification** — Within-cluster LCA, outcome-covariate independence tests
8. **Sensitivity analysis** — Adjusted model (cohort + sex as fixed effects), cohort confounding assessment

## Key Results

- 10 clusters identified, 6 with significantly elevated risk (OR 2.20–4.95)
- High-risk profiles characterized by elevated air pollution, traffic noise, and low vegetation
- Clusters are >90% mono-cohort, reflecting geographic structuring of exposures
- Conditional independence assumptions satisfied (LCA: all clusters best K=1; chi-square: 3.4% significant)
- Adjusted analysis shows 3.5x credible interval inflation due to cohort-cluster quasi-collinearity

## Requirements

```r
install.packages(c("PReMiuM", "coda", "poLCA", "ggplot2", "dplyr", "tidyr"))
```

R ≥ 4.1.0

## Data

The HELIX dataset is not publicly available. Access can be requested through the [HELIX project](https://www.projecthelix.eu/).

## References

- Molitor et al. (2010). Bayesian profile regression. *Biostatistics*, 11(3), 484–498.
- Liverani et al. (2015). PReMiuM: An R package for profile regression mixture models. *Journal of Statistical Software*, 64(7).
- Malsiner-Walli et al. (2025). Without pain – clustering categorical data using Bayesian MFM of LCA. *Advances in Data Analysis and Classification*.
