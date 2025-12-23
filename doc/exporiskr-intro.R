## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  fig.width = 8,
  fig.height = 5,
  dpi = 96
)

# Use a single seed for reproducibility across chunks.
set.seed(1)

## -----------------------------------------------------------------------------
library(ExpoRiskR)

## -----------------------------------------------------------------------------
d <- generate_dummy_exporisk(seed = 1)
str(d, max.level = 1)

## -----------------------------------------------------------------------------
al <- align_omics(
  microbiome = d$microbiome,
  metabolome = d$metabolome,
  exposures  = d$exposures,
  meta       = d$meta,
  id_col     = "sample_id",
  strict     = TRUE
)

## -----------------------------------------------------------------------------
pr <- prep_omics(al$microbiome, al$metabolome, al$exposures)

X <- pr$X
Y <- pr$Y
E <- pr$E

outcome <- d$meta$outcome
names(outcome) <- d$meta$sample_id

c(
  n_samples = nrow(X),
  n_microbes = ncol(X),
  n_metabolites = ncol(Y),
  n_exposures = ncol(E)
)

## -----------------------------------------------------------------------------
net <- build_exposure_network(
  X = X, Y = Y, E = E,
  fdr = 0.8,
  max_pairs = 1500,
  seed = 1
)

nrow(net$edges)

## ----fig.cap="Figure 1: Exposure-adjusted microbe--metabolite network. Nodes represent microbes and metabolites (bipartite graph). Edges represent microbe effects on metabolites after adjusting for exposures; edge width scales with |beta| and color indicates sign."----
plot_exposure_network(net)

## ----fig.cap="Figure 2: Exposure ranking. Exposures are ordered by their aggregate perturbation score computed from exposure-adjusted microbe--metabolite associations in the simulated dataset."----
scores <- exposure_perturbation_score(
  X = X, Y = Y, E = E,
  fdr = 0.8,
  max_pairs = 1500,
  seed = 1
)

print(plot_exposure_ranking(scores, top_n = 15))

## ----fig.cap="Figure 3: Network stability under bootstrap resampling. For vignette speed, a small number of bootstrap replicates is used; increase `n_boot` for more stable estimates in real analyses."----
print(plot_network_stability(
  X = X, Y = Y, E = E,
  n_boot = 8,
  fdr = 0.8,
  max_pairs = 1000,
  seed = 1
))

## ----fig.cap="Figure 4: Risk ROC curve for a simple risk model derived from exposure-adjusted associations. This is a demonstration on simulated data; results should not be interpreted as biological findings."----
print(plot_risk_roc(
  X = X, Y = Y, E = E,
  outcome = outcome,
  edges = net$edges,
  top_edges = 80
))

## ----fig.cap="Figure 5: Exposure feature importance for the simulated outcome. This is a lightweight illustration on simulated data; increase sample size and use domain-specific covariates for real analyses."----
print(plot_feature_importance(
  E = E,
  outcome = outcome,
  top_n = 15
))

## ----fig.cap="Figure 6: Individual risk profile for a single sample (based on exposures). The plot highlights the top exposure features contributing to the sample's risk score in the simulated dataset."----
sid <- rownames(E)[1]
print(plot_individual_risk_profile(
  sample_id = sid,
  E = E,
  outcome = outcome,
  top_n = 12
))

## -----------------------------------------------------------------------------
sessionInfo()

