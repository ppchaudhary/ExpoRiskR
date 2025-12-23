# Run all ExpoRiskR plots interactively (RStudio Plots pane)

# ---- setup ----
if (!requireNamespace("ExpoRiskR", quietly = TRUE)) {
  devtools::load_all(".")
}

library(ExpoRiskR)

# ---- reproducible data ----
d <- generate_dummy_exporisk(seed = 1)

al <- align_omics(
  d$microbiome,
  d$metabolome,
  d$exposures,
  d$meta,
  id_col = "sample_id",
  strict = TRUE
)

pr <- prep_omics(al$microbiome, al$metabolome, al$exposures)
X <- pr$X
Y <- pr$Y
E <- pr$E
outcome <- d$meta$outcome
names(outcome) <- d$meta$sample_id

# ---- build network (relaxed FDR so edges exist) ----
net <- build_exposure_network(
  X, Y, E,
  fdr = 0.8,
  max_pairs = 5000,
  seed = 1
)

# ========== PLOT 1 ==========
# Network plot
plot_exposure_network(net)
message("Plot 1: Exposure network")
readline("Press <Enter> for next plot")

# ========== PLOT 2 ==========
# Exposure ranking
scores <- exposure_perturbation_score(
  X, Y, E,
  fdr = 0.8,
  max_pairs = 5000,
  seed = 1
)
p_rank <- plot_exposure_ranking(scores, top_n = 20)
print(p_rank)
message("Plot 2: Exposure ranking")
readline("Press <Enter> for next plot")

# ========== PLOT 3 ==========
# Network stability
p_stab <- plot_network_stability(
  X, Y, E,
  n_boot = 30,
  fdr = 0.8,
  max_pairs = 3000,
  seed = 1
)
print(p_stab)
message("Plot 3: Network stability")
readline("Press <Enter> for next plot")

# ========== PLOT 4 ==========
# Risk ROC
p_roc <- plot_risk_roc(
  X, Y, E,
  outcome = outcome,
  edges = net$edges,
  top_edges = 100
)
print(p_roc)
message("Plot 4: Risk ROC")
readline("Press <Enter> for next plot")

# ========== PLOT 5 ==========
# Feature importance (exposures only)
p_imp <- plot_feature_importance(E = E, outcome = outcome, top_n = 20)
print(p_imp)
message("Plot 5: Feature importance")
readline("Press <Enter> for next plot")


# ========== PLOT 6 ==========
# Individual risk profile (exposures only)
sid <- rownames(E)[1]
p_ind <- plot_individual_risk_profile(
  sample_id = sid,
  E = E,
  outcome = outcome,
  top_n = 15
)
print(p_ind)
message("Plot 6: Individual risk profile")

message("✅ All plots shown successfully")

