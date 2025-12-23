# plotting_risk.R
# ------------------------------------------------------------

# ---- internal helpers (not exported) ----

.roc_curve <- function(y, score) {
  o <- order(score, decreasing = TRUE)
  y <- y[o]
  score <- score[o]
  
  P <- sum(y == 1)
  N <- sum(y == 0)
  if (P == 0 || N == 0) {
    stop("Need both classes (0/1) in outcome for ROC.", call. = FALSE)
  }
  
  tpr <- c(0, cumsum(y == 1) / P, 1)
  fpr <- c(0, cumsum(y == 0) / N, 1)
  data.frame(fpr = fpr, tpr = tpr)
}

.auc_trapz <- function(df) {
  x <- df$fpr
  y <- df$tpr
  sum((x[-1] - x[-length(x)]) * (y[-1] + y[-length(y)]) / 2)
}

.network_feature_per_sample <- function(X, Y, edges, top_edges = 200) {
  if (is.null(edges) || !is.data.frame(edges) || nrow(edges) == 0) {
    return(rep(0, nrow(X)))
  }
  
  req <- c("microbe", "metabolite", "weight")
  if (!all(req %in% colnames(edges))) {
    stop("`edges` must contain: microbe, metabolite, weight.", call. = FALSE)
  }
  
  ed <- edges[order(abs(edges$weight), decreasing = TRUE), , drop = FALSE]
  ed <- utils::head(ed, min(top_edges, nrow(ed)))
  ed <- ed[
    ed$microbe %in% colnames(X) &
      ed$metabolite %in% colnames(Y),
    ,
    drop = FALSE
  ]
  
  if (nrow(ed) == 0) return(rep(0, nrow(X)))
  
  s <- numeric(nrow(X))
  for (i in seq_len(nrow(ed))) {
    s <- s + (X[, ed$microbe[i]] *
                Y[, ed$metabolite[i]] *
                sign(ed$weight[i]))
  }
  
  as.vector(scale(s))

}

.fit_glm_prob <- function(df, outcome_col = "outcome") {
  f <- stats::as.formula(paste(outcome_col, "~ ."))
  fit <- stats::glm(f, data = df, family = stats::binomial())
  p <- stats::predict(fit, type = "response")
  list(fit = fit, prob = p)
}

# ------------------------------------------------------------
#' Plot disease risk stratification ROC curves (MVP)
#'
#' @param X Microbiome matrix (samples x features)
#' @param Y Metabolome matrix (samples x features)
#' @param E Exposures matrix (samples x features)
#' @param outcome Binary vector (0/1)
#' @param edges Network edges data.frame
#' @param top_edges Number of strongest edges for network feature
#' @examples
#' d <- generate_dummy_exporisk(seed = 1, n = 25, p_micro = 8, p_metab = 10, p_expo = 4)
#' al <- align_omics(d$microbiome, d$metabolome, d$exposures, d$meta,
#'                  id_col = "sample_id", strict = TRUE)
#' pr <- prep_omics(al$microbiome, al$metabolome, al$exposures)
#' net <- build_exposure_network(pr$X, pr$Y, pr$E, fdr = 0.95, max_pairs = 150, seed = 1)
#' outcome <- d$meta$outcome
#' names(outcome) <- d$meta$sample_id
#' p <- plot_risk_roc(pr$X, pr$Y, pr$E, outcome = outcome, edges = net$edges, top_edges = 30)
#' print(p)
#' @return A ggplot object
#'
#' @export
plot_risk_roc <- function(X, Y, E, outcome, edges, top_edges = 200) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_risk_roc().", call. = FALSE)
  }
  
  X <- .as_numeric_matrix(X, "X")
  Y <- .as_numeric_matrix(Y, "Y")
  E <- .as_numeric_matrix(E, "E")
  .check_same_rownames(X, Y, E)
  
  y <- as.numeric(outcome)
  if (length(y) != nrow(X)) {
    stop("`outcome` length must equal nrow(X).", call. = FALSE)
  }
  if (!all(y %in% c(0, 1))) {
    stop("`outcome` must be binary 0/1.", call. = FALSE)
  }
  
  net_feat <- .network_feature_per_sample(X, Y, edges, top_edges)
  
  dfE   <- data.frame(outcome = y, E)
  dfO   <- data.frame(outcome = y, X, Y)
  dfN   <- data.frame(outcome = y, net_feat = net_feat)
  dfAll <- data.frame(outcome = y, E, X, Y, net_feat = net_feat)
  
  mE   <- .fit_glm_prob(dfE)
  mO   <- .fit_glm_prob(dfO)
  mN   <- .fit_glm_prob(dfN)
  mAll <- .fit_glm_prob(dfAll)
  
  rocE   <- .roc_curve(y, mE$prob);   aucE   <- .auc_trapz(rocE)
  rocO   <- .roc_curve(y, mO$prob);   aucO   <- .auc_trapz(rocO)
  rocN   <- .roc_curve(y, mN$prob);   aucN   <- .auc_trapz(rocN)
  rocAll <- .roc_curve(y, mAll$prob); aucAll <- .auc_trapz(rocAll)
  
  rocE$model   <- sprintf("Exposures only (AUC=%.2f)", aucE)
  rocO$model   <- sprintf("Omics only (AUC=%.2f)", aucO)
  rocN$model   <- sprintf("Network feat (AUC=%.2f)", aucN)
  rocAll$model <- sprintf("Combined (AUC=%.2f)", aucAll)
  
  df <- rbind(rocE, rocO, rocN, rocAll)
  
  ggplot2::ggplot(df, ggplot2::aes(x = fpr, y = tpr, color = model)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "False Positive Rate",
      y = "True Positive Rate",
      title = "Disease risk stratification (ROC)",
      color = "Model"
    )
}

