#' Plot feature importance for exposures (logistic regression)
#'
#' @description
#' Fits a logistic regression outcome ~ exposures and ranks exposures by the
#' absolute standardized coefficient magnitude.
#'
#' @param E Numeric matrix (samples x exposures).
#' @param outcome Binary vector (0/1), length = nrow(E).
#' @param top_n Number of top exposures to show.
#' @examples
#' d <- generate_dummy_exporisk(seed = 1, n = 20, p_micro = 6, p_metab = 8, p_expo = 4)
#' outcome <- d$meta$outcome
#' names(outcome) <- d$meta$sample_id
#' p <- plot_feature_importance(E = d$exposures, outcome = outcome, top_n = 10)
#' print(p)

#' @return A ggplot object.
#'
#' @export
plot_feature_importance <- function(E, outcome, top_n = 25) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_feature_importance().", call. = FALSE)
  }
  
  E <- .as_numeric_matrix(E, "E")
  y <- as.numeric(outcome)
  if (length(y) != nrow(E)) stop("`outcome` length must equal nrow(E).", call. = FALSE)
  if (!all(y %in% c(0, 1))) stop("`outcome` must be binary 0/1.", call. = FALSE)
  
  coln <- colnames(E)
  if (is.null(coln)) coln <- paste0("expo_", seq_len(ncol(E)))
  
  # standardize predictors for comparable coefficients
  Ez <- scale(E)
  df <- data.frame(outcome = y, Ez)
  
  fit <- stats::glm(outcome ~ ., data = df, family = stats::binomial())
  coefs <- stats::coef(fit)
  
  # drop intercept
  coefs <- coefs[names(coefs) != "(Intercept)"]
  imp <- abs(coefs)
  nm <- names(imp)
  # nice names (remove leading "X" from data.frame conversions if any)
  nm <- sub("^X", "", nm)
  
  out <- data.frame(
    feature = nm,
    importance = as.numeric(imp),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$importance, decreasing = TRUE), , drop = FALSE]
  out <- utils::head(out, top_n)
  
  ggplot2::ggplot(out, ggplot2::aes(x = stats::reorder(feature, importance), y = importance)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Abs. standardized coefficient",
      title = "Feature importance (exposures)"
    )
}
