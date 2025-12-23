#' Plot individual risk profile from exposure model
#'
#' @description
#' Fits outcome ~ exposures and shows per-exposure contribution for one sample
#' based on standardized coefficients and standardized exposure values.
#'
#' @param sample_id Sample ID (must be in rownames(E)).
#' @param E Numeric matrix (samples x exposures) with rownames.
#' @param outcome Binary vector (0/1), named by sample IDs or same row order as E.
#' @param top_n Number of top contributing exposures to display.
#' 
#' @examples
#' d <- generate_dummy_exporisk(seed = 1, n = 20, p_micro = 6, p_metab = 8, p_expo = 4)
#' outcome <- d$meta$outcome
#' names(outcome) <- d$meta$sample_id
#' sid <- rownames(d$exposures)[1]
#' p <- plot_individual_risk_profile(sample_id = sid, E = d$exposures, outcome = outcome, top_n = 10)
#' print(p)

#' @return A ggplot object.
#'
#' @export
plot_individual_risk_profile <- function(sample_id, E, outcome, top_n = 20) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_individual_risk_profile().", call. = FALSE)
  }
  
  E <- .as_numeric_matrix(E, "E")
  if (is.null(rownames(E))) stop("E must have rownames (sample IDs).", call. = FALSE)
  if (!sample_id %in% rownames(E)) stop("sample_id not found in rownames(E).", call. = FALSE)
  
  # align outcome
  y <- outcome
  if (!is.null(names(y))) {
    if (!all(rownames(E) %in% names(y))) stop("Named outcome must cover all E rownames.", call. = FALSE)
    y <- as.numeric(y[rownames(E)])
  } else {
    y <- as.numeric(y)
    if (length(y) != nrow(E)) stop("Un-named outcome length must equal nrow(E).", call. = FALSE)
  }
  if (!all(y %in% c(0, 1))) stop("`outcome` must be binary 0/1.", call. = FALSE)
  
  coln <- colnames(E)
  if (is.null(coln)) coln <- paste0("expo_", seq_len(ncol(E)))
  
  Ez <- scale(E)
  df <- data.frame(outcome = y, Ez)
  fit <- stats::glm(outcome ~ ., data = df, family = stats::binomial())
  
  coefs <- stats::coef(fit)
  b0 <- coefs["(Intercept)"]
  b <- coefs[names(coefs) != "(Intercept)"]
  
  # sample standardized values
  z <- Ez[rownames(E) == sample_id, ]
  z <- as.numeric(z)
  names(z) <- sub("^X", "", names(b))
  
  # contributions in log-odds space
  contrib <- b * z
  contrib <- contrib[is.finite(contrib)]
  contrib <- sort(contrib, decreasing = TRUE)
  
  # predicted probability
  eta <- b0 + sum(contrib, na.rm = TRUE)
  prob <- 1 / (1 + exp(-eta))
  
  out <- data.frame(
    feature = names(contrib),
    contribution = as.numeric(contrib),
    stringsAsFactors = FALSE
  )
  
  # show top positive and top negative
  top_pos <- utils::head(out[order(out$contribution, decreasing = TRUE), , drop = FALSE], top_n)
  top_neg <- utils::head(out[order(out$contribution, decreasing = FALSE), , drop = FALSE], top_n)
  out2 <- rbind(top_pos, top_neg)
  out2 <- out2[!duplicated(out2$feature), , drop = FALSE]
  
  ggplot2::ggplot(out2, ggplot2::aes(x = stats::reorder(feature, contribution), y = contribution)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Contribution (log-odds)",
      title = sprintf("Individual risk profile: %s (pred=%.2f)", sample_id, prob)
    )
}
