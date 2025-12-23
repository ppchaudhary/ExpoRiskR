#' Plot network stability by bootstrap edge overlap
#'
#' @description
#' Builds a reference network using all samples, then repeatedly bootstraps
#' samples with replacement, rebuilds the network, and computes Jaccard overlap
#' between edge sets.
#'
#' @param X Numeric matrix (samples x microbes).
#' @param Y Numeric matrix (samples x metabolites).
#' @param E Numeric matrix (samples x exposures).
#' @param n_boot Number of bootstrap resamples.
#' @param fdr FDR threshold passed to build_exposure_network().
#' @param max_pairs Maximum pairs passed to build_exposure_network().
#' @param seed Optional seed controlling bootstrap resampling only.
#'
#' @examples
#' d <- generate_dummy_exporisk(seed = 1, n = 20, p_micro = 8, p_metab = 10, p_expo = 4)
#' al <- align_omics(d$microbiome, d$metabolome, d$exposures, d$meta,
#'                  id_col = "sample_id", strict = TRUE)
#' pr <- prep_omics(al$microbiome, al$metabolome, al$exposures)
#' p <- plot_network_stability(pr$X, pr$Y, pr$E, n_boot = 2, fdr = 0.95, max_pairs = 120, seed = 1)
#' print(p)
#' @return A ggplot object.
#'
#' @export
plot_network_stability <- function(X, Y, E,
                                   n_boot = 50,
                                   fdr = 0.2,
                                   max_pairs = 2000,
                                   seed = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_network_stability().", call. = FALSE)
  }
  if (!requireNamespace("withr", quietly = TRUE)) {
    stop("Package 'withr' is required for plot_network_stability().", call. = FALSE)
  }
  
  X <- .as_numeric_matrix(X, "X")
  Y <- .as_numeric_matrix(Y, "Y")
  E <- .as_numeric_matrix(E, "E")
  .check_same_rownames(X, Y, E)
  
  n <- nrow(X)
  if (n < 5) stop("Need at least 5 samples for bootstrap stability.", call. = FALSE)
  
  # reference network on full data
  ref <- build_exposure_network(X, Y, E, fdr = fdr, max_pairs = max_pairs, seed = seed)
  ref_edges <- ref$edges
  if (is.null(ref_edges) || nrow(ref_edges) == 0) {
    warning("Reference network has 0 edges; stability will be 0.", call. = FALSE)
    ref_set <- character(0)
  } else {
    ref_set <- paste(ref_edges$microbe, ref_edges$metabolite, sep = "||")
  }
  
  jaccard <- function(a, b) {
    a <- unique(a); b <- unique(b)
    if (length(a) == 0 && length(b) == 0) return(1)
    if (length(a) == 0 || length(b) == 0) return(0)
    length(intersect(a, b)) / length(union(a, b))
  }
  
  overlaps <- numeric(n_boot)
  
  boot_fun <- function() {
    idx <- sample.int(n, size = n, replace = TRUE)
    Xb <- X[idx, , drop = FALSE]
    Yb <- Y[idx, , drop = FALSE]
    Eb <- E[idx, , drop = FALSE]
    rownames(Xb) <- rownames(Yb) <- rownames(Eb) <- rownames(X)[idx]
    
    nb <- build_exposure_network(Xb, Yb, Eb, fdr = fdr, max_pairs = max_pairs, seed = NULL)
    if (is.null(nb$edges) || nrow(nb$edges) == 0) return(0)
    b_set <- paste(nb$edges$microbe, nb$edges$metabolite, sep = "||")
    jaccard(ref_set, b_set)
  }
  
  if (!is.null(seed)) {
    overlaps <- withr::with_seed(seed, replicate(n_boot, boot_fun()))
  } else {
    overlaps <- replicate(n_boot, boot_fun())
  }
  
  df <- data.frame(overlap = overlaps)
  
  ggplot2::ggplot(df, ggplot2::aes(x = overlap)) +
    ggplot2::geom_histogram(bins = 20) +
    ggplot2::labs(
      x = "Jaccard overlap vs reference",
      y = "Bootstrap count",
      title = "Network stability (bootstrap edge overlap)"
    )
}
