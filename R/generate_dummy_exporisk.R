# generate_dummy_exporisk.R
# ------------------------------------------------------------

#' Generate simulated exposure + multi-omics data with a binary outcome
#'
#' @description
#' Creates a reproducible toy dataset for demonstrating ExpoRiskR workflows:
#' exposures (E), microbiome-like positive features (X), metabolome-like positive
#' features (Y), and a binary disease outcome.
#'
#' If `seed` is provided, reproducibility is ensured locally without modifying
#' the global RNG state.
#'
#' @param n Number of samples.
#' @param p_micro Number of microbiome features.
#' @param p_metab Number of metabolomics features.
#' @param p_expo Number of exposure variables.
#' @param n_signal Number of truly associated features per block.
#' @param seed Optional random seed for reproducible simulation.
#'
#' @return A list with matrices: microbiome, metabolome, exposures; and meta data.frame.
#'
#' @examples
#' d <- generate_dummy_exporisk(n = 20, p_micro = 6, p_metab = 8, p_expo = 3, seed = 1)
#' str(d)
#'
#' @export
generate_dummy_exporisk <- function(n = 120,
                                    p_micro = 50,
                                    p_metab = 80,
                                    p_expo = 10,
                                    n_signal = 6,
                                    seed = NULL) {
  
  if (!is.null(seed)) {
    if (!requireNamespace("withr", quietly = TRUE)) {
      stop("Package 'withr' is required when `seed` is provided.", call. = FALSE)
    }
  }
  
  sim_fun <- function() {
    
    # Exposures (can be negative; that's fine)
    E <- matrix(stats::rnorm(n * p_expo), nrow = n, ncol = p_expo)
    colnames(E) <- paste0("expo_", seq_len(p_expo))
    
    # Latent exposure-driven factor
    z <- base::scale(
      E[, 1] * 0.7 +
        E[, 2] * 0.5 +
        stats::rnorm(n, sd = 0.7)
    )[, 1]
    
    # Positive omics blocks (log-normal)
    X <- matrix(
      exp(stats::rnorm(n * p_micro, mean = 0, sd = 1)),
      nrow = n, ncol = p_micro
    )
    colnames(X) <- paste0("micro_", seq_len(p_micro))
    
    Y <- matrix(
      exp(stats::rnorm(n * p_metab, mean = 0, sd = 1)),
      nrow = n, ncol = p_metab
    )
    colnames(Y) <- paste0("metab_", seq_len(p_metab))
    
    # Inject signal multiplicatively so data stays positive
    n_signal_use <- min(n_signal, p_micro, p_metab)
    mult_x <- exp(z * 0.35)
    mult_y <- exp(z * 0.40)
    
    if (n_signal_use > 0) {
      X[, seq_len(n_signal_use)] <-
        X[, seq_len(n_signal_use), drop = FALSE] * mult_x
      Y[, seq_len(n_signal_use)] <-
        Y[, seq_len(n_signal_use), drop = FALSE] * mult_y
    }
    
    # Outcome depends on latent factor + a couple of omics signals
    linpred <- 0.8 * z +
      0.15 * log1p(X[, 1]) +
      0.15 * log1p(Y[, 1]) +
      stats::rnorm(n, sd = 0.6)
    
    p <- 1 / (1 + exp(-linpred))
    outcome <- stats::rbinom(n, size = 1, prob = p)
    
    meta <- data.frame(
      sample_id = paste0("S", seq_len(n)),
      outcome = outcome
    )
    
    rownames(X) <- rownames(Y) <- rownames(E) <- meta$sample_id
    
    list(
      microbiome = X,
      metabolome = Y,
      exposures = E,
      meta = meta
    )
  }
  
  if (!is.null(seed)) {
    withr::with_seed(seed, sim_fun())
  } else {
    sim_fun()
  }
}
