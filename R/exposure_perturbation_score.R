#' Score exposures by network perturbation (leave-one-exposure-out)
#'
#' @description
#' Builds a reference network using all exposures, then for each exposure j builds
#' a network leaving out exposure j, and computes a perturbation score based on
#' differences in edge weights for a subset of tested pairs.
#'
#' This is an MVP perturbation metric designed to be interpretable and fast enough
#' for simulated/demo datasets.
#'
#' @param X Microbiome matrix (samples x microbes).
#' @param Y Metabolome matrix (samples x metabolites).
#' @param E Exposures matrix (samples x exposures).
#' @param covar Optional covariates data.frame.
#' @param fdr FDR threshold passed to build_exposure_network().
#' @param max_pairs Number of pairs to test per network build (speed control).
#' @param seed Random seed.
#'
#' @return A data.frame with exposure, perturbation_score, n_edges_ref, n_edges_drop.
#' 
#' @examples
#' set.seed(2)
#' d <- generate_dummy_exporisk(n = 30, p_micro = 10, p_metab = 12, p_expo = 4)
#' al <- align_omics(d$microbiome, d$metabolome, d$exposures, d$meta,
#'                  id_col = "sample_id", strict = TRUE)
#' pr <- prep_omics(al$microbiome, al$metabolome, al$exposures)
#' scores <- exposure_perturbation_score(pr$X, pr$Y, pr$E, fdr = 0.5, max_pairs = 120, seed = 1)
#' scores
#' @export
exposure_perturbation_score <- function(X, Y, E, covar = NULL,
                                        fdr = 0.2, max_pairs = 3000, seed = 1) {
  E <- if (is.data.frame(E)) as.matrix(E) else E
  if (is.null(colnames(E))) colnames(E) <- paste0("expo_", seq_len(ncol(E)))
  
  ref <- build_exposure_network(X, Y, E, covar = covar, fdr = fdr, max_pairs = max_pairs, seed = seed)
  
  # Use all ref edges (already subset by FDR)
  ref_edges <- ref$edges
  if (nrow(ref_edges) == 0) {
    out <- data.frame(
      exposure = colnames(E),
      perturbation_score = 0,
      n_edges_ref = 0,
      n_edges_drop = 0,
      stringsAsFactors = FALSE
    )
    return(out)
  }
  
  # Key for matching edges
  ref_key <- paste(ref_edges$microbe, ref_edges$metabolite, sep = "||")
  ref_w <- ref_edges$weight
  names(ref_w) <- ref_key
  
  res <- vector("list", ncol(E))
  
  for (j in seq_len(ncol(E))) {
    Ej <- E[, -j, drop = FALSE]
    drop <- build_exposure_network(X, Y, Ej, covar = covar, fdr = fdr, max_pairs = max_pairs, seed = seed)
    
    drop_edges <- drop$edges
    if (nrow(drop_edges) == 0) {
      score <- sum(abs(ref_w))  # maximal change (all ref edges "gone")
      res[[j]] <- data.frame(
        exposure = colnames(E)[j],
        perturbation_score = score,
        n_edges_ref = nrow(ref_edges),
        n_edges_drop = 0,
        stringsAsFactors = FALSE
      )
      next
    }
    
    drop_key <- paste(drop_edges$microbe, drop_edges$metabolite, sep = "||")
    drop_w <- drop_edges$weight
    names(drop_w) <- drop_key
    
    # Compare only common edges; penalize lost edges by abs(ref)
    common <- intersect(names(ref_w), names(drop_w))
    lost <- setdiff(names(ref_w), names(drop_w))
    
    score_common <- if (length(common) > 0) sum(abs(ref_w[common] - drop_w[common])) else 0
    score_lost <- if (length(lost) > 0) sum(abs(ref_w[lost])) else 0
    
    res[[j]] <- data.frame(
      exposure = colnames(E)[j],
      perturbation_score = score_common + score_lost,
      n_edges_ref = nrow(ref_edges),
      n_edges_drop = nrow(drop_edges),
      stringsAsFactors = FALSE
    )
  }
  
  out <- do.call(rbind, res)
  out <- out[order(out$perturbation_score, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}
