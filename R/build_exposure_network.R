# build_exposure_network.R

#' Build an exposure-adjusted microbe-metabolite association network
#'
#' @description
#' For each (microbe, metabolite) pair, fits a linear model:
#' \deqn{metabolite ~ microbe + exposures + covariates}
#' and uses the microbe coefficient as the edge weight.
#'
#' This is an MVP, interpretable approach suitable for Bioconductor submission.
#'
#' @param X Numeric matrix (samples x microbes).
#' @param Y Numeric matrix (samples x metabolites).
#' @param E Numeric matrix (samples x exposures).
#' @param covar Optional data.frame of sample-level covariates (rows = samples).
#' @param fdr FDR threshold for keeping edges (BH adjusted p-value).
#' @param max_pairs Max number of (microbe, metabolite) pairs to test (for speed).
#'   If NULL, tests all pairs (may be slow).
#' @param seed Optional random seed used only when `max_pairs` is not NULL and
#'   sampling is required. If `NULL`, the current RNG state is used.
#'
#' @return A list with:
#' \itemize{
#'   \item edges: data.frame of significant edges (microbe, metabolite, weight, p_value, fdr)
#'   \item graph: igraph object (bipartite)
#'   \item meta: list of settings and counts
#' }
#'
#' @examples
#' set.seed(1)
#' d <- generate_dummy_exporisk(n = 30, p_micro = 10, p_metab = 12, p_expo = 4)
#' al <- align_omics(d$microbiome, d$metabolome, d$exposures, d$meta,
#'                  id_col = "sample_id", strict = TRUE)
#' pr <- prep_omics(al$microbiome, al$metabolome, al$exposures)
#' net <- build_exposure_network(pr$X, pr$Y, pr$E, fdr = 0.5, max_pairs = 120, seed = 1)
#' utils::head(net$edges)
#'
#' @export
build_exposure_network <- function(X, Y, E, covar = NULL, fdr = 0.1,
                                   max_pairs = 5000, seed = NULL) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for build_exposure_network().", call. = FALSE)
  }
  
  X <- .as_numeric_matrix(X, "X")
  Y <- .as_numeric_matrix(Y, "Y")
  E <- .as_numeric_matrix(E, "E")
  
  .check_same_rownames(X, Y, E)
  
  if (!is.null(covar)) {
    if (!is.data.frame(covar)) stop("`covar` must be a data.frame.", call. = FALSE)
    if (nrow(covar) != nrow(X)) stop("`covar` must have same number of rows as X.", call. = FALSE)
  }
  
  p_micro <- ncol(X)
  p_metab <- ncol(Y)
  
  if (p_micro < 1 || p_metab < 1) {
    stop("X and Y must each have at least 1 feature.", call. = FALSE)
  }
  
  all_pairs <- expand.grid(m = seq_len(p_micro), k = seq_len(p_metab))
  n_all <- nrow(all_pairs)
  
  ## --- REPLACED BLOCK (no set.seed) ---
  if (!is.null(max_pairs) && max_pairs < n_all) {
    
    if (!is.null(seed)) {
      if (!requireNamespace("withr", quietly = TRUE)) {
        stop("Package 'withr' is required when `seed` is provided.", call. = FALSE)
      }
      idx <- withr::with_seed(
        seed,
        sample.int(n_all, size = max_pairs, replace = FALSE)
      )
    } else {
      idx <- sample.int(n_all, size = max_pairs, replace = FALSE)
    }
    
    pairs <- all_pairs[idx, , drop = FALSE]
    
  } else {
    pairs <- all_pairs
  }
  
  # Build design matrix for exposures + covariates (same across models)
  Z <- E
  if (!is.null(covar) && ncol(covar) > 0) {
    covar_m <- stats::model.matrix(~ ., data = covar)
    covar_m <- covar_m[, colnames(covar_m) != "(Intercept)", drop = FALSE]
    if (ncol(covar_m) > 0) Z <- cbind(Z, covar_m)
  }
  
  # Pre-allocate results
  out_micro <- character(nrow(pairs))
  out_metab <- character(nrow(pairs))
  out_beta  <- numeric(nrow(pairs))
  out_p     <- numeric(nrow(pairs))
  
  micro_names <- colnames(X)
  metab_names <- colnames(Y)
  
  if (is.null(micro_names)) micro_names <- paste0("micro_", seq_len(p_micro))
  if (is.null(metab_names)) metab_names <- paste0("metab_", seq_len(p_metab))
  
  for (i in seq_len(nrow(pairs))) {
    mi <- pairs$m[i]
    ki <- pairs$k[i]
    
    x <- X[, mi]
    y <- Y[, ki]
    
    df <- data.frame(y = y, x = x)
    
    if (!is.null(Z) && ncol(Z) > 0) {
      dfZ <- as.data.frame(Z)
      df <- cbind(df, dfZ)
    }
    
    # remove incomplete rows and skip if too few remain
    df <- df[stats::complete.cases(df), , drop = FALSE]
    if (nrow(df) < 3) {
      out_beta[i] <- NA_real_
      out_p[i] <- NA_real_
      out_micro[i] <- micro_names[mi]
      out_metab[i] <- metab_names[ki]
      next
    }
    
    fit <- stats::lm(y ~ ., data = df)
    sm <- summary(fit)$coefficients
    
    if (!"x" %in% rownames(sm)) {
      out_beta[i] <- NA_real_
      out_p[i] <- NA_real_
    } else {
      out_beta[i] <- sm["x", "Estimate"]
      out_p[i]    <- sm["x", "Pr(>|t|)"]
    }
    
    out_micro[i] <- micro_names[mi]
    out_metab[i] <- metab_names[ki]
  }
  
  edges <- data.frame(
    microbe = out_micro,
    metabolite = out_metab,
    weight = out_beta,
    p_value = out_p,
    stringsAsFactors = FALSE
  )
  
  edges <- edges[is.finite(edges$p_value) & is.finite(edges$weight), , drop = FALSE]
  
  if (nrow(edges) == 0) {
    g <- igraph::make_empty_graph()
    return(list(
      edges = edges,
      graph = g,
      meta = list(n_pairs_tested = nrow(pairs), n_edges = 0, fdr = fdr)
    ))
  }
  
  edges$fdr <- stats::p.adjust(edges$p_value, method = "BH")
  edges <- edges[edges$fdr <= fdr, , drop = FALSE]
  
  if (nrow(edges) == 0) {
    g <- igraph::make_empty_graph()
    return(list(
      edges = edges,
      graph = g,
      meta = list(n_pairs_tested = nrow(pairs), n_edges = 0, fdr = fdr)
    ))
  }
  
  # Build bipartite graph
  verts <- unique(c(edges$microbe, edges$metabolite))
  g <- igraph::graph_from_data_frame(
    d = data.frame(from = edges$microbe, to = edges$metabolite, weight = edges$weight),
    directed = FALSE,
    vertices = data.frame(name = verts, stringsAsFactors = FALSE)
  )
  
  # bipartite type: TRUE for metabolite nodes
  igraph::V(g)$type <- igraph::V(g)$name %in% unique(edges$metabolite)
  
  list(
    edges = edges,
    graph = g,
    meta = list(n_pairs_tested = nrow(pairs), n_edges = nrow(edges), fdr = fdr)
  )
}
