# prep_omics.R
# ------------------------------------------------------------
# BiocCheck-friendly preprocessing for ExpoRiskR
# - No duplicated helper names across files
# - Uses internal dot helpers from utils.R:
#     .as_numeric_matrix(), .check_same_rownames()
# - Avoids stats::scale() (not exported). Uses base::scale().
# ------------------------------------------------------------

#' Preprocess exposures and multi-omics blocks for modeling
#'
#' @description
#' Lightweight preprocessing for MVP and Bioconductor-friendly workflows.
#' Converts inputs to numeric matrices, checks sample alignment, optionally
#' imputes missing values, applies log1p transforms, and scales features.
#'
#' @param microbiome Matrix/data.frame of samples x microbes.
#' @param metabolome Matrix/data.frame of samples x metabolites.
#' @param exposures  Matrix/data.frame of samples x exposures.
#' @param log1p_micro If TRUE (default), apply log1p to microbiome.
#' @param log1p_metab If TRUE (default), apply log1p to metabolome.
#' @param z_expo If TRUE (default), z-score exposures.
#' @param scale_omics If TRUE (default), center/scale microbiome and metabolome features.
#' @param na_action What to do with NA values: "error" (default) or "impute".
#'
#' @return A list with processed matrices: \code{X}, \code{Y}, \code{E}.
#'
#' @examples
#' set.seed(1)
#' d <- generate_dummy_exporisk(n = 20, p_micro = 6, p_metab = 8, p_expo = 3)
#' al <- align_omics(d$microbiome, d$metabolome, d$exposures, d$meta,
#'                  id_col = "sample_id", strict = TRUE)
#' pr <- prep_omics(al$microbiome, al$metabolome, al$exposures)
#' str(pr)
#'
#' @export
prep_omics <- function(microbiome, metabolome, exposures,
                       log1p_micro = TRUE,
                       log1p_metab = TRUE,
                       z_expo = TRUE,
                       scale_omics = TRUE,
                       na_action = c("error", "impute")) {
  
  na_action <- match.arg(na_action)
  
  # Coerce + validate
  X <- .as_numeric_matrix(microbiome, "microbiome")
  Y <- .as_numeric_matrix(metabolome, "metabolome")
  E <- .as_numeric_matrix(exposures,  "exposures")
  
  .check_same_rownames(X, Y, E)
  
  # Handle missing values
  if (na_action == "impute") {
    X <- .impute_col_means(X)
    Y <- .impute_col_means(Y)
    E <- .impute_col_means(E)
  } else {
    if (anyNA(X) || anyNA(Y) || anyNA(E)) {
      stop("Missing values detected; use `na_action = 'impute'` or remove NAs.", call. = FALSE)
    }
  }
  
  # Optional transforms
  if (log1p_micro) X <- log1p(X)
  if (log1p_metab) Y <- log1p(Y)
  
  # Optional scaling
  if (scale_omics) {
    X <- .safe_scale(X)
    Y <- .safe_scale(Y)
  }
  if (z_expo) {
    E <- .safe_scale(E)
  }
  
  list(X = X, Y = Y, E = E)
}

# ---- internal helpers (not exported) ----

.impute_col_means <- function(m) {
  if (!anyNA(m)) return(m)
  
  for (j in seq_len(ncol(m))) {
    v <- m[, j]
    if (anyNA(v)) {
      mu <- mean(v, na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      v[is.na(v)] <- mu
      m[, j] <- v
    }
  }
  m
}

# IMPORTANT:
# - Use base::scale(), NOT stats::scale() (stats does not export scale()).
# - Keep deterministic behavior; no randomness here.
.safe_scale <- function(m) {
  # SD per column; if SD==0 then only center (avoid division by zero)
  s <- apply(m, 2, stats::sd)
  keep <- is.finite(s) & s > 0
  
  out <- m
  
  if (any(keep)) {
    out[, keep] <- base::scale(m[, keep, drop = FALSE])
  }
  
  if (any(!keep)) {
    out[, !keep] <- m[, !keep, drop = FALSE] -
      apply(m[, !keep, drop = FALSE], 2, mean)
  }
  
  out
}
