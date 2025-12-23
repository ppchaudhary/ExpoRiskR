# align_omics_se.R
# ------------------------------------------------------------
# Create aligned SummarizedExperiment objects from (X, Y, E, meta)
# ------------------------------------------------------------
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
NULL

#' Align two SummarizedExperiment objects and extract exposures from colData
#'
#' @description
#' Convenience wrapper to (i) align microbiome, metabolome, and exposures by sample ID
#' and (ii) return two \code{SummarizedExperiment} objects (microbiome + metabolome)
#' that share the same \code{colData} (meta + exposures). This is useful for
#' Bioconductor-style workflows.
#'
#' Inputs \code{microbiome}, \code{metabolome}, \code{exposures} are expected to be
#' sample-by-feature matrices (or coercible to matrices). Sample IDs are taken from
#' rownames when present; otherwise from \code{meta[[id_col]]}.
#'
#' @param microbiome Matrix/data.frame (samples x microbes).
#' @param metabolome Matrix/data.frame (samples x metabolites).
#' @param exposures  Matrix/data.frame (samples x exposures).
#' @param meta Data.frame with sample metadata including \code{id_col}.
#' @param id_col Column name in \code{meta} holding sample IDs (default "sample_id").
#' @param strict If TRUE, require that all blocks contain the same sample IDs;
#'   otherwise subset to the intersection (default TRUE).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{se_microbiome}: SummarizedExperiment for microbiome (features x samples)
#'   \item \code{se_metabolome}: SummarizedExperiment for metabolome (features x samples)
#'   \item \code{exposures}: aligned numeric matrix (samples x exposures)
#'   \item \code{meta}: aligned meta data.frame
#'   \item \code{sample_ids}: character vector of aligned sample IDs
#' }
#'
#' @examples
#' set.seed(7)
#' d <- generate_dummy_exporisk(n = 12, p_micro = 5, p_metab = 6, p_expo = 3)
#' out <- align_omics_se(
#'   d$microbiome, d$metabolome, d$exposures, d$meta,
#'   id_col = "sample_id", strict = TRUE
#' )
#' out$se_microbiome
#' out$se_metabolome
#'
#' @export
align_omics_se <- function(microbiome, metabolome, exposures, meta,
                           id_col = "sample_id",
                           strict = TRUE) {
  
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required for align_omics_se().", call. = FALSE)
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required for align_omics_se().", call. = FALSE)
  }
  
  if (!is.data.frame(meta)) {
    stop("`meta` must be a data.frame.", call. = FALSE)
  }
  if (!id_col %in% colnames(meta)) {
    stop(sprintf("`meta` must contain column '%s'.", id_col), call. = FALSE)
  }
  
  X <- .as_numeric_matrix(microbiome, "microbiome")
  Y <- .as_numeric_matrix(metabolome, "metabolome")
  E <- .as_numeric_matrix(exposures,  "exposures")
  
  # Determine IDs for each block
  ids_meta <- as.character(meta[[id_col]])
  
  ids_X <- rownames(X)
  if (is.null(ids_X)) ids_X <- ids_meta
  
  ids_Y <- rownames(Y)
  if (is.null(ids_Y)) ids_Y <- ids_meta
  
  ids_E <- rownames(E)
  if (is.null(ids_E)) ids_E <- ids_meta
  
  # Safety: ensure meta IDs are unique and non-missing
  if (anyNA(ids_meta) || any(ids_meta == "")) {
    stop("`meta[[id_col]]` contains missing/empty sample IDs.", call. = FALSE)
  }
  if (any(duplicated(ids_meta))) {
    stop("`meta[[id_col]]` must contain unique sample IDs.", call. = FALSE)
  }
  
  # Align IDs
  if (isTRUE(strict)) {
    if (!setequal(ids_X, ids_meta) || !setequal(ids_Y, ids_meta) || !setequal(ids_E, ids_meta)) {
      stop("Sample IDs do not match across blocks (strict=TRUE).", call. = FALSE)
    }
    common <- sort(ids_meta)
  } else {
    common <- Reduce(intersect, list(ids_X, ids_Y, ids_E, ids_meta))
    common <- sort(unique(common))
    if (length(common) == 0) {
      stop("No overlapping sample IDs found across blocks.", call. = FALSE)
    }
  }
  
  # Helper to reorder/subset by ID
  .subset_by_id <- function(mat, ids_mat, target_ids) {
    # Map target_ids to positions in ids_mat
    idx <- match(target_ids, ids_mat)
    if (anyNA(idx)) {
      stop("Internal alignment error: missing IDs while subsetting.", call. = FALSE)
    }
    mat[idx, , drop = FALSE]
  }
  
  X2 <- .subset_by_id(X, ids_X, common)
  Y2 <- .subset_by_id(Y, ids_Y, common)
  E2 <- .subset_by_id(E, ids_E, common)
  
  meta2 <- meta[match(common, ids_meta), , drop = FALSE]
  
  # Prepare colData: meta + exposures (as DataFrame)
  E_df <- as.data.frame(E2)
  colData <- S4Vectors::DataFrame(meta2, E_df, row.names = common)
  
  # SummarizedExperiment assays are features x samples, so transpose:
  # input matrices are samples x features
  X_assay <- t(X2)
  Y_assay <- t(Y2)
  
  # Ensure colnames are sample IDs
  colnames(X_assay) <- common
  colnames(Y_assay) <- common
  
  se_micro <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = X_assay),
    colData = colData
  )
  
  se_metab <- SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = Y_assay),
    colData = colData
  )
  
  list(
    se_microbiome = se_micro,
    se_metabolome = se_metab,
    exposures = E2,
    meta = meta2,
    sample_ids = common
  )
}
