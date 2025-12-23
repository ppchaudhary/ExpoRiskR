# align_omics.R
# ----------------------------

#' Align exposures and multi-omics blocks by sample ID
#'
#' @description
#' Ensures that microbiome, metabolome, exposures, and metadata all refer to the same
#' set of samples in the same order. Sample IDs are taken from rownames of matrices/
#' data.frames, or from a column in `meta` if `id_col` is provided.
#'
#' @param microbiome Matrix/data.frame of samples x microbes.
#' @param metabolome Matrix/data.frame of samples x metabolites.
#' @param exposures  Matrix/data.frame of samples x exposures.
#' @param meta       data.frame of sample-level metadata (must include outcome later).
#' @param id_col     Optional column name in `meta` containing sample IDs.
#'                  If NULL, rownames(meta) are used (if present).
#' @param strict     If TRUE, errors if any block has samples not found in others.
#'                  If FALSE, intersects common samples and drops others.
#'
#' @return A list with aligned `microbiome`, `metabolome`, `exposures`, `meta`, and `sample_id`.
#'
#' @examples
#' set.seed(4)
#' d <- generate_dummy_exporisk(n = 20, p_micro = 6, p_metab = 8, p_expo = 3)
#' aligned <- align_omics(d$microbiome, d$metabolome, d$exposures, d$meta,
#'                        id_col = "sample_id", strict = TRUE)
#' names(aligned)
#'
#' @export
align_omics <- function(microbiome, metabolome, exposures, meta,
                        id_col = NULL, strict = TRUE) {
  
  X <- .as_matrix(microbiome, "microbiome")
  Y <- .as_matrix(metabolome, "metabolome")
  E <- .as_matrix(exposures, "exposures")
  
  if (!is.data.frame(meta)) {
    stop("`meta` must be a data.frame.", call. = FALSE)
  }
  
  ids_meta <- .get_meta_ids(meta, id_col)
  
  ids_X <- rownames(X)
  ids_Y <- rownames(Y)
  ids_E <- rownames(E)
  
  if (is.null(ids_X) || anyNA(ids_X) || any(ids_X == "")) {
    stop("`microbiome` must have non-empty rownames as sample IDs.", call. = FALSE)
  }
  if (is.null(ids_Y) || anyNA(ids_Y) || any(ids_Y == "")) {
    stop("`metabolome` must have non-empty rownames as sample IDs.", call. = FALSE)
  }
  if (is.null(ids_E) || anyNA(ids_E) || any(ids_E == "")) {
    stop("`exposures` must have non-empty rownames as sample IDs.", call. = FALSE)
  }
  
  if (any(duplicated(ids_X))) stop("Duplicate sample IDs found in microbiome rownames.", call. = FALSE)
  if (any(duplicated(ids_Y))) stop("Duplicate sample IDs found in metabolome rownames.", call. = FALSE)
  if (any(duplicated(ids_E))) stop("Duplicate sample IDs found in exposures rownames.", call. = FALSE)
  if (any(duplicated(ids_meta))) stop("Duplicate sample IDs found in meta.", call. = FALSE)
  
  all_sets <- list(microbiome = ids_X, metabolome = ids_Y, exposures = ids_E, meta = ids_meta)
  
  if (strict) {
    ref <- sort(ids_meta)
    for (nm in names(all_sets)) {
      if (!identical(sort(all_sets[[nm]]), ref)) {
        missing_in_block <- setdiff(ref, all_sets[[nm]])
        extra_in_block <- setdiff(all_sets[[nm]], ref)
        msg <- paste0(
          "Sample ID mismatch between meta and ", nm, ".\n",
          if (length(missing_in_block) > 0) paste0("Missing in ", nm, ": ", paste(head(missing_in_block, 10), collapse = ", "),
                                                   if (length(missing_in_block) > 10) " ..." else "", "\n") else "",
          if (length(extra_in_block) > 0) paste0("Extra in ", nm, ": ", paste(head(extra_in_block, 10), collapse = ", "),
                                                 if (length(extra_in_block) > 10) " ..." else "", "\n") else ""
        )
        stop(msg, call. = FALSE)
      }
    }
    keep_ids <- ids_meta
  } else {
    keep_ids <- Reduce(intersect, all_sets)
    if (length(keep_ids) < 2) {
      stop("After intersecting, fewer than 2 common samples remain.", call. = FALSE)
    }
  }
  
  # Keep order consistent with meta
  keep_ids <- keep_ids[keep_ids %in% ids_meta]
  
  ord_X <- match(keep_ids, ids_X)
  ord_Y <- match(keep_ids, ids_Y)
  ord_E <- match(keep_ids, ids_E)
  ord_M <- match(keep_ids, ids_meta)
  
  X2 <- X[ord_X, , drop = FALSE]
  Y2 <- Y[ord_Y, , drop = FALSE]
  E2 <- E[ord_E, , drop = FALSE]
  meta2 <- meta[ord_M, , drop = FALSE]
  
  rownames(X2) <- rownames(Y2) <- rownames(E2) <- keep_ids
  if (!is.null(id_col)) {
    meta2[[id_col]] <- keep_ids
    rownames(meta2) <- keep_ids
  } else {
    rownames(meta2) <- keep_ids
  }
  
  list(
    microbiome = X2,
    metabolome = Y2,
    exposures  = E2,
    meta       = meta2,
    sample_id  = keep_ids
  )
}

# ---- internal helpers (not exported) ----
# IMPORTANT: keep this helper defined in ONE place in the package
# to avoid silent overwrites across R/ files.
.as_matrix <- function(x, nm) {
  if (is.matrix(x)) return(x)
  if (is.data.frame(x)) return(as.matrix(x))
  stop("`", nm, "` must be a matrix or data.frame.", call. = FALSE)
}

.get_meta_ids <- function(meta, id_col) {
  if (!is.null(id_col)) {
    if (!is.character(id_col) || length(id_col) != 1) {
      stop("`id_col` must be a single character string.", call. = FALSE)
    }
    if (!id_col %in% names(meta)) {
      stop("`id_col` not found in meta: ", id_col, call. = FALSE)
    }
    ids <- meta[[id_col]]
  } else {
    ids <- rownames(meta)
    if (is.null(ids)) {
      stop("`meta` must have rownames as sample IDs, or provide `id_col`.", call. = FALSE)
    }
  }
  
  if (anyNA(ids) || any(ids == "")) {
    stop("Sample IDs in meta are missing/empty.", call. = FALSE)
  }
  as.character(ids)
}
