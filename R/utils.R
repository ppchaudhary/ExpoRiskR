# utils.R
# ================================
# Internal utility helper functions
# ================================

# NOTE (Bioconductor-friendly design):
# - The "dot" helpers (.as_numeric_matrix, .check_same_rownames) are internal API
#   used throughout the package.
# - The non-dot wrappers are kept for backward compatibility but are NOT exported.
#
# IMPORTANT:
# - Do NOT document internal helpers with roxygen blocks that create man pages.
#   This avoids generating man/dot-*.Rd files that trigger BiocCheck example rules.

.as_numeric_matrix <- function(x, name = "x") {
  if (is.null(x)) {
    stop(sprintf("'%s' is NULL.", name), call. = FALSE)
  }
  
  m <- as.matrix(x)
  if (!is.matrix(m)) {
    stop(sprintf("'%s' cannot be coerced to a matrix.", name), call. = FALSE)
  }
  had_warning <- FALSE
  withCallingHandlers({
    storage.mode(m) <- "numeric"
  }, warning = function(w) {
    had_warning <- TRUE
    invokeRestart("muffleWarning")
  })

  if (isTRUE(had_warning) || anyNA(m)) {
    warning(
      sprintf("'%s' contains non-numeric values that became NA.", name),
      call. = FALSE
    )
  }
  
  m
}

.check_same_rownames <- function(..., name = "inputs") {
  xs <- list(...)
  
  if (length(xs) < 2) {
    stop("At least two objects are required to compare rownames.", call. = FALSE)
  }
  
  rn <- lapply(xs, rownames)
  
  if (any(vapply(rn, is.null, logical(1)))) {
    stop(sprintf("All %s must have rownames.", name), call. = FALSE)
  }
  
  ref <- rn[[1]]
  ok <- all(vapply(rn[-1], identical, logical(1), ref))
  
  if (!ok) {
    stop(sprintf("Rownames do not match across %s.", name), call. = FALSE)
  }
  
  invisible(TRUE)
}

# ---- Backward-compatible wrappers (NOT exported) ----
# These are plain internal functions (no roxygen).
as_numeric_matrix <- function(x, name = "x") .as_numeric_matrix(x, name = name)
check_same_rownames <- function(..., name = "inputs") .check_same_rownames(..., name = name)
