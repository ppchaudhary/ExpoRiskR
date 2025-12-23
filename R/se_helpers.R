#' @importFrom SummarizedExperiment assay colData
#' @importFrom S4Vectors DataFrame
NULL

.is_se <- function(x) inherits(x, "SummarizedExperiment")

.get_se_assay <- function(se, assay_name = NULL) {
  if (!.is_se(se)) stop("Expected a SummarizedExperiment.", call. = FALSE)
  if (is.null(assay_name)) {
    assay_name <- SummarizedExperiment::assayNames(se)[1]
  }
  SummarizedExperiment::assay(se, assay_name)
}

.get_se_coldata_df <- function(se) {
  as.data.frame(SummarizedExperiment::colData(se))
}
