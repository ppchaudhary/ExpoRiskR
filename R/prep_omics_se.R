#' Preprocess SummarizedExperiment-based omics blocks and exposures
#'
#' @param aligned Output from align_omics_se() or align_omics().
#' @param assay_micro Assay name for microbiome SE (default: first assay).
#' @param assay_metab Assay name for metabolome SE (default: first assay).
#' @param ... Passed to prep_omics().
#'
#' @return A list with preprocessed matrices: X, Y, E.
#'
#' @examples
#' set.seed(8)
#' d <- generate_dummy_exporisk(n = 12, p_micro = 5, p_metab = 6, p_expo = 3)
#' aligned <- align_omics_se(d$microbiome, d$metabolome, d$exposures, d$meta,
#'                           id_col = "sample_id", strict = TRUE)
#' se2 <- prep_omics_se(aligned)
#' se2
#'
#' @export
prep_omics_se <- function(aligned,
                          assay_micro = NULL,
                          assay_metab = NULL,
                          ...) {
  
  if (!is.list(aligned)) {
    stop("`aligned` must be a list returned by align_omics_se() or align_omics().", call. = FALSE)
  }
  
  # Case 1: output of align_omics_se()
  if (all(c("se_microbiome", "se_metabolome", "exposures") %in% names(aligned))) {
    
    se_micro <- aligned$se_microbiome
    se_metab <- aligned$se_metabolome
    E <- aligned$exposures
    
    if (!inherits(se_micro, "SummarizedExperiment") ||
        !inherits(se_metab, "SummarizedExperiment")) {
      stop("`se_microbiome` and `se_metabolome` must be SummarizedExperiment objects.", call. = FALSE)
    }
    
    # pick assays (default: first)
    if (is.null(assay_micro)) assay_micro <- SummarizedExperiment::assayNames(se_micro)[1]
    if (is.null(assay_metab)) assay_metab <- SummarizedExperiment::assayNames(se_metab)[1]
    
    # SummarizedExperiment assays are features x samples; prep_omics expects samples x features
    X <- t(SummarizedExperiment::assay(se_micro, assay_micro))
    Y <- t(SummarizedExperiment::assay(se_metab, assay_metab))
    
    return(prep_omics(X, Y, E, ...))
  }
  
  # Case 2: output of align_omics()
  if (all(c("microbiome", "metabolome", "exposures") %in% names(aligned))) {
    return(prep_omics(aligned$microbiome, aligned$metabolome, aligned$exposures, ...))
  }
  
  stop("`aligned` must be the output of align_omics_se() or align_omics().", call. = FALSE)
}
