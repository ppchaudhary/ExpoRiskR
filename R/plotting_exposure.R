# plotting_exposure.R
# ------------------------------------------------------------

#' Plot exposure perturbation ranking
#'
#' @param scores A data.frame from exposure_perturbation_score().
#' @param top_n Show only top N exposures (default 20). Use NULL for all.
#'
#' @return A ggplot object.
#'
#' @examples
#' d <- generate_dummy_exporisk(n = 30, p_micro = 10, p_metab = 12, p_expo = 4)
#' al <- align_omics(d$microbiome, d$metabolome, d$exposures, d$meta,
#'                  id_col = "sample_id", strict = TRUE)
#' pr <- prep_omics(al$microbiome, al$metabolome, al$exposures)
#' scores <- exposure_perturbation_score(pr$X, pr$Y, pr$E,
#'                                      fdr = 0.5, max_pairs = 120, seed = 1)
#' plot_exposure_ranking(scores)
#'
#' @export
plot_exposure_ranking <- function(scores, top_n = 20) {
  
  # ---- dependency guard ----
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_exposure_ranking().", call. = FALSE)
  }
  
  if (!is.data.frame(scores)) {
    stop("`scores` must be a data.frame.", call. = FALSE)
  }
  
  req <- c("exposure", "perturbation_score")
  if (!all(req %in% colnames(scores))) {
    stop(
      "`scores` must contain columns: ",
      paste(req, collapse = ", "),
      call. = FALSE
    )
  }
  
  df <- scores
  df <- df[order(df$perturbation_score, decreasing = TRUE), , drop = FALSE]
  
  if (!is.null(top_n)) {
    df <- utils::head(df, top_n)
  }
  
  df$exposure <- factor(df$exposure, levels = rev(df$exposure))
  
  ggplot2::ggplot(df, ggplot2::aes(x = exposure, y = perturbation_score)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Exposure",
      y = "Perturbation score (higher = more network change)",
      title = "Exposure perturbation ranking"
    )
}
