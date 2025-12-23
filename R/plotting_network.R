#' Plot exposure-adjusted multi-omics network (bipartite)
#'
#' @description
#' Plots a bipartite igraph network returned by \code{build_exposure_network()}.
#' Uses base igraph plotting (no extra dependencies).
#'
#' @param net A list returned by \code{build_exposure_network()} with elements
#'   \code{$graph} and \code{$edges}.
#' @param file Optional output filename. If provided, saves a PNG (recommended).
#' @param width,height Plot device size (in inches) when saving.
#' @param dpi DPI when saving PNG.
#' @param layout Layout function name passed to igraph. Default
#'   \code{"layout_with_fr"}.
#' @param max_label_nodes Max nodes to label (largest by degree). Default 30.
#' @examples
#' d <- generate_dummy_exporisk(seed = 1, n = 12, p_micro = 5, p_metab = 6, p_expo = 3)
#' al <- align_omics(d$microbiome, d$metabolome, d$exposures, d$meta,
#'                  id_col = "sample_id", strict = TRUE)
#' pr <- prep_omics(al$microbiome, al$metabolome, al$exposures)
#' net <- build_exposure_network(pr$X, pr$Y, pr$E, fdr = 0.95, max_pairs = 120, seed = 1)
#' plot_exposure_network(net)

#'
#' @return Invisibly returns \code{net$graph}.
#'
#' @export
plot_exposure_network <- function(net,
                                  file = NULL,
                                  width = 10, height = 7, dpi = 300,
                                  layout = "layout_with_fr",
                                  max_label_nodes = 30) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for plot_exposure_network().", call. = FALSE)
  }
  
  if (is.null(net) || !is.list(net) || is.null(net$graph)) {
    stop("`net` must be a list returned by build_exposure_network() containing `$graph`.",
         call. = FALSE)
  }
  
  g <- net$graph
  
  if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
    warning("Network has no edges; nothing to plot.", call. = FALSE)
    return(invisible(g))
  }
  
  # ---- weights ----
  w <- igraph::E(g)$weight
  if (is.null(w)) w <- rep(1, igraph::ecount(g))
  
  # edge width by abs(weight)
  ew <- abs(w)
  if (all(!is.finite(ew)) || max(ew, na.rm = TRUE) == 0) {
    ew <- rep(1, length(ew))
  } else {
    ew <- ew / max(ew, na.rm = TRUE)
    ew[!is.finite(ew)] <- 0
    ew <- 1 + 4 * ew
  }
  
  # edge color by sign
  ec <- ifelse(w >= 0, "#2C7BB6", "#D7191C")
  
  # node size by degree
  deg <- igraph::degree(g)
  if (all(deg == 0)) {
    vsz <- rep(6, length(deg))
  } else {
    vsz <- 4 + 6 * (deg / max(deg))
  }
  
  # node color by bipartite type: TRUE=metabolite per builder
  vtype <- igraph::V(g)$type
  if (is.null(vtype)) vtype <- rep(FALSE, igraph::vcount(g))
  vcol <- ifelse(vtype, "#4DAF4A", "#984EA3")
  
  # label only top-degree nodes
  ord <- order(deg, decreasing = TRUE)
  lab_idx <- ord[seq_len(min(max_label_nodes, length(ord)))]
  vlab <- rep("", igraph::vcount(g))
  vlab[lab_idx] <- igraph::V(g)$name[lab_idx]
  
  # ---- Layout (FIX: FR layout requires positive weights) ----
  lay_fun <- get(layout, envir = asNamespace("igraph"), inherits = TRUE)
  
  # Use positive weights for layout ONLY (do not change scientific weights)
  w_layout <- abs(w)
  w_layout[!is.finite(w_layout)] <- 1
  w_layout[w_layout <= 0] <- 1
  
  coords <- tryCatch(
    {
      # layout_with_fr supports weights; passing positive weights avoids error
      lay_fun(g, weights = w_layout)
    },
    error = function(e) {
      # fall back to unweighted layout if the chosen layout doesn't accept weights
      lay_fun(g)
    }
  )
  
  # ---- Save if requested ----
  if (!is.null(file)) {
    grDevices::png(filename = file, width = width, height = height, units = "in", res = dpi)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  
  graphics::par(mar = c(0, 0, 2, 0))
  igraph::plot.igraph(
    g,
    layout = coords,
    vertex.size = vsz,
    vertex.color = vcol,
    vertex.label = vlab,
    vertex.label.cex = 0.7,
    vertex.label.color = "black",
    edge.width = ew,
    edge.color = ec,
    main = "Exposure-adjusted Microbe-Metabolite Network"
  )
  
  invisible(g)
}
