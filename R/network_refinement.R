#' Refine network using ARACNE DPI (wrapper)
#'
#' Convenience wrapper around \code{\link{prune_network_dpi}}.
#' Maintained for backward compatibility with the original function name.
#'
#' @param adj_matrix Matrix. Gene regulatory weight matrix to refine.
#' @param eps Numeric. DPI tolerance threshold (0-1), default 0.1.
#'
#' @return Refined matrix after DPI filtering.
#' @export
refine_network_aracne <- function(adj_matrix, eps = 0.1) {
  prune_network_dpi(adj_matrix, eps = eps)
}
