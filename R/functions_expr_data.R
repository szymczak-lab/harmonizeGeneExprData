#' Prepare count data
#'
#' Remove genes with only zero counts across samples and perform voom
#' transformation (as implemented in \code{\link[limma]{voom}})
#'
#' @param counts [data.frame or matrix] raw counts with genes in rows and
#' samples in columns
#'
#' @return [list] assays (counts, voom and voom.weights) to be used in XXX
#' @export
prepare_count_data <- function(counts) {

  counts = as.matrix(counts)

  ## remove genes with only zero counts
  total = rowSums(counts)
  counts = counts[total > 0, ]

  ## voom transformation
  v = limma::voom(counts)
  weights = v$weights
  dimnames(weights) = dimnames(v$E)
  assays = list(
    counts = counts,
    voom = v$E,
    voom.weights = weights)

  return(assays)
}
