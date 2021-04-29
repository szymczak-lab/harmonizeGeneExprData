
#' Wrapper of the GEOquery::getGEO() function
#'
#' The wrapper function adds some basic validity checks including format of the
#' study identifier and number of different platforms used in the study.
#'
#' This functions results in an error if the study identifier does not start
#' with GSE and if more than one platform is used in a study. If the latter is
#' the case, the desired platform can be set using the platform parameter.
#'
#' @keywords internal
wrapper_getGEO <- function(
  gse,
  GSEMatrix = TRUE,
  AnnotGPL = FALSE,
  getGPL = FALSE,
  platform = NULL,
  temp.dir = tempdir()) {

  ## check that accession is a Series record
  if (!base::grepl("^GSE", gse)) {
    stop("accession needs to be a Series record (GSExxx)!")
  }

  ## download from GEO
  temp.l = GEOquery::getGEO(
    GEO = gse,
    GSEMatrix = GSEMatrix,
    AnnotGPL = AnnotGPL,
    getGPL = getGPL,
    destdir = temp.dir)

  ## some studies contain samples from different platforms which are stored as
  ## separate ExpressionSet objects
  if (length(temp.l) > 1) {
    if (is.null(platform)) {
      stop(paste("more than one ExpressionSet returned for accession",
                 gse, "and no platform specified!"))
    } else {
      info.platform = vapply(
        X = temp.l,
        FUN = BiocGenerics::annotation,
        FUN.VALUE = character(1))
      ind = which(info.platform == platform)
      if (length(ind) == 0) {
        stop(paste("platform", platform, "does not match information",
                   "about downloaded platforms:",
                   paste(info.platform, collapse = ", ")))
      } else {
        eset = temp.l[[ind]]
      }
    }

  } else {
    eset = temp.l[[1]]
  }
  return(eset)
}

#' Extract phenotype data from GEO
#'
#' Extract the pData() slot from the ExpressionSet of a GEO study
#'
#' @keywords internal
extract_pheno_data_GEO <- function(
  gse,
  temp.dir = tempdir(),
  platform = NULL) {

  eset = wrapper_getGEO(
    gse = gse,
    GSEMatrix = TRUE,
    AnnotGPL = FALSE,
    getGPL = FALSE,
    platform = platform,
    temp.dir = temp.dir)
  return(Biobase::pData(eset))
}
