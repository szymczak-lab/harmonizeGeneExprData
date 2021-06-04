#' Extract meta data
#'
#' Extract meta data from GEO or ArrayExpress depending on study
#' identifier.
#'
#' Different meta data information is available at GEO and ArrayExpress.
#'
#' @param study.id [character(1)] study identifier (GSEXXX for GEO, E-MTAB-XXX
#' for ArrayExpress)
#' @param temp.dir [character(1)] directory where temporary files should be
#' stored
#' @param platform [character(1)] platform to be selected in GEO
#' studies which use several platforms
#'
#' @return [list] meta data (content depending on information available at GEO
#' or ArrayExpress)
#' @export
#'
#' @examples
#' meta.data.l = extract_meta_data(study.id = "GSE67785")

extract_meta_data <- function(
  study.id,
  temp.dir = tempdir(),
  platform = NULL) {

  if (grepl("^GSE", study.id)) {
    meta.data.l = extract_meta_data_GEO(
      gse = study.id,
      temp.dir = temp.dir,
      platform = platform)

  }  else if (grepl("^E-MTAB", study.id)) {
    meta.data.l = extract_meta_data_AE(
      ae.id = study.id,
      temp.dir = temp.dir)
  }

  return(meta.data.l)
}
