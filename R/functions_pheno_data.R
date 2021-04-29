
#' Extract phenotype data
#'
#' Extract phenotype data from GEO or ArrayExpress depending on study
#' identifier.
#'
#' All phenotype data available at GEO or AE is returned without any
#' modifications.
#'
#' @param study.id [character(1)] study identifier (GSEXXX for GEO, E-MTAB-XXX
#' for ArrayExpress)
#' @param temp.dir [character(1)] directory where temporary files should be
#' stored
#' @param platform [character(1)] platform to be selected in GEO
#' studies which use several platforms
#'
#' @return [data.frame] phenotype data for each sample and all variables
#' provided by GEO or ArrayExpress
#' @export
extract_pheno_data <- function(
  study.id,
  temp.dir = tempdir(),
  platform = NULL) {

  if (grepl("^GSE", study.id)) {
    pheno.all = extract_pheno_data_GEO(
      gse = study.id,
      temp.dir = temp.dir,
      platform = platform)
  } else if (grepl("^E-MTAB", study.id)) {
#    pheno.all = extract_pheno_data_AE(
#      ae.id = study.id,
#      temp.dir = temp.dir)
  } else {
    stop(paste(study.id, "unknown! Use valid GEO or ArrayExpress study
               identifier starting with GSE or E-MTAB"))
  }

  return(pheno.all)
}
