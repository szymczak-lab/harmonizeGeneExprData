
#' Download files from ArrayExpress
#'
#' Download files from ArrayExpress and save in temporary directory. If files
#' are already available, they will not be downloaded again.
#'
#' @keywords internal
download_files_AE <- function(
  ae.id,
  temp.dir = tempdir()) {
  files = dir(
    temp.dir,
    pattern = ae.id)
  if (length(files) == 0) {
    info.ae = ArrayExpress::getAE(
      accession = ae.id,
      path = temp.dir,
      type = "processed")
  } else {
    info.ae = list(
      sdrf = grep("sdrf", files, value = TRUE),
      idf = grep("idf", files, value = TRUE))
  }

  return(info.ae)
}

#' Extract phenotype data from ArrayExpress
#'
#' Download the SDRF file from ArrayExpress and return it as a data.frame
#'
#' @keywords internal
extract_pheno_data_AE <- function(
  ae.id,
  temp.dir = tempdir()) {

  info.ae = download_files_AE(
    ae.id = ae.id,
    temp.dir = temp.dir)

  sdrf.file = file.path(
    temp.dir,
    info.ae$sdrf)
  sdrf = rio::import(
    sdrf.file,
    check.names = TRUE)

  return(sdrf)
}

## platform can be extracted from xml file at
#https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments/E-MTAB-6556/protocols
#' Extract meta data from ArrayExpress
#'
#' Download the IDF file from ArrayExpress and extract some of the information
#'
#' @keywords internal
extract_meta_data_AE <- function(
  ae.id,
  temp.dir = temp.dir) {

  info.ae = download_files_AE(
    ae.id = ae.id,
    temp.dir = temp.dir)

  idf = ArrayExpress:::readExperimentData(
    idf = info.ae$idf,
    path = temp.dir)

  if (idf@pubMedIds != "") {
    pubmed = idf@pubMedIds
  } else {
    pubmed = NULL
  }
  meta.data.l = list(ae.id = ae.id,
                     contact_institute = idf@lab,
                     contact_email = idf@contact,
                     pubmed_id = NULL)

  return(meta.data.l)
}


