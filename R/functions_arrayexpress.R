
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

  idf = scan(
    file.path(temp.dir, info.ae$idf),
    character(),
    sep = "\n",
    encoding = "UTF-8")

  meta.data.l = list(
    study_id = ae.id,
    contact_institute = extract_info_from_idf(
      idf = idf,
      label = "Person Affiliation"),
    contact_email = extract_info_from_idf(
      idf = idf,
      label = "Person Email"),
    date = extract_info_from_idf(
      idf = idf,
      label = "Public Release Date"),
    platform_name = extract_info_from_idf(
      idf = idf,
      label = "Protocol Hardware"),
    pubmed_id = extract_info_from_idf(
      idf = idf,
      label = "PubMed ID"))

  return(meta.data.l)
}


#' Extract information from IDF
#'
#' Extract specific information from IDF file
#'
#' @keywords internal
extract_info_from_idf <- function(idf, label) {
  l = unlist(strsplit(
    grep(label, idf, value = TRUE),
    "\t"))
  l = l[l != ""]
  final = ifelse(
    is.null(l), NA,
    ifelse(length(l) == 1, NA,
           paste(l[-1], collapse = "|")))
  return(final)
}

#' Download CEL file (Affy) from ArrayExpress
#'
#' @keywords internal
download_cel_file_AE <- function(
  sample,
  temp.dir = tempdir()) {

  cel.file = file.path(temp.dir, paste0(sample, ".CEL"))

  if (!file.exists(cel.file)) {
    ## download raw data (all and unzip)
    info.ae = ArrayExpress::getAE(
      accession = sample,
      path = temp.dir,
      type = "raw")

    if (!file.exists(cel.file)) {
      stop(paste("no CEL file found for sample", sample))
    }
  }

  return(cel.file)
}

