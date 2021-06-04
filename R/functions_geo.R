
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


#' Extract meta data from GEO
#'
#' Extract platform information and information from the experimentData()
#' slot from the ExpressionSet of a GEO study
#'
#' @keywords internal
extract_meta_data_GEO <- function(
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
  if (is.null(platform)) platform = unique(eset$platform_id)
  gpl = GEOquery::getGEO(platform)

  temp = Biobase::experimentData(eset)@other

  ## convert dates
  prev = Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
  temp$submission_date = as.Date(
    x = temp$submission_date,
    format = c("%b %d %Y"))
  temp$last_update_date = as.Date(
    x = temp$last_update_date,
    format = c("%b %d %Y"))
  Sys.setlocale("LC_TIME", prev)

  pubmed.id = temp$pubmed_id
  if (is.null(pubmed.id)) {
    pubmed.id = NA
  } else {
    pubmed.id = gsub("\\n", "|", temp$pubmed_id)
  }

  meta.data.l = list(
    study_id = gse,
    contact_institute = temp$contact_institute,
    contact_email = temp$contact_email,
    date = temp$submission_date,
    last_update_date = temp$last_update_date,
    type = temp$type,
    platform_id = temp$platform_id,
    platform_name = gpl@header$title,
    pubmed_id = pubmed.id,
    sra_id = stringr::str_extract(string = temp$relation,
                         pattern = "SRP[0-9]*"))

  meta.data.l$technology = ifelse(
    grepl("sequencing", meta.data.l$type),
    "rnaseq",
    ifelse(
      grepl("agilent",
            meta.data.l$platform_name,
            ignore.case = TRUE),
      "array.agilent",
      ifelse(
        grepl("illumina",
              meta.data.l$platform_name,
              ignore.case = TRUE),
        "array.illumina", "array.affy")))

  return(meta.data.l)
}

#' Download CEL file (Affy) from GEO
#'
#' @keywords internal
download_cel_file_GEO <- function(
  sample,
  temp.dir = tempdir()) {

  info.cel.file = NULL
  #  count = 0
  #  while (is.null(info.cel.file)) {
  #    if (count > 0 && count %% 10 == 0) {
  #      print(count/10)
  #    }
  info.cel.file = GEOquery::getGEOSuppFiles(
    GEO = sample,
    makeDirectory = FALSE,
    baseDir = temp.dir,
    fetch_files = FALSE,
    filter_regex = "CEL|cel")
  #    count = count + 1
  #  }
  if (nrow(info.cel.file) > 1) {
    stop(paste("more than one CEL file found for sample", sample))
  }
  if (nrow(info.cel.file) == 0) {
    stop(paste("no CEL file found for sample", sample))
  }

  cel.file = file.path(temp.dir, info.cel.file$fname)
  if (!file.exists(cel.file)) {
    utils::download.file(
      url = info.cel.file$url,
      destfile = cel.file)
  }
  return(cel.file)
}

