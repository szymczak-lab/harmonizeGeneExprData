#' Prepare count data
#'
#' Remove genes with only zero counts across samples and perform voom
#' transformation (as implemented in \code{\link[limma]{voom}})
#'
#' @param counts [data.frame or matrix] raw counts with genes in rows and
#' samples in columns
#'
#' @return [list] assays (counts, voom and voom.weights) to be used in
#' \code{\link{make_se_object}}
#' @export
#'
#' @examples
#' # define temporary directory for storing file from recount
#' temp.dir = tempdir()
#'
#' # get counts from recount
#' library(recount)
#' download_study(
#'   project = "SRP057087",
#'   type = "counts-gene",
#'   outdir = temp.dir)
#' counts = read.table(
#'   file.path(temp.dir, "counts_gene.tsv.gz"),
#'   header = TRUE,
#'   row.names = 29)
#'
#' assays = prepare_count_data(counts)

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
    expr.voom = v$E,
    voom.weights = weights)

  return(assays)
}

#' Convert SRR identifiers to GSM identifiers (GEO)
#'
#' Maps SRR identifiers to GSM identifiers using information from the SRA ftp
#' server (see Details). Note that this function assumes that the command line
#' tool grep has been installed to extract relevant information for the study.
#'
#' The global.mapping.file SRA_Accessions.tab has to be downloaded from
#' https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab
#' If study.mapping.file exists, it will be loaded without repeating the
#' extraction process.
#'
#' @param counts [data.frame or matrix] raw counts with genes in rows and
#' samples in columns
#' @param global.mapping.file [character(1)] path to local copy of
#' SRA_Accessions.tab
#' @param study.mapping.file [character(1)] name of file to save study specific
#' mapping information
#'
#' @return [data.frame or matrix] raw counts with colnames replaced by GSM
#' identifiers
#' @export
#'
#' @examples
#' \dontrun{
#' # define temporary directory for storing file from recount
#' temp.dir = tempdir()
#'
#' # get counts from recount
#' library(recount)
#' recount::download_study(
#'   project = "SRP057087",
#'   type = "counts-gene",
#'   outdir = temp.dir)
#' counts = read.table(
#'   file.path(temp.dir, "counts_gene.tsv.gz"),
#'   header = TRUE,
#'   row.names = 29)
#'
#' # mapping
#' counts.new = SRR_2_GSM(
#'   counts = counts,
#'   global.mapping.file = system.file(
#'     "extdata",
#'     "SRA_Accessions_example.tab",
#'     package = "harmonizeGeneExprData"),
#'   study.mapping.file = file.path(temp.dir, "mapping.rds"))
#'   }

SRR_2_GSM <- function(
  counts,
  global.mapping.file,
  study.mapping.file) {

  if (file.exists(study.mapping.file)) {
    print("using existing file ...")
    info = readRDS(study.mapping.file)
  } else {
    if (!file.exists(global.mapping.file)) {
      stop(paste("global.mapping.file", global.mapping.file, "does not exist"))
    }
    print("extracting mapping ...")
    info = NULL
    for (s in colnames(counts)) {
      cmd = paste("grep", s, global.mapping.file)
      match = try(system(cmd, intern = TRUE))
      ind = grep(paste0(s, "\t"), match)
      if (length(ind) != 1) {
        stop(paste("ID", s, "not found or multiple matches!"))
      }
      match = match[ind]
      info = rbind(info,
                   unlist(strsplit(match, "\t")))
    }
    names = unlist(strsplit(
      readLines(global.mapping.file,
                n = 1),
      "\t"))
    info = data.frame(
      info,
      stringsAsFactors = FALSE)
    dimnames(info) = list(
      info[, 1],
      names)
    info$Alias = gsub("_r[1:9]$", "", info$Alias)
    saveRDS(info,
            file = study.mapping.file)

  }

  if (nrow(info) < ncol(counts)) {
    stop("some SRR ids not found!")
  }

  colnames(counts) = info[colnames(counts), "Alias"]

  return(counts)
}


#' Prepare array data
#'
#' Download raw data from GEO or ArrayExpress and normalize it using the SCAN
#' and UPC approaches (as implemented in \code{\link[SCAN.UPC]{SCANfast}}
#' and \code{\link[SCAN.UPC]{UPCfast}}).
#'
#' @param study.id [character(1)] study identifier (GSEXXX for GEO, E-MTAB-XXX
#' for ArrayExpress)
#' @param pheno [data.frame] phenotype data (needs to contain column with
#' sample identifiers (name contains "accession"))
#' @param array.type [character(1)] type of array ("affy" (default), "agilent"
#' or "illumina")
#' @param platform [character(1)] array platform (needs to be specified for
#' Illumina arrays (?))
#' @param temp.dir [character(1)] directory where temporary files should be
#' stored
#'
#' @return [list] assays (expr and upc) to be used in
#' \code{\link{make_se_object}}
#' @export

prepare_array_data <- function(
  study.id,
  pheno,
  array.type = "affy",
  platform = NULL,
  temp.dir = tempdir()) {

  if (!(array.type %in% c("affy", "agilent", "illumina"))) {
    stop("array.type needs to be 'affy', 'agilent' or 'illumina'")
  }

  ## get SCAN normalized expression data
  expr = get_expr_data(
    study.id = study.id,
    pheno = pheno,
    type = "expr",
    array.type = array.type,
    temp.dir = temp.dir,
    platform = platform)

  ## get UPC values
  upc = get_expr_data(
    study.id = study.id,
    pheno = pheno,
    type = "upc",
    array.type = array.type,
    temp.dir = temp.dir,
    platform = platform)

  if (any(rownames(expr) != rownames(upc))) {
    stop("some rownames in expr and upc do not match")
  }

  assays = list(expr = expr, upc = upc)
  return(assays)
}


## platform: Illumina array (name used to load corresponding Bioconductor annotation package)
#' Download and normalize raw array data
#'
#' Use specific functions for each array type.
#'
#' @keywords internal
get_expr_data <- function(
  study.id,
  pheno,
  type = "expr",
  array.type = "affy",
  platform = "illuminaHumanv4",
  temp.dir = tempdir()) {

  if (type == "expr") {
    file.name = file.path(
      temp.dir,
      paste0(study.id, "_norm_scan.txt"))
  } else if (type == "upc") {
    file.name = file.path(
      temp.dir,
      paste0(study.id, "_upc.txt"))
  } else {
    stop(paste("type", type, "unknown"))
  }

  # if (array.type == "agilent") {
  #   if (!file.exists(file.name)) {
  #     if (type == "expr") {
  #       print("running SCAN normalization ...")
  #       eset = SCAN_TwoColor(gse,
  #                            outFilePath = file.name)
  #     } else if (type == "upc") {
  #       print("running UPC ...")
  #       eset = UPC_TwoColor(gse,
  #                           outFilePath = file.name)
  #     }
  #   }
  #   dat = as.matrix(utils::read.table(file = file.name))
  #
  #   ## keep only channel 1
  #   dat = dat[, grepl("Channel1$", colnames(dat))]
  #
  #   ## set colnames to GSM
  #   colnames(dat) = sapply(colnames(dat), function(x) {
  #     unlist(strsplit(x, "_"))[1]
  #   })
  #
  # } else if (array.type == "affy") {
    if (!file.exists(file.name)) {
      run_norm_scan_affy(study.id = study.id,
                         pheno = pheno,
                         type = type,
                         temp.dir = temp.dir,
                         file.name = file.name)
    }
    dat = as.matrix(utils::read.table(file = file.name))

  # } else if (array.type == "illumina") {
  #   if (!file.exists(file.name)) {
  #     run_norm_scan_illumina(gse = gse,
  #                            expr.dir = expr.dir,
  #                            type = type,
  #                            res.file = file.name,
  #                            platform = platform,
  #                            numCores = 1)
  #   }
  #   dat = utils::read.table(file = file.name,
  #                    header = TRUE,
  #                    as.is = TRUE,
  #                    check.names = FALSE,
  #                    sep = "\t")
  #   dat = as.matrix(dat)
  #
  # }

  # if (array.type != "illumina") {
  #   if (!all(rownames(pheno) %in% colnames(dat))) {
  #     stop(paste("data for some subjects is missing in", gse))
  #   }
  #
  #   dat = dat[, rownames(pheno)]
  #
  # }
  return(dat)
}


# ## calling SCANfast with GSE identifier does not work for some studies (e.g. GSE102725)
# ## since in addition to CEL files additional files are stored in tar file with raw dat
#' Download and normalize CEL files (Affy)
#'
#' For each sample download and normalize CEL file using BrainArray annotation
#' for which relevant package will be installed. Keep only genes with Ensembl
#' gene identifier and remove suffix '_at'.
#'
#' @keywords internal
run_norm_scan_affy <- function(
  study.id,
  pheno,
  type,
  temp.dir = tempdir(),
  file.name) {

  ## get sample ids
  col.ids = "geo_accession"
  if (!(col.ids %in% colnames(pheno))) {
    col.ids = grep("_accession", colnames(pheno), value = TRUE)
    if (length(col.ids) == 0) {
      stop(paste("pheno has no column with sample ids",
                 "(name containing 'accession'"))
    }
    if (length(col.ids) > 1) {
      stop(paste("pheno has more than one column with sample ids",
           "(name containing 'accession'"))
    }
  }
  sample.ids = pheno[, col.ids]

  res.all = NULL
  pkg.name = NULL
  for (s in sample.ids) {
    sample.file = file.path(temp.dir,
                            paste0(s, "_", type, ".txt"))

    ## download CEL file
    if (grepl("GSE", study.id)) {
      cel.file = download_cel_file_GEO(
        sample = s,
        temp.dir = temp.dir)
    }

    if (is.null(pkg.name) & s == sample.ids[1]) {
      pkg.name = SCAN.UPC::InstallBrainArrayPackage(
        celFilePath = cel.file,
        version = "25.0.0",
        organism = "hs",
        annotationSource = "ensg")
    }

    if (!file.exists(sample.file)) {
      if (type == "expr") {
        temp = SCAN.UPC::SCANfast(
          celFilePattern = cel.file,
          probeSummaryPackage = pkg.name,
          outFilePath = sample.file)
      } else if (type == "upc") {
        temp = SCAN.UPC::UPCfast(
          celFilePattern = cel.file,
          probeSummaryPackage = pkg.name,
          outFilePath = sample.file)
      }
      rm(temp)
      gc = gc()
    }
    res = utils::read.table(sample.file)

    if (is.null(res.all)) {
      res.all = res
    } else {
      if (any(rownames(res.all) != rownames(res))) {
        stop(paste("different rownames for sample", sample))
      }
      res.all = data.frame(res.all, res)
    }
  }

  ## rename colnames
  colnames(res.all) = sample.ids

  ## extract and rename Ensembl gene ids
  res.all = res.all[grepl("^ENSG", rownames(res.all)), ]
  rownames(res.all) = gsub("_at", "", rownames(res.all))

  utils::write.table(
    res.all,
    sep = "\t",
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE,
    file = file.name)
}

