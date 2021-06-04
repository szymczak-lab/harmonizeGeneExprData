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
