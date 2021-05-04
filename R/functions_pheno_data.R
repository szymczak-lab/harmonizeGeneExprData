
#' Extract phenotype data
#'
#' Extract phenotype data from GEO or ArrayExpress depending on study
#' identifier.
#'
#' All phenotype data available at GEO or ArrayExpress is returned without any
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
   pheno.all = extract_pheno_data_AE(
     ae.id = study.id,
     temp.dir = temp.dir)
  } else {
    stop(paste(study.id, "unknown! Use valid GEO or ArrayExpress study
               identifier starting with GSE or E-MTAB"))
  }

  ## check organism
  col.org = grep(
    "organism",
    colnames(pheno.all),
    ignore.case = TRUE)
  col.org = setdiff(
    col.org,
    grep(
      "organism.part",
      colnames(pheno.all),
      ignore.case = TRUE))

  if (length(col.org) == 0) {
    stop("information about organism missing!")
  } else if (length(col.org) > 1) {
    stop("more than one column for organism found!")
  }
  if (any(pheno.all[, col.org] != "Homo sapiens")) {
    stop("some samples have wrong organism!")
  }

  ## other values for NA
  na.strings = c("NA", "N/A",
                 "unknown",
                 "not reported")
  for (c in 1:ncol(pheno.all)) {
    pheno.all[grepl(
      paste(na.strings, collapse = "|"),
      pheno.all[, c]), c] = NA
  }

  return(pheno.all)
}


#' Harmonize phenotype data
#'
#' Extend original phenotype data by new variables with harmonized names and
#' values including sample and subject identifiers. Values of numeric variables
#' are converted to numeric. Values of character variables are mapped to
#' pre-specified values. For time variables numbers are extracted and
#' converted to days. If a study contains only patients of a specific
#' disease, this variable can be set globally.
#'
#' info.var contains information about all variables that should be harmonized
#' within a project, i.e. across studies. The list needs to be named by the
#' names of the harmonized variables and each element contains:
#' \itemize{
#' \item type: either "character", "numeric" or "time"
#' \item values: list named by final value and original values that should be
#' mapped (can include regular expressions)
#' }
#'
#' @param project [character(1)] name of project (used as prefix for all
#' harmonized variables)
#' @param pheno [data.frame] original phenotype data (e.g. as returned by
#' \code{\link{extract_pheno_data}})
#' @param info.var [list] project level information about variables that should
#' be harmonized (see Details)
#' @param col.id [character(1)] column with information about subject
#' identifiers
#' @param ind.use.id [numeric(1)] part of subject identifier that should be
#' kept after splitting using " ", "_" or "-"
#' @param cols.use [vector(n)] vector of columns used for harmonization, named
#' by variable names given in info.var
#' @param disease [character(1)] disease that should be set for all samples
#'
#' @return [data.frame] harmonized phenotype data
#' @export
#'
#' @examples
#' info.var = list(
#'   age = list(type = "numeric"),
#'   sex = list(
#'     type = "character",
#'     values = list(
#'       female = c("^f$", "female", "woman"),
#'       male = c("^m$", "^male", "^man"))
#'   ))
harmonize_pheno_data <- function(
  project,
  pheno,
  info.var,
  col.id,
  ind.use.id = NULL,
  cols.use,
  disease = NULL) {

  ## checks
  cols.not = setdiff(cols.use, colnames(pheno))
  if (length(cols.not) > 0) {
    stop(paste("some variables not available in pheno:",
               paste(cols.not, collapse = ", ")))
  }

  wrong.names = setdiff(names(cols.use), names(info.var))
  if (length(wrong.names) > 0) {
    stop(paste("some variables not available in info.var:",
               paste(wrong.names, collapse = ", ")))
  }

  ## harmonize each variable
  new = lapply(seq_len(length(cols.use)), function(i) {
    col = cols.use[i]
    print("===========================================")
    print(col)
    x = pheno[, col]
    info = info.var[[names(cols.use)[i]]]
    if (info$type == "numeric") {
      x.new = convert_numeric_variable(
        x = x)
    } else if (info$type == "character") {
      if (!is.null(info$values)) {
        x.new = map_variable(
          x = x,
          values = info$values)
      } else {
        x.new = tolower(gsub(" ", "_", x))
      }
    } else if (info$type == "time") {
      x.new = convert_time_variable(
        x = x,
        values = info$values)
    } else {
      stop("variable type must be 'character' or 'numeric'")
    }
  })
  new = do.call(data.frame, new)
  colnames(new) = names(cols.use)

  ## add information about age group
  if ("age" %in% colnames(new)) {
    new$age_group = ifelse(
      is.numeric(new$age),
      ifelse(new$age >= 18, "adult", "pediatric"),
      NA)
  }

  ## add information about disease
  if (!is.null(disease)) {
    new$disease = rep(disease, nrow(new))
  }

  ## set lesional to "healthy_control" for controls
  if ("lesional" %in% colnames(new) &
      "disease" %in% colnames(new)) {
    ind.ctrl = which(new$disease == "healthy_control")
    new$lesional[ind.ctrl] = "healthy_control"
  }

  ## set tissue to skin if only information about lesional is available
  if ("lesional" %in% colnames(new) &
      !("tissue" %in% colnames(new))) {
    tissue = rep(NA, nrow(new))
    tissue[grepl("lesional", new$lesional)] = "skin"
    new$tissue = tissue
  }

  ## add columns without information
  var.add = setdiff(names(info.var), colnames(new))
  info.add = matrix(
    NA,
    nrow = nrow(new),
    ncol = length(var.add))
  colnames(info.add) = var.add
  new = data.frame(
    new,
    info.add,
    stringsAsFactors = FALSE)

  ## identifier
  info.id = extract_identifier(pheno = pheno,
                               col.id = col.id,
                               ind.use.id = ind.use.id)
  new = data.frame(
    info.id,
    new,
    stringsAsFactors = FALSE)
  colnames(new) = paste(project, colnames(new), sep = "_")

  ## combine
  pheno.new = data.frame(new,
                         pheno)
  return(pheno.new)
}

#' Harmonize values of a character variable
#'
#' Map values of a character variable to pre-specified values. Print table of
#' mapping result.
#'
#' @keywords internal
map_variable <- function(
  x,
  values) {

  ## mapping to each category
  info.map = sapply(
    values,
    function(v) {
      base::grepl(paste(v, collapse = "|"), x, ignore.case = TRUE)
    })
  info.map[which(is.na(x)), ] = NA

  ## check that mapping is unique
  sum = apply(info.map, 1, sum)
  if (all(sum == 0)) {
    print("none of the values could be mapped:")
    print(sort(unique(x)))
    stop()
    #return(rep(NA, length(x)))
  }
  if (any(sum == 0, na.rm = TRUE)) {
    print("the following values were not mapped:")
    print(sort(unique(x[which(sum == 0)])))
    stop()
  }
  if (any(sum > 1, na.rm = TRUE)) {
    print("the following values were not uniquely mapped:")
    print(sort(unique(x[which(sum > 1)])))
    stop()
  }

  x.new = apply(info.map, 1, function(y) {
    if (all(is.na(y))) {
      return(NA)
    } else {
      names(values)[which(y)[1]]
    }
  })

  print(table(x, x.new, exclude = NULL))

  x.new = as.character(gsub(" ", "_", x.new))

  return(x.new)
}

#' Convert variable to numeric
#'
#' Extract numbers if necessary and convert to numeric.
#'
#' @keywords internal
convert_numeric_variable <- function(x) {

  if (is.numeric(x)) return(x)

  x.new = rep(NA, length(x))
  ind.not.na = which(!is.na(x))
  x.new[ind.not.na] = unlist(filesstrings::str_extract_numbers(
    string = x[ind.not.na],
    decimals = TRUE,
    negs = TRUE,
    leave_as_string = FALSE))

  ## check if only a single number was extracted for each element
  if (length(x.new) != length(x)) {
    stop("different length after conversion to numeric")
  }

  ## check if conversion worked for each element
  ind.na.new = which(is.na(x.new) & !is.na(x))
  if (length(ind.na.new) > 0) {
    print(sort(unique(x[ind.na.new])))
    stop("some values not converted to numeric")
  }

  return(x.new)
}


#' Convert variable to time in days
#'
#' Extract numbers if necessary, convert to numeric and transform to days.
#'
#' @keywords internal
convert_time_variable <- function(
  x,
  values) {

  ## replace baseline with 0
  if (!is.null(values)) {
    for (i in seq_len(length(values))) {
      x = gsub(
        paste(values[[i]], collapse = "|"),
        names(values)[i],
        x,
        ignore.case = TRUE)
    }
  }

  ## extract numbers
  x.num = convert_numeric_variable(x = x)

  ## identify unit
  units = list(
#    day = c("day", "d"),
    week = c("week", "wk", "w"),
    month = c("month", "m"))

  info.unit = sapply(
    units,
    function(u) {
      base::grepl(paste(u, collapse = "|"), x, ignore.case = TRUE)
    })
  info.unit[which(is.na(x)), ] = NA

  ## convert
  ind.week = which(info.unit[, "week"])
  x.num[ind.week] = x.num[ind.week] * 7
  ind.month = which(info.unit[, "month"])
  x.num[ind.month] = x.num[ind.month] * 30.5

  print(table(x, x.num))

  return(x.num)
}






#' Extract sample and subject identifiers
#'
#' @keywords internal
extract_identifier <- function(
  pheno,
  col.id,
  ind.use.id = NULL) {

  ## extract sample identifiers depending on repository
  if ("geo_accession" %in% colnames(pheno)) { ## GEO
    info = pheno[, c("geo_accession", "title")]
    colnames(info) = c("accession", "name")

  } else if ("Source.Name" %in% colnames(pheno)) { ## AE
    info = pheno[, c("Source.Name"), drop = FALSE]
    colnames(info) = "accession"
  }

  ## check if sample identifiers are unique
  if (any(duplicated(info$accession))) {
    warning("some accessions are duplicated!")
  } else {
    rownames(info) = info$accession
  }

  ## extract subject identifier
  if (!(col.id %in% colnames(pheno))) {
    stop(paste("identifier column", col.id, "not found"))
  }
  id = pheno[, col.id]
  if (!is.null(ind.use.id)) {
    id = extract_subject_id(id = id,
                            ind.use = ind.use.id)
  }
  info$id = id

  return(info)
}

#' Extract part of subject identifiers
#'
#' @keywords internal
extract_subject_id <- function(id,
                               ind.use) {

  temp = strsplit(id, split = " |_|-")
  id.new = sapply(temp, function(x) {
    if (length(x) == 1) {
      return(x)
    } else {
      return(paste(x[ind.use], collapse = "_"))
    }

  })
  return(id.new)
}


#' Check subject level information
#'
#' Check if subject information is consistent across all samples.
#'
#' @param pheno [data.frame] harmonized phenotype data
#' @param cols.subject [vector(n)] column names to be used for check
#'
#' @export
check_subject_info <- function(pheno,
                               cols.subject) {

  col.id = grep("_id$", colnames(pheno), value = TRUE)
  col.id = setdiff(col.id, "platform_id")
  if (length(col.id) == 0) {
    stop("no subject identifier column found")
  }
  if (length(col.id) > 1) {
    stop("more than one subject identifier column found")
  }

  cols.miss = setdiff(
    cols.subject,
    colnames(pheno))
  if (length(cols.miss) > 0) {
    stop(paste(
      "some subject specific columns missing:",
      paste(cols.miss, collapse = ",")))
  }

  ids = unique(pheno[, col.id])
  if (length(ids) < nrow(pheno)) {
    cols.use = unlist(sapply(cols.subject, function(x) {
      if (all(is.na(pheno[, x]))) {
        return(NULL)
      } else {
        return(x)
      }}))
    for (id in ids) {
      temp = pheno[pheno[, col.id] == id, cols.use, drop = FALSE]
      if (nrow(unique(temp)) > 1) {
        warning(paste("different characteristics for",
                      "subject", id, "\n"))
      }
    }
  }
}
