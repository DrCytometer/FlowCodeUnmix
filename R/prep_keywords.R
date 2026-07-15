# prep_keywords.r

#' @title Prepare Keywords
#'
#' @description
#' Updates FCS file keywords after unmixing to define the new parameters. Tracks
#' existing keywords from the input FCS file for metadata compatibility.
#'
#' @param fcs.keywords The input keywords obtained from the read FCS file.
#' @param final.matrix The expression data, containing both unmixed data and any
#' retained parameters such as scatter and time.
#' @param original.param The original parameter (column) names of the input FCS
#' expression data.
#' @param spectra A matrix containing the spectral data.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns.
#' @param flow.control A list containing flow cytometry control parameters.
#' @param asp The AutoSpectral parameter list.
#' @param method A character string specifying the unmixing method used.
#' @param file.name The name of the FCS file to be written.
#' @param combos A character vector of the FlowCode combination tag names.
#' @param include.imaging A logical value indicating whether to include imaging
#' parameters in the written FCS file. Default is `FALSE`.
#'
#' @return The updated keyword list for writing the FCS file.
#'
#' @export


prep.keywords <- function(
    fcs.keywords,
    final.matrix,
    original.param,
    spectra,
    af.spectra,
    flow.control,
    asp,
    method,
    file.name,
    combos,
    include.imaging = FALSE
) {

  # Identify non-parameter keywords
  non.param.keys <- fcs.keywords[!grepl("^\\$?P\\d+", names(fcs.keywords))]
  if (asp$cytometer == "Mosaic") {
    non.param.keys <- non.param.keys[!grepl("^\\$?CH\\d+", names(non.param.keys))]
  }

  # Build parameter lookup from original file
  pN.keys <- grep("^\\$?P\\d+N$", names(fcs.keywords), value = TRUE)
  param.lookup <- lapply(pN.keys, function(k) {
    p.idx <- sub("\\$?P(\\d+)N", "\\1", k)
    matches <- grep(paste0("^\\$?P", p.idx, "(?:[A-Z]+)$"), names(fcs.keywords), value = TRUE)
    stats::setNames(fcs.keywords[matches], matches)
  })
  names(param.lookup) <- sapply(pN.keys, function(k) fcs.keywords[[k]])

  # Whitelist: only channels genuinely carried over from the original file may
  # match param.lookup. BD instruments (A5SE, Discover S8/A8) write
  # pre-unmixed FCS files whose $PnN values are fluorophore names (e.g.
  # "FITC-A", "PE-A") identical to FlowCodeUnmix's unmixed output column
  # names. Without this filter, param.lookup would match those names and
  # inject metadata from the wrong original parameter index into the new
  # file. Mirrors AutoSpectral::define.keywords().
  if (!is.null(asp$time.and.scatter)) {
    whitelist.channels <- asp$time.and.scatter
    if (isTRUE(include.imaging) && !is.null(asp$non.spectral.channel)) {
      ns.pattern <- paste0(asp$non.spectral.channel, collapse = "|")
      ns.matches <- grep(ns.pattern, names(param.lookup), value = TRUE)
      whitelist.channels <- union(whitelist.channels, ns.matches)
    }
  } else if (!is.null(asp$non.spectral.channel)) {
    ns.pattern <- paste0(asp$non.spectral.channel, collapse = "|")
    whitelist.channels <- grep(ns.pattern, names(param.lookup), value = TRUE)
  } else {
    whitelist.channels <- character(0)
  }
  param.lookup <- param.lookup[names(param.lookup) %in% whitelist.channels]

  af.n <- nrow(af.spectra)

  # Build new parameter keywords
  param.keywords <- list()
  n.param <- ncol(final.matrix)
  p.names <- colnames(final.matrix)

  for (i in seq_len(n.param)) {
    p.name <- p.names[i]
    # FCS 3.1 3.2.23: commas are not allowed in $PnN values as they conflict
    # with keywords such as $SPILLOVER and $TR. Replace any with underscores.
    if (grepl(",", p.name, fixed = TRUE)) {
      warning(sprintf(
        "Parameter name '%s' contains a comma, which is invalid in $PnN (FCS 3.1 3.2.23). Replacing with '_'.",
        p.name
      ), call. = FALSE)
      p.name <- gsub(",", "_", p.name, fixed = TRUE)
    }
    p_prefix <- paste0("$P", i)

    # 1. AF Index Logic
    if (p.name == "AF Index") {
      param.keywords[[paste0(p_prefix, "N")]] <- "AF Index"
      param.keywords[[paste0(p_prefix, "S")]] <- "Autofluorescence Index"
      param.keywords[[paste0(p_prefix, "B")]] <- "32"
      param.keywords[[paste0(p_prefix, "E")]] <- "0,0"
      param.keywords[[paste0(p_prefix, "R")]] <- as.character(af.n)
      param.keywords[[paste0(p_prefix, "G")]] <- "1"
      param.keywords[[paste0(p_prefix, "DISPLAY")]] <- "LIN"

      # 2. Main FlowCode Column
    } else if (p.name == "FlowCode") {
      param.keywords[[paste0(p_prefix, "N")]] <- p.name
      param.keywords[[paste0(p_prefix, "S")]] <- "FlowCode Intensity"
      param.keywords[[paste0(p_prefix, "B")]] <- "32"
      param.keywords[[paste0(p_prefix, "E")]] <- "0,0"
      param.keywords[[paste0(p_prefix, "R")]] <- as.character(asp$expr.data.max)
      param.keywords[[paste0(p_prefix, "G")]] <- "1"
      param.keywords[[paste0(p_prefix, "DISPLAY")]] <- "LOG"

      # 3. Individual Combo Tags
    } else if (p.name %in% combos) {
      param.keywords[[paste0(p_prefix, "N")]] <- p.name
      param.keywords[[paste0(p_prefix, "S")]] <- "Tag"
      param.keywords[[paste0(p_prefix, "B")]] <- "32"
      param.keywords[[paste0(p_prefix, "E")]] <- "0,0"
      param.keywords[[paste0(p_prefix, "R")]] <- as.character(asp$expr.data.max)
      param.keywords[[paste0(p_prefix, "G")]] <- "1"
      param.keywords[[paste0(p_prefix, "DISPLAY")]] <- "LOG"

      # 4. Autofluorescence
    } else if (p.name == "AF-A") {
      param.keywords[[paste0(p_prefix, "N")]] <- p.name
      param.keywords[[paste0(p_prefix, "S")]] <- "Autofluorescence"
      param.keywords[[paste0(p_prefix, "B")]] <- "32"
      param.keywords[[paste0(p_prefix, "E")]] <- "0,0"
      param.keywords[[paste0(p_prefix, "R")]] <- as.character(asp$expr.data.max)
      param.keywords[[paste0(p_prefix, "G")]] <- "1"
      param.keywords[[paste0(p_prefix, "DISPLAY")]] <- "LOG"

      # 5. Whitelisted carry-through parameters (Scatter, Time, and optionally imaging)
    } else if (p.name %in% names(param.lookup)) {
      old.entry <- param.lookup[[p.name]]
      keep.fields <- c("N", "S", "B", "E", "R", "G", "V", "DISPLAY", "TYPE")
      old.names <- names(old.entry)
      field.types <- sub("^\\$?P\\d+", "", old.names)
      keep.idx <- field.types %in% keep.fields

      filtered.entry <- old.entry[keep.idx]
      new_names <- paste0(p_prefix, field.types[keep.idx])
      for (k in seq_along(filtered.entry)) {
        param.keywords[[new_names[k]]] <- filtered.entry[[k]]
      }

      # 6. New Unmixed Fluorophores
    } else {
      bit.depth <- if (!is.null(asp$bit.depth)) asp$bit.depth else "32"
      param.keywords[[paste0(p_prefix, "N")]] <- p.name
      param.keywords[[paste0(p_prefix, "B")]] <- as.character(bit.depth)
      param.keywords[[paste0(p_prefix, "E")]] <- "0,0"
      param.keywords[[paste0(p_prefix, "R")]] <- as.character(asp$expr.data.max)
      param.keywords[[paste0(p_prefix, "G")]] <- "1"
      param.keywords[[paste0(p_prefix, "DISPLAY")]] <- "LOG"

      # Map Marker/Stain from flow.control
      clean.name <- sub("-A$", "", p.name)
      f.idx <- match(clean.name, flow.control$fluorophore)
      marker <- if (!is.na(f.idx)) as.character(flow.control$antigen[f.idx]) else p.name
      param.keywords[[paste0(p_prefix, "S")]] <- marker
    }
  }

  # FCS 3.1 3.2.18: $PnE and $PnR are required for every parameter.
  # 3.2.22: when $DATATYPE/F/, all parameters must have $PnE/0,0/.
  # Carry-through parameters copied from source files may be missing either;
  # apply safe fallbacks here.
  for (i in seq_len(n.param)) {
    p_prefix <- paste0("$P", i)
    e.key <- paste0(p_prefix, "E")
    r.key <- paste0(p_prefix, "R")
    if (is.null(param.keywords[[e.key]])) {
      param.keywords[[e.key]] <- "0,0"
    }
    if (is.null(param.keywords[[r.key]])) {
      param.keywords[[r.key]] <- as.character(asp$expr.data.max)
    }
  }

  # Format Spectra for Keywords
  format.matrix.string <- function(m) {
    vals <- as.vector(t(m))
    formatted <- formatC(vals, digits = 8, format = "fg", flag = "#")
    paste(c(nrow(m), ncol(m), rownames(m), colnames(m), formatted), collapse = ",")
  }

  # Combine everything
  new.keywords <- utils::modifyList(non.param.keys, param.keywords)

  # Package versions (custom keywords: no "$" prefix, which is reserved for
  # FCS-standard keys)
  asp.ver <- as.character(utils::packageVersion("AutoSpectral"))
  asp.rcpp.ver <- if (requireNamespace("AutoSpectralRcpp", quietly = TRUE)) {
    as.character(utils::packageVersion("AutoSpectralRcpp"))
  } else {
    "0"
  }
  fc.rcpp.ver <- if (requireNamespace("FlowCodeUnmixRcpp", quietly = TRUE)) {
    as.character(utils::packageVersion("FlowCodeUnmixRcpp"))
  } else {
    "0"
  }

  new.keywords <- utils::modifyList(new.keywords, list(
    "$FIL" = file.name,
    "$PAR" = as.character(n.param),
    "$TOT" = as.character(nrow(final.matrix)),
    "UNMIXINGMETHOD" = method,
    "$BYTEORD" = "1,2,3,4",
    "$DATATYPE" = "F",
    "SPECTRA" = format.matrix.string(spectra),
    "FLUOROCHROMES" = paste(rownames(spectra), collapse = ","),
    "AUTOSPECTRAL" = asp.ver,
    "AUTOSPECTRALRCPP" = asp.rcpp.ver,
    "FLOWCODEUNMIXRCPP" = fc.rcpp.ver,
    '$LAST_MODIFIED' = toupper(format(Sys.time(), "%d-%b-%Y %H:%M:%OS2")),
    '$LAST_MODIFIER' = sprintf("FlowCodeUnmixRcpp_%s", fc.rcpp.ver),
    '$ORIGINALITY' = "DataModified"
  ))

  if (!is.null(af.spectra)) {
    new.keywords[["AUTOFLUORESCENCE"]] <- format.matrix.string(af.spectra)
  }

  return(new.keywords)
}
