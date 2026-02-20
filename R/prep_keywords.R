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
    file.name
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

  # Build new parameter keywords
  param.keywords <- list()
  n.param <- ncol(final.matrix)
  p.names <- colnames(final.matrix)

  for (i in seq_len(n.param)) {
    p.name <- p.names[i]
    if (p.name %in% original.param) {
      old.entry <- param.lookup[[p.name]]
      if (!is.null(old.entry)) {
        names(old.entry) <- sub("^\\$?P\\d+", paste0("$P", i), names(old.entry))
        param.keywords <- c(param.keywords, old.entry)
      } else {
        param.keywords[[paste0("$P", i, "N")]] <- p.name
      }
    } else {
      bit.depth <- if (!is.null(asp$bit.depth)) asp$bit.depth else "32"
      param.keywords[[paste0("$P", i, "N")]] <- p.name
      param.keywords[[paste0("$P", i, "B")]] <- as.character(bit.depth)
      param.keywords[[paste0("$P", i, "E")]] <- "0,0"
      param.keywords[[paste0("$P", i, "R")]] <- as.character(asp$expr.data.max)
      param.keywords[[paste0("$P", i, "DISPLAY")]] <- ifelse(p.name == "AF.Index", "LIN", "LOG")
      param.keywords[[paste0("$P", i, "TYPE")]] <- ifelse(p.name == "AF.Index", "AF_Index", "Fluorescence")

      # Map Marker/Stain from flow.control
      clean.name <- sub("-A$", "", p.name)
      f.idx <- match(clean.name, flow.control$fluorophore)
      param.keywords[[paste0("$P", i, "S")]] <- if (!is.na(f.idx)) as.character(flow.control$antigen[f.idx]) else ""
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
  new.keywords <- utils::modifyList(new.keywords, list(
    "$FIL" = file.name,
    "$PAR" = as.character(n.param),
    "$TOT" = as.character(nrow(final.matrix)),
    "$UNMIXINGMETHOD" = method,
    "$BYTEORD" = "1,2,3,4",
    "$DATATYPE" = "F",
    "$SPECTRA" = format.matrix.string(spectra),
    "$FLUOROCHROMES" = paste(rownames(spectra), collapse = ","),
    "$AUTOSPECTRAL" = as.character(utils::packageVersion("AutoSpectral"))
  ))

  if (!is.null(af.spectra)) {
    new.keywords[["$AUTOFLUORESCENCE"]] <- format.matrix.string(af.spectra)
  }

  return(new.keywords)
}
