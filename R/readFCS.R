# readFCS.r

#' @title Read FCS File
#'
#' @description
#' A light weight FCS file binary reader based on `flowstate`. Lower memory usage
#' and faster than `flowCore::read.FCS()`.
#'
#' @param fcs.path A character string specifying the file path (directory and
#' file name) for the .fcs file to be read.
#' @param keywords Logical, default `FALSE`. Controls whether keywords are
#' returned as well as the expression data matrix (useful for writing a new file).
#'
#' @return If `keywords = TRUE`, a list containing two elements:
#' 1) the expression data in a matrix, and 2) the keywords.
#' If `keywords = FALSE`, only the expression data is returned as a matrix.
#'
#' @export

readFCS <- function(fcs.path, keywords = FALSE) {
  con <- file(fcs.path, open = "rb")
  on.exit(close(con))

  # Header & Text Offsets
  header <- readChar(con, 58)
  txt.st <- as.numeric(trimws(substr(header, 11, 18)))
  txt.en <- as.numeric(trimws(substr(header, 19, 26)))

  # Extract Keywords
  seek(con, txt.st)
  txt_raw <- readBin(con, "raw", n = txt.en - txt.st + 1)
  txt <- rawToChar(txt_raw[txt_raw != as.raw(0)])

  delim <- substr(txt, 1, 1)
  kv <- strsplit(substr(txt, 2, nchar(txt)), delim, fixed = TRUE)[[1]]
  keywords <- stats::setNames(as.list(kv[seq(2, length(kv), by = 2)]), kv[seq(1, length(kv), by = 2)])

  # Data Offsets & Metadata
  data.st <- as.numeric(keywords[["$BEGINDATA"]])
  data.en <- as.numeric(keywords[["$ENDDATA"]])
  n.par <- as.numeric(keywords[["$PAR"]])
  endian <- if (keywords[["$BYTEORD"]] == "1,2,3,4") "little" else "big"

  # Read Matrix
  seek(con, data.st)
  n.vals <- (data.en - data.st + 1) / 4

  data.mat <- matrix(
    readBin(con, "numeric", n = n.vals, size = 4, endian = endian),
    ncol = n.par, byrow = TRUE
  )

  colnames(data.mat) <- unname(sapply(1:n.par, function(i) keywords[[paste0("$P", i, "N")]]))
  return(list(data = data.mat, keywords = keywords))
}
