# writeFCS.r

#' @title Write FCS File
#'
#' @description
#' A light weight FCS file binary writer based on `flowstate`. Lower memory usage
#' and faster than `flowCore::write.FCS()`.
#'
#' @param mat The expression data, containing both unmixed data and any
#' retained parameters such as scatter and time.
#' @param keys The keywords to be written to the file.
#' @param file.name The name of the FCS file to be written.
#' @param output.dir A character string specifying the directory to save the
#' unmixed FCS file.
#'
#' @export

writeFCS <- function(mat, keys, file.name, output.dir) {
  delim <- "/"
  # Placeholder offsets
  keys[['$BEGINDATA']] <- "0000000000"
  keys[['$ENDDATA']]   <- "0000000000"

  text_seg <- paste0(delim, paste0(names(keys), delim, keys, delim, collapse = ""))

  # Calculate absolute byte positions
  t.start <- 58
  t.end   <- t.start + nchar(text_seg, "bytes") - 1
  d.start <- t.end + 1
  d.end   <- d.start + (nrow(mat) * ncol(mat) * 4) - 1

  # Patch the string
  text_seg <- sub("0000000000", sprintf("%10d", d.start), text_seg)
  text_seg <- sub("0000000000", sprintf("%10d", d.end), text_seg)

  # Header
  header <- sprintf("FCS3.1          58%8d%8d%8d       0       0", t.end, d.start, d.end)

  # Write
  con <- file(file.path(output.dir, file.name), open = "wb")
  writeBin(charToRaw(header), con)
  writeBin(charToRaw(text_seg), con)
  writeBin(as.vector(t(mat)), con, size = 4, endian = "little")
  close(con)
}
