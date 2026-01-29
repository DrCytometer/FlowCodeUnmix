# debarcode.r

#' @title Debarcode
#'
#' @description
#' Debarcodes FlowCode data, identifying valid and invalid combinations.
#'
#' @param unmixed Matrix of unmixed flow cytometry expression data. Cells in
#' rows and fluorophores in columns. Fluorophores must include `flowcode.fluors`.
#' @param flowcode.fluors Named character vector of fluorophores used to identify
#' the FlowCodes. Names should be FlowCode epitope tags, values should be the
#' corresponding fluorophores.
#' @param flowcode.thresholds Named numeric vector of thresholds above which
#' the cell (event) may be considered positive for a given FlowCode. Names should
#' be fluorophores, values should be the threshold in raw, untransformed unmixed
#' scale.
#' @param combo.df Data frame describing the valid combinations of FlowCodes.
#' Structure: One row per combination. Columns are `Id`, `Procode.tag1`,
#' `Procode.tag2` and `Procode.tag3`, describing the name (e.g., CRISPR target),
#' and three epitopes for the combination, respectively.
#'
#' @return A numeric vector, length `nrow(unmixed)`, linking each event to an `Id`
#' in `combo.df` by row number.
#'
#' @export

debarcode <- function(
    unmixed,
    flowcode.fluors,
    flowcode.thresholds,
    combo.df
  ) {

  # set initial ID as missing
  flowcode.ids <- rep( NA_real_, nrow( unmixed ) )

  # determine which events are above the thresholds
  above.threshold <- sweep(
    unmixed[, flowcode.fluors ],
    2,
    flowcode.thresholds,
    ">"
  )
  above.threshold.count <- rowSums( above.threshold )

  # calculate how much real signal there is (for 3+ positive FlowCode correction)
  signal.over.threshold <- sweep(
    unmixed[ , flowcode.fluors ],
    2,
    flowcode.thresholds,
    FUN = "-"
  )

  # try to rectify >3 tagged events
  multi.idx <- which( above.threshold.count > 3 )

  for( i in multi.idx ) {
    x <- sort( signal.over.threshold[ i, ], decreasing = TRUE )
    if( x[ 3 ] > 2 * x[ 4 ] ) {
      above.threshold[ i, names( x )[ 4:length( x ) ] ] <- 0
    }
  }

  colnames( above.threshold ) <- names( flowcode.fluors )
  above.threshold.count <- rowSums( above.threshold )

  # only search triplet combos
  triplet.idx <- which( above.threshold.count == 3 )
  above.threshold <- above.threshold[ triplet.idx, , drop = FALSE ] > 0

  # generate tag combos
  hits <- which( above.threshold, arr.ind = TRUE )
  combo.list <- split( colnames( above.threshold )[ hits[ , 2 ] ], hits[ , 1 ] )
  flowcode.combos <- sapply( combo.list, function( x )
    paste( sort( toupper( x ) ), collapse = "_" ) )

  combo.df$FlowCode_Combination <- apply(
    combo.df, 1,
    function( x ) paste( sort( toupper( unlist( x[ -1 ] ) ) ), collapse = "_" )
  )

  # number of valid combinations
  flowcode.n <- length( unique( combo.df$FlowCode_Combination ) )

  flowcode.ids[ triplet.idx ] <- match( flowcode.combos, combo.df$FlowCode_Combination )

  # assign zero for untransduced cells
  flowcode.ids[ which( above.threshold.count == 0 ) ] <- 0

  # assign flowcode.n + 1 for cells with "other" combinations
  flowcode.ids[ is.na( flowcode.ids ) ] <- flowcode.n + 1

  # return ID numbers only
  return( flowcode.ids )
}
