# define_gate.r

# asp dependency could be removed:
# gate.data[ , 1 ], gate.data[,2]

#' @title Define Gate
#'
#' @description
#' Defines a gate boundary on forward (FSC) and side (SSC) scatter parameters. It
#' does so by selecting the brightest events in the FlowCode tag fluorescence
#' channels, determining their location on scatter, and defining a density/distance-
#' based region around the centroid of those events. The aim is to focus on the
#' scatter region with the highest signal, without any additional inputs or prior
#' knowledge, as means of determining where the cells are. This is important for
#' identifying positive events correctly, but also critical for matching the
#' background noise/autofluorescence/non-specific staining level on similar events
#' (cells rather than debris).
#'
#' @importFrom tripack tri.mesh convex.hull
#' @importFrom MASS cov.rob
#'
#' @param scatter.data Matrix of FSC and SSC data corresponding to the same
#' events in `unmixed.data`. Used to define the scatter gating region.
#' @param unmixed.data Matrix of unmixed flow cytometry expression data. Cells in
#' rows and fluorophores in columns. Fluorophores must include `flowcode.fluors`.
#' @param flowcode.fluors Named character vector of fluorophores used to identify
#' the FlowCodes. Names should be FlowCode epitope tags, values should be the
#' corresponding fluorophores.
#' @param asp The AutoSpectral parameter list.
#' @param n.cells Numeric, default is `200`. Number of brightest events per
#' channel to use in defining the gate region.
#'
#' @return Gate boundary
#'
#' @export

define.gate <- function(
    scatter.data,
    unmixed.data,
    flowcode.fluors,
    asp,
    n.cells = 200
  ) {

  # Check columns
  stopifnot( all( asp$default.scatter.parameter %in% colnames( scatter.data ) ) )

  # Keep only fluorophores present as columns
  fluor.cols <- intersect( flowcode.fluors, colnames( unmixed.data ) )
  if ( length( fluor.cols ) == 0 )
    stop( "None of flowcode.fluors found in scatter.data column names." )

  top.idx <- unique( unlist(
    lapply( fluor.cols, function( fluor ) {
      x <- unmixed.data[ , fluor ]
      if ( all( is.na( x ) ) ) return( NULL )
      order( x, decreasing = TRUE )[ 1:min( n.cells, length( x ) ) ]
    } )
  ) )

  fsc <- scatter.data[ top.idx, asp$default.scatter.parameter[ 1 ] ]
  ssc <- scatter.data[ top.idx, asp$default.scatter.parameter[ 2 ] ]

  scatter.coords <- cbind( fsc, ssc)

  # robust location and covariance
  cov.rob <- MASS::cov.rob( scatter.coords, method = "mcd" )

  center <- cov.rob$center
  covmat <- cov.rob$cov

  # squared Mahalanobis distance
  d2 <- stats::mahalanobis(
    scatter.coords,
    center = center,
    cov = covmat
  )

  # 70th percentile contour
  cutoff <- stats::quantile( d2, 0.70 )

  scatter.idx <- which( d2 <= cutoff )

  gate <- tripack::convex.hull(
    tripack::tri.mesh(
      fsc[ scatter.idx ],
      ssc[ scatter.idx ]
    )
  )

  return( gate )

}
