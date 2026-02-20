# unmix_flowcode.r

#' @title Unmix FlowCode
#'
#' @description
#' Unmix FlowCode samples, correcting FRET errors and debarcoding the data.
#'
#' @importFrom parallelly availableCores
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#' @param spectra.variants Named list (names are fluorophores) carrying matrices
#' of spectral signature variations for each fluorophore. Prepare using
#' `get.spectral.variants()`. Default is `NULL`.
#' @param flowcode.spectra Structured output from `get.flowcode.spectra()`, which
#' details the combination-level spectral unmixing errors due to FRET-like
#' artefacts.
#' @param asp The AutoSpectral parameter list. Prepare using
#' `get.autospectral.param`.
#' @param thresholds Optional named numeric vector of positivity thresholds for
#' the FlowCode fluorophores. Overrides the thresholds provided by `flowcode.spectra`.
#' Default is `NULL`, which is unused.
#' @param k Numeric, controls the number of variants tested for each fluorophore,
#' autofluorescence and FRET spectrum. Default is `10`, which will be good, `1`
#' is fastest. Values up to `10` provide additional benefit in unmixing quality.
#'  Only used if `optimize=TRUE`.
#' @param parallel Logical, whether to use parallel processing for the per-cell
#' unmixing. Default is `FALSE`.
#' @param threads Numeric. Number of threads to use for parallel processing.
#' Defaults to `1` for sequential processing, or `0` (all cores) if `parallel=TRUE`.
#' @param verbose Logical, whether to send messages to the console.
#' Default is `TRUE`.
#' @param optimize Logical, whether to perform per-cell spectral optimization.
#' Faster without this, usually better with it.
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

unmix.flowcode <- function(
    raw.data,
    spectra,
    af.spectra,
    spectra.variants,
    flowcode.spectra,
    asp,
    thresholds = NULL,
    k = 1,
    parallel = TRUE,
    threads = if ( parallel ) 0 else 1,
    verbose = TRUE,
    optimize = TRUE
  ) {

  ### Setup-----------
  # check for AF in spectra, remove if present
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  fluorophores <- rownames( spectra )

  # check for fatal errors
  if ( nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided." )

  ### TBD implement more checks

  # define data structures
  cell.n <- nrow( raw.data )
  mean.af <- colMeans( af.spectra )
  af.spectra <- rbind( mean.af, af.spectra)
  fluorophore.n <- nrow( spectra )
  detector.n <- ncol( spectra )

  # set positivity thresholds vector
  pos.thresholds <- rep( Inf, fluorophore.n )
  names( pos.thresholds ) <- fluorophores
  # fill with data
  pos.thresholds[ names( spectra.variants$thresholds ) ] <- spectra.variants$thresholds

  # unpack spectral variants
  variants <- spectra.variants$variants
  delta.list <- spectra.variants$delta.list
  delta.norms <- spectra.variants$delta.norms

  if ( is.null( pos.thresholds ) )
    stop( "Check that spectral variants have been calculated using get.spectra.variants" )
  if ( !( length( variants ) > 1 ) )
    stop( "Multiple fluorophore spectral variants must be provided." )

  optimize.fluors <- fluorophores[ fluorophores %in% names( variants ) ]
  if ( !( length( optimize.fluors ) > 0 ) )
    stop( "No matching fluorophores between supplied spectra and spectral variants.
             No spectral optimization performed." )

  # unpack FlowCode FRET and thresholds
  flowcode.thresholds <- flowcode.spectra$Thresholds
  fret.spectra <- flowcode.spectra$FRET
  flowcode.fluors <- flowcode.spectra$Flowcode.fluors
  combo.df <- flowcode.spectra$Combos
  combo.n <- nrow( combo.df ) # note: combo.n + 1 = "other", 0 = "untransduced"
  flowcode.combo.logical <- flowcode.spectra$Logical.combo

  # overwrite thresholds if user has provided new ones
  if ( !is.null( thresholds ) )
    flowcode.thresholds <- thresholds

  if ( is.null( flowcode.thresholds ) )
    stop( "Check that FlowCode thresholds have been supplied correctly (use get.flowcode.spectra)." )
  if ( ! all( names( flowcode.thresholds ) %in% fluorophores ) )
    stop( "FlowCode thresholds don't match the fluorophore names supplied." )
  if ( !( length( fret.spectra ) > 1 ) )
    stop( "Multiple FRET options for FlowCode spectra must be provided." )

  # if delta.list and delta.norms are not provided by AutoSpectral (<v1.0.0), calculate
  # this can be done in a single lapply call
  if ( is.null( delta.list ) ) {
    # calculate deltas for each fluorophore's variants
    delta.list <- lapply( optimize.fluors, function( fl ) {
      variants[[ fl ]] - matrix(
        spectra[ fl, ],
        nrow = nrow( variants[[ fl ]] ),
        ncol = detector.n,
        byrow = TRUE
      )
    } )
    names( delta.list ) <- optimize.fluors

    # precompute delta norms
    delta.norms <- lapply( delta.list, function( d ) {
      sqrt( rowSums( d^2 ) )
    } )
  }


  ### single cell AF extraction------------
  if ( verbose ) message( "Determining initial AF assignments" )

  # calculate pseudoinverse
  unmixing.matrix <- solve.default( tcrossprod( spectra ), spectra )

  # initial unmix (no AF)
  unmixed <- raw.data %*% t( unmixing.matrix )

  result <- fit.af(
    raw.data,
    unmixed,
    unmixing.matrix,
    spectra,
    af.spectra
  )

  unmixed <- result$unmixed
  raw.data <- raw.data - result$fitted.af
  AF <- result$AF
  af.idx <- result$af.idx


  ### debarcoding-----------
  if ( verbose ) message( "Debarcoding FlowCodes..." )

  flowcode.ids <- debarcode(
    unmixed = unmixed,
    flowcode.fluors = flowcode.fluors,
    flowcode.thresholds = flowcode.thresholds,
    combo.df = combo.df
  )

  # set cells that have valid FlowCode combos for FRET correction
  flowcode.pos <- which( flowcode.ids %in% seq_len( combo.n ) )
  has.flowcode <- logical( cell.n )
  has.flowcode[ flowcode.pos ] <- TRUE


  ### FRET correction-----------
  if ( verbose ) message( "Correcting FRET..." )

  # return FRET indices as well
  result <- fit.fret(
    raw.data_subset = raw.data[ has.flowcode, , drop = FALSE ],
    unmixed_subset = unmixed[ has.flowcode, , drop = FALSE ],
    unmixing.matrix = unmixing.matrix,
    spectra = spectra,
    flowcode.ids_subset = flowcode.ids[has.flowcode],
    fret.spectra = fret.spectra,
    flowcode.fluors = flowcode.fluors,
    flowcode.combo.logical = flowcode.combo.logical
  )

  raw.data[ has.flowcode ] <- raw.data[ has.flowcode ] - result$fitted.fret
  unmixed[ has.flowcode, ] <- result$unmixed


  ### fluorophore spectral optimization-----------
  if ( optimize ) {
    if ( verbose ) message( "Optimizing unmix..." )

    # set number of threads to use
    if ( parallel ) {
      if ( is.null( threads ) ) threads <- asp$worker.process.n
      if ( threads == 0 ) threads <- parallelly::availableCores()
    } else {
      threads <- 1
    }

    unmixed <- optimize.flowcode(
      raw.data = raw.data,
      unmixed = unmixed,
      spectra = spectra,
      pos.thresholds = pos.thresholds,
      flowcode.ids = flowcode.ids,
      has.flowcode = has.flowcode,
      flowcode.combo.logical = flowcode.combo.logical,
      flowcode.fluors = flowcode.fluors,
      optimize.fluors = optimize.fluors,
      variants = variants,
      delta.list = delta.list,
      delta.norms = delta.norms,
      asp = asp,
      k = k,
      threads = threads,
      parallel = parallel
    )
  }


  ### create "FlowCode" channels-----------
  FlowCode <- rep( 0L, cell.n )

  # create empty matrix for assigning FlowCode combo data channels
  fc.data <- matrix(
    0L,
    nrow = cell.n,
    ncol = combo.n,
    dimnames = list( NULL, paste( "Tag:", combo.df$Id ) )
  )

  # loop over cells that have identified FlowCodes and assign them to "channels"
  for ( cell in which( has.flowcode ) ) {
    # which tag does this cell have?
    id <- flowcode.ids[ cell ]
    # FlowCode fluors for this Id
    fc.fluors <- flowcode.fluors[ flowcode.combo.logical[ id, ] == 1 ]
    # extract unmixed values for this cell
    vals <- unmixed[ cell, fc.fluors ]
    # assign median to the corresponding tag's channel and "FlowCode" channel
    FlowCode[ cell ] <- stats::median( vals )
    fc.data[ cell, id ] <- FlowCode[ cell ]
  }

  # column names for the output data
  extra.cols <- c( "AF", "AF Index", "FlowCode" )
  fc.colnames <- colnames( fc.data )
  final.colnames <- c( colnames( unmixed ), extra.cols, fc.colnames )

  # pre-allocate the matrix
  final.matrix <- matrix( 0, nrow = nrow( unmixed ), ncol = length( final.colnames ) )
  colnames( final.matrix ) <- final.colnames

  # fill by index
  unmixed.n <- ncol( unmixed )
  fc.n <- ncol( fc.data )

  final.matrix[ , 1:unmixed.n ] <- unmixed
  final.matrix[ , unmixed.n + 1 ] <- AF
  final.matrix[ , unmixed.n + 2 ] <- af.idx
  final.matrix[ , unmixed.n + 3 ] <- FlowCode
  final.matrix[ , ( unmixed.n + 4 ):ncol( final.matrix ) ] <- fc.data

  # clean up
  rm( unmixed, AF, af.idx, FlowCode, fc.data )

  return( unmixed )
}
