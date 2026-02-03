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
#' @param weights Optional numeric vector of weights (one per fluorescent
#' detector). Default is `NULL`, in which case weighting will be done by
#' channel means (Poisson variance).
#' @param cell.weighting Logical, whether to use cell-specific weighting for a
#' more Poisson-like unmixing. Default is `FALSE`.
#' @param cell.weight.regularize Logical, whether to regularize cell-specific
#' weights towards the bulk mean weighting set by `weights`. 50:50 averaging.
#' Default is `TRUE`. Only active if `cell.weighting=TRUE`.
#' @param k Numeric, controls the number of variants tested for each fluorophore,
#' autofluorescence and FRET spectrum. Default is `10`, which will be good, `1`
#' is fastest. Values up to `10` provide additional benefit in unmixing quality.
#' @param parallel Logical, whether to use parallel processing for the per-cell
#' unmixing. Default is `FALSE`.
#' @param threads Numeric. Number of threads to use for parallel processing.
#' Defaults to `1` for sequential processing, or `0` (all cores) if `parallel=TRUE`.
#' @param verbose Logical, whether to send messages to the console.
#' Default is `TRUE`.
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
    weights = NULL,
    cell.weighting = FALSE,
    cell.weight.regularize = TRUE,
    k = 1,
    parallel = TRUE,
    threads = if ( parallel ) 0 else 1,
    verbose = TRUE
  ) {

  # check for AutoSpectral in NAMESPACE
  if ( !requireNamespace( "AutoSpectral", quietly = TRUE ) ) {
    stop(
      "The AutoSpectral package is required but is not installed or not available.",
      call. = FALSE
    )
  }

  # check for FlowCodeUnmixRcpp in NAMESPACE
  if ( !requireNamespace( "FlowCodeUnmixRcpp", quietly = TRUE ) ) {
    warning(
      "The FlowCodeUnmixRcpp package is not installed--unmixing will be slow.",
      call. = FALSE
    )
  }

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
  af.n <- nrow( af.spectra )
  fluorophore.n <- nrow( spectra )
  detector.n <- ncol( spectra )
  combined.spectra <- matrix(
    NA_real_,
    nrow = fluorophore.n + 1,
    ncol = detector.n
  )
  colnames( combined.spectra ) <- colnames( spectra )
  fluors.af <- c( fluorophores, "AF" )
  rownames( combined.spectra ) <- fluors.af
  combined.spectra[ 1:fluorophore.n, ] <- spectra
  af.idx.in.spectra <- which( rownames( combined.spectra ) == "AF" ) # can also just be fluorophore.n + 1

  # set positivity thresholds vector
  pos.thresholds <- rep( Inf, fluorophore.n + 1 )
  names( pos.thresholds ) <- fluors.af
  pos.thresholds[ "AF" ] <- -Inf
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
  combo.fret <- flowcode.spectra$FRET
  flowcode.fluors <- flowcode.spectra$Flowcode.fluors
  combo.df <- flowcode.spectra$Combos
  combo.n <- nrow( combo.df ) # note: combo.n + 1 = "other", 0 = "untransduced"
  fret.delta.list <- flowcode.spectra$Delta
  fret.delta.norms <- flowcode.spectra$Delta.norms
  flowcode.combo.logical <- flowcode.spectra$Logical.combo

  # overwrite thresholds if user has provided new ones
  if ( !is.null( thresholds ) )
    flowcode.thresholds <- thresholds

  if ( is.null( flowcode.thresholds ) )
    stop( "Check that FlowCode thresholds have been supplied correctly (use get.flowcode.spectra)." )
  if ( ! all( names( flowcode.thresholds ) %in% fluorophores ) )
    stop( "FlowCode thresholds don't match the fluorophore names supplied." )
  if ( !( length( combo.fret ) > 1 ) )
    stop( "Multiple FRET options for FlowCode spectra must be provided." )
  ## more checks to be implemented here
  # check for match between flowcode fluors and `fluorophores`


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

  # set unmixing algorithm
  unmix <- AutoSpectral::unmix.wls.fast

  if ( is.null( weights ) ) {
    weights <- abs( colMeans( raw.data ) )
    weights[ weights < 1e-6 ] <- 1e-6
    weights <- 1 / weights
  }

  # initial unmixing without any AF--------
  if ( verbose ) message( "Initializing unmix" )

  unmixed <- unmix( raw.data, spectra, weights )

  # calculate initial fluorophore signal as error
  error <- rowSums( abs( unmixed[ , fluorophores, drop = FALSE ] ) )

  # set-up AF column, set-up AF Index
  initial.af <- matrix( 0, nrow = cell.n, ncol = 2 )
  colnames( initial.af ) <- c( "AF", "AF Index" )
  unmixed <- cbind( unmixed, initial.af )

  ### TBD: test scoring approach for speed-up

  ### single cell AF extraction------------
  if ( verbose ) message( "Extracting AF cell-by-cell..." )

  for ( af in seq_len( af.n ) ) {

    # set this AF as the spectrum to use
    combined.spectra[ fluorophore.n + 1, ] <- af.spectra[ af, , drop = FALSE ]

    # unmix with this AF
    unmixed.af <- unmix( raw.data, combined.spectra, weights )

    error.af <- rowSums( abs( unmixed.af[ , fluorophores, drop = FALSE ] ) )
    improved <- which( error.af < error )

    # track improvements
    if ( length( improved ) > 0 ) {
      # update error and unmixed data for improved cells
      error[ improved ] <- error.af[ improved ]
      unmixed[ improved, fluors.af ] <- unmixed.af[ improved, ]
      unmixed[ improved, "AF Index" ] <- af
    }
  }

  # assign median AF for cells that are unassigned
  zero.af <- which( unmixed[ , "AF Index" ] == 0 )
  combined.spectra[ af.idx.in.spectra, ] <- mean.af
  unmixed[ zero.af, fluors.af ] <- unmix(
    raw.data[ zero.af, , drop = FALSE ],
    combined.spectra,
    weights
  )
  unmixed[ zero.af, "AF Index" ] <- 1

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

  ### per cell fluorophore optimization with FRET correction------------
  if ( verbose ) message( "Optimizing fluorophore unmixing cell-by-cell..." )

  # split out the AF tracking from the unmixed data
  af.idx <- unmixed[ , "AF Index" ]

  # pre-calculate common indices
  fluors.af <- which( rownames( combined.spectra ) == fluors.af )
  fluorophores <- which( rownames( combined.spectra ) %in% fluorophores )

  # set number of threads to use
  if ( parallel ) {
    if ( is.null( threads ) ) threads <- asp$worker.process.n
    if ( threads == 0 ) threads <- parallelly::availableCores()
  } else {
    threads <- 1
  }

  # Use C++ for per-cell unmixing if FlowCodeUnmixRcpp is installed
  if ( requireNamespace( "FlowCodeUnmixRcpp", quietly = TRUE ) ) {
    if ( verbose ) message( "Using FlowCodeUnmixRcpp" )

    flowcode.ids.cpp <- flowcode.ids
    flowcode.ids.cpp[ !has.flowcode ] <- 1L
    af_idx_cpp <- af.idx.in.spectra - 1L

    unmixed[ , fluors.af ] <- tryCatch(
      FlowCodeUnmixRcpp::optimize_flowcode_unmix(
        raw_data = raw.data,
        unmixed = unmixed[ , fluors.af, drop = FALSE ],
        combined_spectra = combined.spectra,
        weights = weights,
        pos_thresholds = pos.thresholds,
        af_idx = af.idx,
        af_spectra = af.spectra,
        flowcode_ids = flowcode.ids.cpp,
        has_flowcode = has.flowcode,
        combo_fret = combo.fret,
        fret_delta_list = fret.delta.list,
        fret_delta_norms = fret.delta.norms,
        flowcode_combo_logical = flowcode.combo.logical,
        flowcode_fluors = flowcode.fluors,
        optimize_fluors = optimize.fluors,
        variants = variants,
        delta_list = delta.list,
        delta_norms = delta.norms,
        all_fluor_names = rownames( combined.spectra ),
        af_idx_in_spectra = af_idx_cpp,
        k = k,
        weighted = TRUE,
        cell_weighting = cell.weighting,
        cell_weight_regularize = cell.weight.regularize,
        nthreads = threads
      ),
      error = function( e ) {
        warning(
          "FlowCodeUnmixRcpp failed, falling back to standard FlowCodeUnmix: ",
          e$message,
          call. = FALSE
        )

        # set unmixing algorithm
        unmix <- FlowCodeUnmix::unmix.wls.fast

        optimize.flowcode(
          raw_data = raw.data,
          unmixed = unmixed[ , fluors.af, drop = FALSE ],
          unmix_fun = unmix,
          combined_spectra = combined.spectra,
          weights = weights,
          pos_thresholds = pos.thresholds,
          af_idx = af.idx,
          af_spectra = af.spectra,
          flowcode_ids = flowcode.ids,
          has_flowcode = has.flowcode,
          combo_fret = combo.fret,
          fret_delta_list = fret.delta.list,
          fret_delta_norms = fret.delta.norms,
          flowcode_combo_logical = flowcode.combo.logical,
          flowcode_fluors = flowcode.fluors,
          optimize_fluors = optimize.fluors,
          variants = variants,
          delta_list = delta.list,
          delta_norms = delta.norms,
          fluorophores = fluorophores,
          af_idx_in_spectra = af.idx.in.spectra,
          k = k,
          cell_weighting = cell.weighting,
          cell_weight_regularize = cell.weight.regularize,
          nthreads = threads,
          parallel = parallel
        )
      }
    )

  } else {
    # fall back to slow R-based unmixing
    if ( verbose )
      message( "FlowCodeUnmixRcpp unavailable, falling back to standard FlowCodeUnmix" )

    # set unmixing algorithm
    unmix <- FlowCodeUnmix::unmix.wls.fast

    unmixed[ , fluors.af ] <- optimize.flowcode(
      raw_data = raw.data,
      unmixed = unmixed[ , fluors.af, drop = FALSE ],
      unmix_fun = unmix,
      combined_spectra = combined.spectra,
      weights = weights,
      pos_thresholds = pos.thresholds,
      af_idx = af.idx,
      af_spectra = af.spectra,
      flowcode_ids = flowcode.ids,
      has_flowcode = has.flowcode,
      combo_fret = combo.fret,
      fret_delta_list = fret.delta.list,
      fret_delta_norms = fret.delta.norms,
      flowcode_combo_logical = flowcode.combo.logical,
      flowcode_fluors = flowcode.fluors,
      optimize_fluors = optimize.fluors,
      variants = variants,
      delta_list = delta.list,
      delta_norms = delta.norms,
      fluorophores = fluorophores,
      af_idx_in_spectra = af.idx.in.spectra,
      k = k,
      cell_weighting = cell.weighting,
      cell_weight_regularize = cell.weight.regularize,
      nthreads = threads,
      parallel = parallel
    )
  }


  # create "FlowCode" channel--for now, just median of expression in the three tag channels
  ### later, use unmixing of normalized combo + FRET spectrum to unmix this per cell
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
    fc.fluors <- flowcode.fluors[
      flowcode.combo.logical[ id, ] == 1
    ]

    # extract unmixed values for this cell
    vals <- unmixed[ cell, fc.fluors ]

    # assign median to the corresponding tag's channel and "FlowCode" channel
    FlowCode[ cell ] <- stats::median( vals )
    fc.data[ cell, id ] <- FlowCode[ cell ]

    ### Note: The "FlowCode" channel should be generated by unmixing each cell's
    # combined 3 tag + FRET signature in the low dimensional space
    # This will need to be left to later.
  }

  unmixed <- cbind( unmixed, FlowCode, fc.data )

  return( unmixed )
}
