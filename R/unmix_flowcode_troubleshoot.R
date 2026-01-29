# unmix_flowcode.r

#' @title Unmix FlowCode
#'
#' @description
#' Unmix FlowCode samples, correcting FRET errors and debarcoding the data.
#'
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
#' @param asp description
#' @param weights Optional numeric vector of weights (one per fluorescent
#' detector). Default is `NULL`, in which case weighting will be done by
#' channel means (Poisson variance).
#' @param cell.weighting Logical, whether to use cell-specific weighting for a
#' more Poisson-like unmixing. Default is `FALSE`.
#' @param cell.weight.regularize Logical, whether to regularize cell-specific
#' weights towards the bulk mean weighting set by `weights`. 50:50 averaging.
#' Default is `FALSE`.
#' @param k Numeric, controls the number of variants tested for each fluorophore,
#' autofluorescence and FRET spectrum. Default is `1`, which will be fastest.
#' Values up to `10` provide additional benefit in unmixing quality.
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

unmix.flowcode.test <- function(
    raw.data,
    spectra,
    af.spectra,
    spectra.variants,
    flowcode.spectra,
    asp,
    weights = NULL,
    cell.weighting = FALSE,
    cell.weight.regularize = FALSE,
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

  # check for AF in spectra, remove if present
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  fluorophores <- rownames( spectra )

  # check for fatal errors
  if ( nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided." )

  # define data structures
  cell.n <- nrow( raw.data )
  af.n <- nrow( af.spectra )
  mean.af <- colMeans( af.spectra )
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

  # unpack spectral variants
  pos.thresholds <- spectra.variants$thresholds
  pos.thresholds <- c( pos.thresholds, AF = -Inf )
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

  if ( is.null( flowcode.thresholds ) )
    stop( "Check that FlowCode thresholds have been supplied correctly (use get.flowcode.spectra)." )
  if ( ! all( names( flowcode.thresholds ) %in% fluorophores ) )
    stop( "FlowCOde thresholds don't match the fluorophore names supplied." )
  if ( !( length( combo.fret ) > 1 ) )
    stop( "Multiple FRET options for FlowCode spectra must be provided." )
  ## more checks to be implemented here
  # check for match between flowcode fluors and `fluorophores`


  # if delta.list and delta.norms are not provided by AutoSpectral (<v1.0.0), calculate
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

  # set baseline for cell-specific weighting
  cell.weight <- weights



  # loop over each cell, optimizing fluorophore spectra
  # also extracting FRET errors if the cell has a FlowCode combination

  for ( cell in seq_len( cell.n ) ) {
    # get cell's data
    cell.raw <- raw.data[ cell, , drop = FALSE ]
    cell.unmixed <- unmixed[ cell, fluors.af, drop = FALSE ]

    # check whether this cell has any fluorophores present
    # if it does not, we skip the rest of the processing
    # identify which fluorophores are present on this cell
    pos.fluors <- stats::setNames(
      as.vector( cell.unmixed >= pos.thresholds[ colnames( cell.unmixed ) ] ),
      colnames( cell.unmixed )
    )

    # if nothing is positive, everything is AF
    # assign spectrum accordingly and re-unmix
    if ( !any( pos.fluors ) ) {
      cell.spectra.final[ "AF", ] <- cell.raw / max( abs( cell.raw ) ) # should safeguard against 0
      cell.unmixed <- unmix( cell.raw, cell.spectra.final, cell.weights )

      return( cell.unmixed )
    }

    # otherwise, proceed
    pos.fluor.names <- names( pos.fluors )[ pos.fluors ]

    # if the cell does have fluorophores, determine which AF has been selected
    cell.af.idx <- unmixed[ cell, "AF Index" ]

    # if no AF was assigned, use mean as starting point for search
    if ( cell.af.idx == 0 ) {
      cell.af <- mean.af
    } else {
      cell.af <- af.spectra[ cell.af.idx, ]
    }

    # calculate delta matrix and norm for this AF spectrum
    af.delta <- af.spectra - matrix(
      cell.af,
      nrow = af.n,
      ncol = detector.n,
      byrow = TRUE
    )
    af.delta.norm <- sqrt( rowSums( af.delta^2 ) )

    # set weights in a cell-specific manner
    cell.weights <- weights
    if ( cell.weighting ) {
      # use cell-specific weighting (Poisson-like)
      cell.weights <- abs( cell.raw )
      cell.weights[ cell.weights < 1e-6 ] <- 1e-6
      cell.weights <- 1 / cell.weights

      if ( cell.weight.regularize ) {
        # regularize weight towards weights for full data
        cell.weights <- ( cell.weights + weights ) / 2
      }
    }

    # set baseline spectra
    cell.spectra.final <- combined.spectra
    cell.spectra.final[ "AF", ] <- cell.af
    cell.spectra.curr <- cell.spectra.final[ pos.fluors, , drop = FALSE ]
    cell.unmixed <- unmix( cell.raw, cell.spectra.curr, cell.weights )

    # set baseline unmixed and residuals
    fitted <- cell.unmixed %*% cell.spectra.curr
    resid <- cell.raw - fitted
    error.final <- sum( abs( resid ) )




    ### check whether this cell has a valid FlowCode combo
    # if it does, remove any non-combo FlowCode fluors (incorrect) for the
    # optimization steps. These will be allowed in later.
    # does this cell have a FlowCode combo?
    if ( has.flowcode[ cell ] ) {
      # which FlowCode combo is this?
      id <- flowcode.ids[ cell ]

      # empty FRET vector
      fitted.fret <- 0

      # use alignment between delta and residuals to select best FRET spectrum
      variants.fr <- combo.fret[[ id ]]
      delta.fr <- fret.delta.list[[ id ]]
      delta.norm <- fret.delta.norms[[ id ]]

      # remove non-combo FlowCode fluors from those considered present
      non.combo.fluors <- colnames( flowcode.combo.logical )[ which(
        flowcode.combo.logical[ id, ] == 0 ) ]
      pos.fluor.names <- pos.fluor.names[ !( pos.fluor.names %in% non.combo.fluors ) ]

      # drop non-combo FlowCode fluors from spectra for testing
      cell.spectra.curr <- cell.spectra.final[ pos.fluor.names, , drop = FALSE ]

      # unmix with median FRET spectrum
      cell.spectra.fret <- rbind(
        variants.fr[ 1, ], # first row is median
        cell.spectra.curr
      )
      rownames( cell.spectra.fret ) <- c( "FRET", rownames( cell.spectra.curr ) )

      trial.unmix <-  unmix(
        cell.raw,
        cell.spectra.fret,
        cell.weights
      )

      # check if this lowers the residual
      trial.resid <- cell.raw - ( trial.unmix %*% cell.spectra.fret )
      trial.error <- sum( abs( trial.resid ) )

      # set this as the target to beat
      resid <- trial.resid
      error.final <- trial.error
      fitted.fret <- trial.unmix[ , "FRET", drop = FALSE ] %*%
        cell.spectra.fret[ "FRET", , drop = FALSE ]

      # score FRET variants based on alignment to the residual
      joint.score <- as.numeric( delta.fr %*% t( resid ) ) * trial.unmix[ , "FRET" ]
      joint.score <- joint.score / delta.norm
      resid.norm <- sqrt( sum( resid^2 ) )
      joint.score <- joint.score / resid.norm

      # select number of high-scoring variants to test
      k.eff <- min( k, length( joint.score ) )
      topK <- order( joint.score, decreasing = TRUE )[ seq_len( k.eff ) ]

      # test the top k scoring variants sequentially
      for ( var in topK ) {
        # supplant the base spectrum with this variant
        cell.spectra.fret[ "FRET", ] <- variants.fr[ var, ]

        # reunmix with this variant
        trial.unmix <- unmix(
          cell.raw,
          cell.spectra.fret,
          cell.weights
        )

        # assess the residual error with this variant
        trial.resid <- cell.raw - ( trial.unmix %*% cell.spectra.fret )
        trial.error <- sum( abs( trial.resid ) )

        # accept change if residual is lower
        if ( trial.error < error.final ) {
          error.final <- trial.error
          resid <- trial.resid
          fitted.fret <- trial.unmix[ , "FRET", drop = FALSE ] %*%
            cell.spectra.fret[ "FRET", , drop = FALSE ]
        }
      }

      # subtract FRET from cell.raw for fluorophore optimization
      cell.raw <- cell.raw - fitted.fret

      # recompute baseline unmix and residuals after FRET correction
      cell.unmixed <- unmix( cell.raw, cell.spectra.curr, cell.weights )
      fitted <- cell.unmixed %*% cell.spectra.curr
      resid <- cell.raw - fitted
      error.final <- sum( abs( resid ) )
    }




    ### assess AF and update
    joint.score <- as.numeric( af.delta %*% t( resid ) ) * cell.unmixed[ , "AF" ]
    joint.score <- joint.score / af.delta.norm
    resid.norm <- sqrt( sum( resid^2 ) )
    joint.score <- joint.score / resid.norm

    k.eff <- min( k, length( joint.score ) )

    topK <- order( joint.score, decreasing = TRUE )[ seq_len( k.eff ) ]

    # test the top k scoring variants
    for ( var in topK ) {
      # supplant the base spectrum with this variant
      cell.spectra.curr[ "AF", ] <- af.spectra[ var, ]

      # reunmix with this variant
      trial.unmix <- unmix(
        cell.raw,
        cell.spectra.curr,
        cell.weights
      )

      # assess the residual error with this variant
      trial.resid <- cell.raw - ( trial.unmix %*% cell.spectra.curr )
      trial.error <- sum( abs( trial.resid ) )

      # accept change if residual is lower
      if ( trial.error < error.final ) {
        error.final <- trial.error
        cell.spectra.final[ "AF", ] <- cell.spectra.curr[ "AF", ]
        resid <- trial.resid
      } else {
        # reject if not
        cell.spectra.curr[ "AF", ] <- cell.spectra.final[ "AF", ]
      }
    }

    # recompute baseline unmix after AF correction
    cell.unmixed <- unmix( cell.raw, cell.spectra.final, cell.weights )

    # reassess fluorophore positivity--this may not be really necessary
    pos.fluors <- stats::setNames(
      as.vector( cell.unmixed >= pos.thresholds[ colnames( cell.unmixed ) ] ),
      colnames( cell.unmixed )
    )

    # if it is necessary, we can again drop out cells with only AF remaining
    if ( !any( pos.fluors ) ) {
      cell.spectra.final[ "AF", ] <- cell.raw / max( abs( cell.raw ) ) # should safeguard against 0
      cell.unmixed <- unmix( cell.raw, cell.spectra.final, cell.weights )

      return( cell.unmixed )
    }

    pos.fluor.names <- names( pos.fluors )[ pos.fluors ]

    # reset baseline spectra
    cell.spectra.curr <- cell.spectra.final[ pos.fluors, , drop = FALSE ]
    cell.unmixed <- unmix( cell.raw, cell.spectra.curr, cell.weights )
    fitted <- cell.unmixed %*% cell.spectra.curr
    resid <- cell.raw - fitted
    error.final <- sum( abs( resid ) )





    ### now optimize fluorophores
    # optimize variants
    fluors.to.sort <- optimize.fluors[ optimize.fluors %in% pos.fluor.names ]
    if ( length( fluors.to.sort ) > 0 ) {
      # sort by abundance to optimize brightest fluors first (error is proportional to signal)
      fluor.order <- sort( cell.unmixed[ , fluors.to.sort ], decreasing = TRUE )

      for ( fl in names( fluor.order ) ) {
        fl.variants <- variants[[ fl ]]
        delta.fl <- delta.list[[ fl ]]
        delta.norm  <- delta.norms[[ fl ]]

        # score variants
        joint.score <- as.numeric( delta.fl %*% t( resid ) ) * cell.unmixed[ , fl ]
        joint.score <- joint.score / delta.norm
        resid.norm <- sqrt( sum( resid^2 ) )
        joint.score <- joint.score / resid.norm

        k.eff <- min( k, length( joint.score ) )

        topK <- order( joint.score, decreasing = TRUE )[ seq_len( k.eff ) ]

        # test the top k scoring variants
        for ( var in topK ) {
          # supplant the base spectrum with this variant
          cell.spectra.curr[ fl, ] <- fl.variants[ var, ]

          # reunmix with this variant
          trial.unmix <- unmix(
            cell.raw,
            cell.spectra.curr,
            cell.weights
          )
          # assess the residual error with this variant
          trial.resid <- cell.raw - ( trial.unmix %*% cell.spectra.curr )
          trial.error <- sum( abs( trial.resid ) )

          # accept change if residual is lower
          if ( trial.error < error.final ) {
            error.final <- trial.error
            cell.spectra.final[ fl, ] <- cell.spectra.curr[ fl, ]
            resid <- trial.resid
          } else {
            # reject if not
            cell.spectra.curr[ fl, ] <- cell.spectra.final[ fl, ]
          }
        }
      }
    }

    # final unmix using optimized spectra
    unmixed[ cell, fluors.af ] <- unmix(
      cell.raw,
      cell.spectra.final,
      cell.weights
    )
  }




  return( unmixed )
}
