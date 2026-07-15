# unmix_flowcode_fcs.r

#' @title Spectral Voltage/Gain Keyword Suffix
#' @description
#' Internal helper. Determines which FCS keyword suffix ($PnV, $PnG, or
#' $PnR) should be compared as the "voltage" for a given cytometer.
#' - Mosaic uses $PnG (gain) instead of $PnV.
#' - ID7000 has no per-channel PMT voltage.
#' - All other supported cytometers use $PnV.
#' @param asp The AutoSpectral parameter list; `asp$cytometer` is used.
#' @param header Parsed FCS header/keyword list for the file being checked.
#' @keywords internal
.spectral.voltage.suffix <- function( asp, header ) {
  if ( grepl( "ID7000", asp$cytometer, ignore.case = TRUE ) ) return( NA_character_ )
  if ( grepl( "Mosaic", asp$cytometer, ignore.case = TRUE ) ) return( "G" )
  "V"
}

#' @title Extract Spectral Voltages/Gains
#' @description
#' Internal helper. Given an FCS header and a keyword suffix from
#' `.spectral.voltage.suffix()`, extracts the voltage/gain value for each
#' named spectral channel (matched by channel name, not position).
#' @param header Parsed FCS header/keyword list for the file being checked.
#' @param spectral.channel Character vector of spectral detector/channel names.
#' @param suffix Character keyword suffix (`"V"`, `"G"`, or `"R"`) from
#' `.spectral.voltage.suffix()`.
#' @keywords internal
.extract.spectral.voltages <- function( header, spectral.channel, suffix ) {
  if ( is.na( suffix ) ) {
    return( stats::setNames(
      rep( NA_character_, length( spectral.channel ) ), spectral.channel
    ) )
  }
  p.names <- unlist( header[ grep( "^\\$P\\d+N$", names( header ) ) ] )
  stats::setNames(
    vapply( spectral.channel, function( ch ) {
      p.idx.key <- names( p.names )[ which( p.names == ch ) ]
      if ( length( p.idx.key ) == 0 ) return( NA_character_ )
      n <- gsub( "[^0-9]", "", p.idx.key )
      val <- header[[ paste0( "$P", n, suffix ) ]]
      if ( is.null( val ) ) NA_character_ else as.character( val )
    }, character( 1 ) ),
    spectral.channel
  )
}

#' @title Unmix FlowCode FCS Data
#'
#' @description
#' This function performs spectral unmixing on FCS data from FlowCode samples.
#'
#' @param fcs.file A character string specifying the path to the FCS file.
#' @param spectra A matrix containing the spectral data.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param flow.control A list containing flow cytometry control parameters.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`. Required for `FlowCodeUnmix` unmixing.
#' @param spectra.variants Named list (names are fluorophores) carrying matrices
#' of spectral signature variations for each fluorophore. Prepare using
#' `get.spectral.variants`.
#' @param flowcode.spectra Structured output from `get.flowcode.spectra()`, which
#' details the combination-level spectral unmixing errors due to FRET-like
#' artefacts. Pass the result of `get.flowcode.spectra()` or use
#' `readRDS( "./flowcode_spectra/FlowCode_Spectra.rds" )` to read the file into R.
#' @param thresholds.file Optional, file name and path to the thresholds CSV file
#' produced using the threshold setting app. Thresholds are provided by default
#' as part of `flowcode.spectra`, but those will have been selected on the
#' FlowCode backbone sample. If the thresholds for positivity of any of the
#' fluorophores used to identify the FlowCodes differs in the fully stained
#' samples being unmixed here, provide a new set of thresholds from the
#' ThresholdApp. Run this on an unmixed sample (OLS from the instrument is fine),
#' and refer to the new thresholds file here. Default is `NULL`, which will be
#' ignored.
#' @param output.dir A character string specifying the directory to save the
#' unmixed FCS file. Default is `NULL`, which reverts to `asp$unmixed.fcs.dir`.
#' @param file.suffix A character string to append to the output file name.
#' Default is `NULL`.
#' @param include.imaging A logical value indicating whether to include imaging
#' parameters in the written FCS file. Default is `FALSE`.
#' @param parallel Logical, default is `TRUE`, which enables parallel processing
#' for per-cell unmixing methods.
#' @param threads Numeric, defaults to all available cores if `parallel=TRUE`.
#' @param verbose Logical, whether to send messages to the console.
#' Default is `TRUE`.
#' @param optimize Logical, whether to perform per-cell spectral optimization.
#' Faster without this, usually better with it.
#' @param n.passes Numeric, default `1`. Rounds of optimization to perform.
#' @param n.af.passes Numeric, default `1`. Rounds of AF extraction to perform.
#' @param refine.af.quantile Numeric in \[0, 1\], default `0.5`. Passed
#' through to `unmix.autospectral.rcpp()`.
#' @param cell.weight Logical, default depends on `asp$cytometer` (`TRUE` for
#' `"ID7000"`, else `FALSE`). Per-cell IRLS-style detector weighting, passed
#' through to the staged pipeline's `AutoSpectralRcpp::unmix.autospectral.rcpp()`
#' calls.
#' @param noise.floor Numeric, default `125`. Passed through to the staged
#' pipeline.
#' @param alpha Numeric in \[0, 1\], default `0.5`. Passed through to the
#' staged pipeline.
#' @param collinear.threshold Numeric in \[0, 1\], default `0.5`. Passed
#' through to the staged pipeline as `collinear.thresh`.
#' @param joint.pair.resolution Logical, default `TRUE`. Passed through to
#' the staged pipeline.
#' @param fret.alpha Numeric in \[0, 1\], default `0.5`. Passed through to
#' the staged pipeline's FRET fit (currently unused there).
#' @param fret.median.only Logical, default `FALSE`. If `TRUE`, skips the
#' per-cell FRET variant scan in the staged pipeline and always fits/subtracts
#' the median FRET row for each combo.
#' @param return.diagnostics Logical, default `FALSE`. If `TRUE`, per-chunk
#' diagnostics from the staged pipeline are collected and returned invisibly
#' as a list after the FCS file is written, instead of returning `NULL`.
#' Note `median.background` is pooled per chunk, not across the whole file.
#' @param chunk.size Numeric, number of events to use per chunk of unmixing. Used
#' to manage memory when processing large FCS files. As a rough guide, you will
#' need approximately 10x the size of the raw FCS file on disk as available
#' memory. Default is set at `2e6` events, assuming ~20GB memory available.
#'
#' @return None. The function writes the unmixed FCS data to a file.
#'
#' @export

unmix.flowcode.fcs <- function(
    fcs.file,
    spectra,
    asp,
    flow.control,
    af.spectra,
    spectra.variants,
    flowcode.spectra,
    thresholds.file = NULL,
    output.dir = NULL,
    file.suffix = NULL,
    include.imaging = FALSE,
    parallel = TRUE,
    threads = if ( parallel ) 0 else 1,
    verbose = TRUE,
    optimize = TRUE,
    n.passes = 1,
    n.af.passes = 1,
    refine.af.quantile = 0.5,
    cell.weight = if ( asp$cytometer == "ID7000" ) TRUE else FALSE,
    noise.floor = 125,
    alpha = 0.5,
    collinear.threshold = 0.5,
    joint.pair.resolution = TRUE,
    fret.alpha = 0.5,
    fret.median.only = FALSE,
    return.diagnostics = FALSE,
    chunk.size = 2e6
) {

  # check for AutoSpectral in NAMESPACE
  if ( !requireNamespace( "AutoSpectral", quietly = TRUE ) ) stop( "AutoSpectral required." )

  # unmix.flowcode.cpp.staged() hardcodes pipeline = "joint" internally when
  # calling AutoSpectralRcpp::unmix.autospectral.rcpp(), which requires
  # AutoSpectralRcpp >= 1.1.0
  if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) ) {
    if ( utils::packageVersion( "AutoSpectralRcpp" ) < package_version( "1.1.0" ) ) {
      stop(
        "FlowCodeUnmix's staged pipeline requires `AutoSpectralRcpp` >= 1.1.0, ",
        "but version ", utils::packageVersion( "AutoSpectralRcpp" ), " is installed.\n",
        "Please update AutoSpectralRcpp:\n",
        "  remotes::install_github(\"DrCytometer/AutoSpectralRcpp\")",
        call. = FALSE
      )
    }
  } else {
    warning(
      "Package `AutoSpectralRcpp` not found. The staged pipeline requires it; ",
      "unmixing will fall back to the (currently outdated) pure-R `unmix.flowcode()`.",
      call. = FALSE
    )
  }

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  # create the output folder if it doesn't exist
  output.dir <- if ( is.null( output.dir ) ) asp$unmixed.fcs.dir else output.dir
  if ( !dir.exists( output.dir ) ) dir.create( output.dir )

  #  import FCS
  if ( verbose ) message( "Reading FCS metadata: ", fcs.file )
  import.meta <- AutoSpectral::readFCS( fcs.file, return.keywords = TRUE, start.row = 1, end.row = 1 )
  fcs.keywords <- import.meta$keywords
  total.events <- as.numeric( fcs.keywords[[ "$TOT" ]] )
  original.param <- colnames( import.meta$data )

  # check for voltage/gain settings consistency with single-stained controls
  if ( !is.null( flow.control$voltages ) && !all( is.na( flow.control$voltages ) ) ) {
    all.keys <- names( fcs.keywords )
    p.name.keys <- grep( "^\\$P\\d+N$", all.keys, value = TRUE )
    p.names <- vapply( p.name.keys, function( k ) fcs.keywords[[ k ]], character( 1 ) )

    pnv.id <- .spectral.voltage.suffix( asp, fcs.keywords )

    for ( ch.name in names( flow.control$voltages ) ) {
      matching.key <- names( p.names )[ p.names == ch.name ]

      if ( length( matching.key ) > 0 ) {
        if ( is.na( pnv.id ) ) next  # no usable keyword for this file/cytometer

        v.key <- gsub( "N$", pnv.id, matching.key )
        current.v <- fcs.keywords[[ v.key ]]
        ref.v <- flow.control$voltages[[ ch.name ]]

        if ( !identical( as.character( current.v ), as.character( ref.v ) ) ) {
          warning( sprintf(
            "Voltage/gain mismatch for channel %s! Controls: %s, Current: %s. Unmixing may be inaccurate.",
            ch.name, ref.v, current.v
          ),
          call. = FALSE )
        }
      } else {
        stop( paste( "Spectral channel", ch.name, "not found in target FCS file." ),
              call. = FALSE )
      }
    }
  }

  # determine base filename
  method <- "AutoSpectral"
  file.name <- if ( !is.null( fcs.keywords$`$FIL` ) ) fcs.keywords$`$FIL` else basename( fcs.file )

  # handle manufacturer/cytometer specifics
  if ( asp$cytometer %in% c( "ID7000", "Mosaic" ) ) {
    # replace "Raw" with "AutoSpectral" (e.g., "Experiment_Raw.fcs" -> "Experiment_AutoSpectral.fcs")
    file.name <- sub("([ _])Raw(\\.fcs$|\\s|$)", paste0("\\1", method, "\\2"),
                     file.name, ignore.case = TRUE)

  } else if ( grepl("Discover", asp$cytometer ) ) {
    file.name <- if ( !is.null( fcs.keywords$FILENAME ) ) fcs.keywords$FILENAME else file.name
    file.name <- basename( file.name )
    file.name <- sub( "\\.fcs$", paste0(" ", method, ".fcs" ), file.name, ignore.case = TRUE)
  } else {
    # append method to filename
    file.name <- sub( "\\.fcs$", paste0(" ", method, ".fcs"), file.name, ignore.case = TRUE )
  }

  # add suffix if provided (before the .fcs extension)
  if ( !is.null( file.suffix ) ) {
    file.name <- sub( "\\.fcs$", paste0(" ", file.suffix, ".fcs"), file.name, ignore.case = TRUE )
  }

  # channel selection for unmixing
  spectral.channel <- colnames( spectra )
  other.channels <- setdiff( original.param, spectral.channel )

  # remove height and width if present
  for ( ch in spectral.channel[ grepl( "-A$", spectral.channel ) ] ) {
    base <- sub( "-A$", "", ch )
    other.channels <- setdiff( other.channels, paste0( base, c( "-H", "-W" ) ) )
  }
  if ( grepl( "Discover", asp$cytometer ) && !include.imaging ) {
    other.channels <- intersect( other.channels, asp$time.and.scatter )
  }

  # if user has provided a new threshold file, extract the thresholds
  if ( !is.null( thresholds.file ) ) {
    if ( verbose ) message( "Reading thresholds file" )
    flowcode.thresholds.df <- utils::read.csv( thresholds.file )
    flowcode.thresholds <- flowcode.thresholds.df$Threshold_Raw
    names( flowcode.thresholds ) <- flowcode.thresholds.df$Fluor
  } else {
    flowcode.thresholds <- NULL
  }

  # set multithreading
  if ( parallel & is.null( threads ) )
    threads <- asp$worker.process.n
  else if ( parallel & threads == 0 )
    threads <- parallelly::availableCores()


  ######### Unmixing ##############
  # unmix the data, correcting FRET, autofluorescence and spillover errors where possible

  # pre-allocate the full results matrix
  combo.n <- nrow( flowcode.spectra$Combos )
  cols.n <- length( other.channels ) + nrow( spectra ) + 3 + combo.n
  final.matrix <- matrix( 0, nrow = total.events, ncol = cols.n )

  chunk.n <- ceiling( total.events / chunk.size )
  diagnostics.list <- if ( isTRUE( return.diagnostics ) ) vector( "list", chunk.n ) else NULL

  # unmix in chunks for big files
  for ( i in 1:chunk.n ) {
    s.row <- ( (i - 1) * chunk.size ) + 1
    e.row <- min( i * chunk.size, total.events )

    if ( verbose )
      message( sprintf( "Processing chunk %d/%d (Events %d to %d)", i, chunk.n, s.row, e.row ) )

    # read in only events from this chunk
    chunk.data <- AutoSpectral::readFCS( fcs.file, return.keywords = FALSE, start.row = s.row, end.row = e.row )
    chunk.spectral <- chunk.data[ , spectral.channel, drop = FALSE ]
    chunk.other <- chunk.data[ , other.channels, drop = FALSE ]

    # Unmix Chunk
    # check whether FlowCodeUnmixRcpp in installed
    if ( requireNamespace("FlowCodeUnmixRcpp", quietly = TRUE ) &&
         "fit_flowcode_fret_lowd_cpp" %in% ls( getNamespace( "FlowCodeUnmixRcpp" ) ) ) {
      if ( verbose ) message( "Unmixing using FlowCodeUnmixRcpp" )
      chunk.result <- tryCatch(
        unmix.flowcode.cpp.staged(
          raw.data = chunk.spectral,
          spectra = spectra,
          af.spectra = af.spectra,
          spectra.variants = spectra.variants,
          flowcode.spectra = flowcode.spectra,
          asp = asp,
          thresholds = flowcode.thresholds,
          parallel = parallel,
          threads = threads,
          verbose = verbose,
          optimize = optimize,
          cell.weight = cell.weight,
          noise.floor = noise.floor,
          alpha = alpha,
          collinear.thresh = collinear.threshold,
          joint.pair.resolution = joint.pair.resolution,
          n.passes = n.passes,
          n.af.passes = n.af.passes,
          refine.af.quantile = refine.af.quantile,
          fret.alpha = fret.alpha,
          fret.median.only = fret.median.only,
          return.diagnostics = return.diagnostics
        ),
        error = function( e ) {
          warning(
            "FlowCodeUnmixRcpp unmixing failed, falling back to standard FlowCodeUnmix: ",
            e$message,
            call. = FALSE
          )
          if ( isTRUE( return.diagnostics ) )
            warning(
              "Diagnostics were requested but are not available from the pure-R fallback.",
              call. = FALSE
            )
          # NB: unmix.flowcode() does not accept cell.weight/noise.floor/alpha/
          # collinear.threshold/joint.pair.resolution/fret.alpha/fret.median.only
          # yet -- it predates the staged pipeline design. Revisit once the
          # pure-R fallback is rewritten to match.
          unmix.flowcode(
            raw.data = chunk.spectral,
            spectra = spectra,
            af.spectra = af.spectra,
            spectra.variants = spectra.variants,
            flowcode.spectra = flowcode.spectra,
            asp = asp,
            thresholds = flowcode.thresholds,
            parallel = parallel,
            threads = threads,
            verbose = verbose,
            optimize = optimize
          )
        }
      )
    } else {
      warning( "The FlowCodeUnmixRcpp package is not installed: unmixing will be slow.",
               call. = FALSE )
      if ( isTRUE( return.diagnostics ) )
        warning(
          "Diagnostics were requested but are not available from the pure-R fallback.",
          call. = FALSE
        )
      chunk.result <- unmix.flowcode(
        raw.data = chunk.spectral,
        spectra = spectra,
        af.spectra = af.spectra,
        spectra.variants = spectra.variants,
        flowcode.spectra = flowcode.spectra,
        asp = asp,
        thresholds = flowcode.thresholds,
        parallel = parallel,
        threads = threads,
        verbose = verbose,
        optimize = optimize
      )
    }

    # Normalize chunk.result: unmix.flowcode.cpp.staged() returns a plain
    # matrix unless return.diagnostics = TRUE, in which case it returns
    # list(unmixed = ..., diagnostics = ...). unmix.flowcode() (the pure-R
    # fallback) always returns a plain matrix.
    if ( isTRUE( return.diagnostics ) && is.list( chunk.result ) &&
         !is.null( chunk.result$unmixed ) ) {
      unmixed.chunk <- chunk.result$unmixed
      diagnostics.list[[ i ]] <- chunk.result$diagnostics
    } else {
      unmixed.chunk <- chunk.result
    }

    # Store in Pre-allocated Matrix
    final.matrix[ s.row:e.row, 1:ncol( chunk.other ) ] <- chunk.other
    final.matrix[ s.row:e.row, ( ncol( chunk.other ) + 1 ):cols.n ] <- unmixed.chunk

    # Cleanup memory immediately
    rm( chunk.data, chunk.spectral, chunk.other, unmixed.chunk )
    if (i %% 5 == 0) gc()
  }

  # fix any NA values (e.g., plate location with S8)
  if ( anyNA( final.matrix ) ) final.matrix[ is.na( final.matrix ) ] <- 0

  # set colnames
  combo.df <- flowcode.spectra$Combos
  unmixed.names <- c( rownames( spectra ), "AF", "AF Index", "FlowCode", combo.df$Id )
  colnames( final.matrix ) <- c( other.channels, unmixed.names )

  # append "-A" to fluorophore and AF channel names
  colnames( final.matrix ) <- ifelse(
    colnames( final.matrix ) %in% c( rownames( spectra ), "AF" ),
    paste0( colnames( final.matrix ), "-A" ),
    colnames( final.matrix )
  )

  # update keywords
  new.keywords <- prep.keywords(
    fcs.keywords = fcs.keywords,
    final.matrix = final.matrix,
    original.param = original.param,
    spectra = spectra,
    af.spectra = af.spectra,
    flow.control = flow.control,
    asp = asp,
    method = method,
    file.name = file.name,
    combos = combo.df$Id,
    include.imaging = include.imaging
  )

  # save file ---------
  if ( verbose ) message( paste( "Writing:", file.name ) )
  AutoSpectral::writeFCS( final.matrix, new.keywords, file.name, output.dir )

  if ( isTRUE( return.diagnostics ) ) {
    keep <- !vapply( diagnostics.list, is.null, logical( 1 ) )
    combined.diagnostics <- list(
      # pooled per chunk, not across the whole file -- see @param return.diagnostics
      median.background = lapply( diagnostics.list[ keep ], `[[`, "median.background" ),
      flowcode.ids  = do.call( c, lapply( diagnostics.list[ keep ], `[[`, "flowcode.ids" ) ),
      fret.k        = do.call( c, lapply( diagnostics.list[ keep ], `[[`, "fret.k" ) ),
      fret.index    = do.call( c, lapply( diagnostics.list[ keep ], `[[`, "fret.index" ) ),
      resid.ratio   = do.call( c, lapply( diagnostics.list[ keep ], `[[`, "resid.ratio" ) ),
      leakage.ratio = do.call( c, lapply( diagnostics.list[ keep ], `[[`, "leakage.ratio" ) )
    )
    return( invisible( combined.diagnostics ) )
  }

  invisible( NULL )
}



