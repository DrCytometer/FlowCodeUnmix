# unmix_flowcode_fcs.r

#' @title Unmix FlowCode FCS Data
#'
#' @description
#' This function performs spectral unmixing on FCS data from FlowCode samples.
#'
#' @importFrom flowCore read.FCS keyword exprs flowFrame parameters
#' @importFrom flowCore write.FCS parameters<- keyword<-
#'
#' @param fcs.file A character string specifying the path to the FCS file.
#' @param spectra A matrix containing the spectral data.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param flow.control A list containing flow cytometry control parameters.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`. Required for `FlowCodeUnmixe` unmixing.
#' @param spectra.variants Named list (names are fluorophores) carrying matrices
#' of spectral signature variations for each fluorophore. Prepare using
#' `get.spectral.variants`.
#' @param flowcode.spectra Structured output from `get.flowcode.spectra()`, which
#' details the combination-level spectral unmixing errors due to FRET-like
#' artefacts.
#' @param thresholds.file Optional, file name and path to the thresholds CSV file
#' produced using the threshold setting app. Thresholds are provided by default
#' as part of `flowcode.spectra`, but those will have been selected on the
#' FlowCode backbone sample. If the thresholds for positivity of any of the
#' fluorophores used to identify the FlowCodes differs in the fully stained
#' samples being unmixed here, provide a new set of thresholds from the
#' ThresholdApp. Run this on an unmixed sample (OLS from the instrument is fine),
#' and refer to the new thresholds file here. Default is `NULL`, which will be
#' ignored.
#' @param weights Optional numeric vector of weights (one per fluorescent
#' detector). Default is `NULL`, in which case weighting will be done by
#' channel means (Poisson variance).
#' @param cell.weighting Logical, whether to use cell-specific weighting for a
#' more Poisson-like unmixing. Default is `FALSE`.
#' @param cell.weight.regularize Logical, whether to regularize cell-specific
#' weights towards the bulk mean weighting set by `weights`. 50:50 averaging.
#' Default is `TRUE`. Only active if `cell.weighting=TRUE`.
#' @param k Numeric, controls the number of variants tested for each fluorophore,
#' autofluorescence and FRET spectrum. Default is `10`. Values up to `10` provide
#' additional benefit in unmixing quality, `1` will be fastest.
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
    weights = NULL,
    cell.weighting = FALSE,
    cell.weight.regularize = TRUE,
    k = 10,
    output.dir = NULL,
    file.suffix = NULL,
    include.imaging = FALSE,
    parallel = TRUE,
    threads = if ( parallel ) 0 else 1,
    verbose = TRUE
  ) {

  # create the output folder if it doesn't exist
  if ( is.null( output.dir ) )
    output.dir <- asp$unmixed.fcs.dir
  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

  # the unmixing method is "AutoSpectral"
  method <- "AutoSpectral"

  # import fcs, without warnings for fcs 3.2
  fcs.data <- suppressWarnings(
    flowCore::read.FCS(
      fcs.file,
      transformation = FALSE,
      truncate_max_range = FALSE,
      emptyValue = FALSE
    )
  )

  fcs.keywords <- flowCore::keyword( fcs.data )
  file.name <- flowCore::keyword( fcs.data, "$FIL" )

  # deal with manufacturer peculiarities in writing fcs files
  if ( asp$cytometer %in% c( "ID7000", "Mosaic" ) ) {
    file.name <- sub(
      "([ _])Raw(\\.fcs$|\\s|$)",
      paste0( "\\1", method, "\\2" ),
      file.name,
      ignore.case = TRUE
    )
  } else if ( grepl( "Discover", asp$cytometer ) ) {
    file.name <- fcs.keywords$FILENAME
    file.name <- sub( ".*\\/", "", file.name )
    file.name <- sub( ".fcs", paste0( " ", method, ".fcs" ), file.name )
  } else {
    file.name <- sub( ".fcs", paste0( " ", method, ".fcs" ), file.name )
  }

  if ( !is.null( file.suffix ) )
    file.name <- sub( ".fcs", paste0( " ", file.suffix, ".fcs" ), file.name )

  # extract exprs
  fcs.exprs <- flowCore::exprs( fcs.data )
  rm( fcs.data )
  original.param <- colnames( fcs.exprs )

  spectral.channel <- colnames( spectra )
  spectral.exprs <- fcs.exprs[ , spectral.channel, drop = FALSE ]

  other.channels <- setdiff( colnames( fcs.exprs ), spectral.channel )

  # remove height and width if present
  suffixes <- c( "-H", "-W" )
  for ( ch in spectral.channel[ grepl( "-A$", spectral.channel ) ] ) {
    base <- sub( "-A$", "", ch )
    other.channels <- setdiff( other.channels, paste0( base, suffixes ) )
  }
  other.exprs <- fcs.exprs[ , other.channels, drop = FALSE ]

  rm( fcs.exprs )

  if ( grepl( "Discover", asp$cytometer ) & !include.imaging )
    other.exprs <- other.exprs[ , asp$time.and.scatter ]

  # define weights
  if ( is.null( weights ) ) {
    # weights are inverse of channel variances (mean if Poisson)
    weights <- abs( colMeans( spectral.exprs ) )
    weights[ weights < 1e-6 ] <- 1e-6
    weights <- 1 / weights
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

  # unmix the data, correcting FRET, autofluorescence and spillover errors where possible
  unmixed.data <- unmix.flowcode(
    raw.data = spectral.exprs,
    spectra = spectra,
    af.spectra = af.spectra,
    spectra.variants = spectra.variants,
    flowcode.spectra = flowcode.spectra,
    asp = asp,
    thresholds = flowcode.thresholds,
    weights = weights,
    cell.weighting = cell.weighting,
    cell.weight.regularize = cell.weight.regularize,
    k = k,
    parallel = parallel,
    threads = threads,
    verbose = verbose
  )

  # combine with non-fluorescence data (scatter, time, etc.)
  unmixed.data <- cbind( other.exprs, unmixed.data )

  # fix any NA values (e.g., plate location with S8)
  if ( anyNA( unmixed.data ) )
    unmixed.data[ is.na( unmixed.data ) ] <- 0

  # update keywords----------
  # identify non-parameter keywords
  non.param.keys <- fcs.keywords[ !grepl( "^\\$?P\\d+", names( fcs.keywords ) ) ]
  if ( asp$cytometer == "Mosaic" )
    non.param.keys <- non.param.keys[ !grepl( "^\\$?CH\\d+", names( non.param.keys ) ) ]

  # build lookup
  pN.keys <- grep( "^\\$?P\\d+N$", names( fcs.keywords ), value = TRUE )
  param.lookup <- lapply( pN.keys, function( k ) {
    p.idx <- sub( "\\$?P(\\d+)N", "\\1", k )
    matches <- grep( paste0( "^\\$?P", p.idx, "(?:[A-Z]+)$" ), names( fcs.keywords ),
                     value = TRUE )
    stats::setNames( fcs.keywords[ matches], matches )
  } )
  # name the list by parameter name
  names( param.lookup ) <- sapply( pN.keys, function( k ) fcs.keywords[[ k ]])

  # keywords for new parameters
  param.keywords <- list()
  n.param <- ncol( unmixed.data )

  # check all parameters and update as needed
  for ( i in seq_len( n.param ) ) {
    p.name <- colnames( unmixed.data )[ i ]

    if ( p.name %in% original.param ) {
      # retain keywords from original file if present
      old.entry <- param.lookup[[ p.name ]]
      if ( !is.null( old.entry ) ) {
        # update index to current parameter number
        names( old.entry ) <- sub( "^\\$P\\d+", paste0( "$P", i ), names( old.entry ) )
        param.keywords <- c( param.keywords, old.entry )

      } else {
        # fallback if missing
        param.keywords[[ paste0( "$P", i, "N" )]] <- p.name
        param.keywords[[ paste0( "$P", i, "S" )]] <- p.name
      }

    } else {
      # keywords for new unmixed parameters
      bit.depth <- if ( !is.null( asp$bit.depth ) ) asp$bit.depth else "32"

      param.keywords[[ paste0( "$P", i, "N" ) ]] <- p.name
      param.keywords[[ paste0( "$P", i, "B" ) ]] <- as.character( bit.depth )
      param.keywords[[ paste0( "$P", i, "E" ) ]] <- "0,0"
      param.keywords[[ paste0( "$P", i, "R" ) ]] <- as.character( asp$expr.data.max )
      param.keywords[[ paste0( "$P", i, "DISPLAY" ) ]] <- "LOG"
      param.keywords[[ paste0( "$P", i, "TYPE" ) ]] <- "Fluorescence"

      # exception for AF.Index
      if ( p.name == "AF.Index" ) {
        param.keywords[[ paste0( "$P", i, "DISPLAY" ) ]] <- "LIN"
        param.keywords[[ paste0( "$P", i, "TYPE" ) ]] <- "AF_Index"
      }

      # assign $PnS (stain) based on flow.control
      f.idx <- match( p.name, flow.control$fluorophore )
      marker <- if ( !is.na( f.idx ) ) flow.control$antigen[ f.idx ] else ""
      param.keywords[[ paste0( "$P", i, "S" ) ]] <- as.character( marker )
    }
  }

  # combine new keywords with original keywords
  new.keywords <- utils::modifyList(
    utils::modifyList( non.param.keys, param.keywords ),
    list(
      "$FIL" = file.name,
      "$PAR" = as.character( n.param ),
      "$UNMIXINGMETHOD" = method,
      "$AUTOSPECTRAL" = as.character( utils::packageVersion( "AutoSpectral" ) )
    )
  )

  # add weights used to a new keyword
  if ( !is.null( weights ) ) {
    # add weights to a new keyword in correct format
    weights.str <- paste(
      c(
        length( spectral.channel ),
        spectral.channel,
        formatC(
          weights,
          digits = 8,
          format = "fg"
        )
      ),
      collapse = ","
    )
    new.keywords[[ "$WEIGHTS" ]] <- weights.str
  }

  # add spectra used to a new keyword
  fluor.n <- nrow( spectra )
  detector.n <- ncol( spectra )
  vals <- as.vector( t( spectra ) )

  # correctly format spillover/spectra
  formatted.vals <- formatC( vals, digits = 8, format = "fg", flag = "#" )
  spill.string <- paste(
    c( fluor.n, detector.n, rownames( spectra ), colnames( spectra ), formatted.vals ),
    collapse = ","
  )
  new.keywords[[ "$SPECTRA" ]] <- spill.string
  new.keywords[[ "$FLUOROCHROMES" ]] <- paste( rownames( spectra ), collapse = "," )

  # add AF spectra as a keyword
  if ( !is.null( af.spectra ) ) {
    af.n <- nrow( af.spectra )
    vals <- as.vector( t( af.n ) )
    formatted.vals <- formatC( vals, digits = 8, format = "fg", flag = "#" )
    af.string <- paste(
      c( af.n, detector.n, rownames( af.n ), colnames( af.spectra ), formatted.vals ),
      collapse = ","
    )
    new.keywords[[ "$AUTOFLUORESCENCE" ]] <- af.string
  }

  ### define new FCS file
  # append "-A" to fluorophore and AF channel names
  fluor.orig <- colnames( unmixed.data )
  colnames( unmixed.data ) <-
    ifelse( fluor.orig %in% c( rownames( spectra ), "AF" ),
            paste0( fluor.orig, "-A" ),
            fluor.orig )

  # create the flowFrame for writing the FCS file
  flow.frame <- suppressWarnings( flowCore::flowFrame( unmixed.data ) )
  param.desc <- flowCore::parameters( flow.frame )@data$desc

  # add marker names to description
  for ( i in seq_len( n.param ) ) {
    orig.name <- fluor.orig[ i ]
    # get the marker from flow.control
    f.idx <- match( orig.name, flow.control$fluorophore )
    if ( !is.na( f.idx ) )
      param.desc[ i ] <- as.character( flow.control$antigen[ f.idx ] )
  }

  # write the parameter to the flowFrame
  flowCore::parameters( flow.frame )@data$desc <- param.desc
  keyword( flow.frame ) <- new.keywords

  # save file ---------
  write.FCS(
    flow.frame,
    filename = file.path( output.dir, file.name )
  )

}



