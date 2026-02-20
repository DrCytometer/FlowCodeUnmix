# unmix_flowcode_fcs.r

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
#' @param k Numeric, controls the number of variants tested for each fluorophore,
#' autofluorescence and FRET spectrum. Default is `10`. Values up to `10` provide
#' additional benefit in unmixing quality, `1` will be fastest. Only used if
#' `optimize=TRUE`.
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
    k = 10,
    output.dir = NULL,
    file.suffix = NULL,
    include.imaging = FALSE,
    parallel = TRUE,
    threads = if ( parallel ) 0 else 1,
    verbose = TRUE,
    optimize = TRUE
  ) {

  # check for AutoSpectral in NAMESPACE
  if ( !requireNamespace( "AutoSpectral", quietly = TRUE ) ) stop( "AutoSpectral required." )

  # create the output folder if it doesn't exist
  output.dir <- if ( is.null( output.dir ) ) asp$unmixed.fcs.dir else output.dir
  if ( !dir.exists( output.dir ) ) dir.create( output.dir )

  #  import FCS
  if ( verbose ) message( "Reading FCS: ", fcs.file )
  import <- readFCS( fcs.file, keywords = TRUE )
  fcs.exprs <- import$data
  fcs.keywords <- import$keywords
  original.param <- colnames( fcs.exprs )

  # determine base filename
  method <- "AutoSpectral"
  file.name <- if ( !is.null( fcs.keywords$`$FIL` ) ) fcs.keywords$`$FIL` else basename( fcs.file )

  # handle manufacturer/cytometer specifics
  if ( asp$cytometer %in% c( "ID7000", "Mosaic" ) ) {
    # replace "Raw" with "AutoSpectral" (e.g., "Experiment_Raw.fcs" -> "Experiment_AutoSpectral.fcs")
    file.name <- sub("([ _])Raw(\\.fcs$|\\s|$)", paste0("\\1", method, "\\2"),
                     file.name, ignore.case = TRUE)

  } else if (grepl("Discover", asp$cytometer)) {
    # Discover uses a full path in the FILENAME keyword; strip to basename
    file.name <- if ( !is.null( fcs.keywords$FILENAME ) ) fcs.keywords$FILENAME else file.name
    file.name <- basename( file.name )
    file.name <- sub("\\.fcs$", paste0(" ", method, ".fcs"), file.name, ignore.case = TRUE)

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
  spectral.exprs <- fcs.exprs[ , spectral.channel, drop = FALSE ]
  other.channels <- setdiff( colnames( fcs.exprs ), spectral.channel )

  # remove height and width if present
  for ( ch in spectral.channel[ grepl( "-A$", spectral.channel ) ] ) {
    base <- sub( "-A$", "", ch )
    other.channels <- setdiff( other.channels, paste0( base, c( "-H", "-W" ) ) )
  }
  if (grepl("Discover", asp$cytometer) && !include.imaging) {
    other.channels <- intersect(other.channels, asp$time.and.scatter)
  }
  other.exprs <- fcs.exprs[ , other.channels, drop = FALSE ]
  rm( fcs.exprs ); gc()

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

  # check whether FlowCodeUnmixRcpp in installed
  if ( requireNamespace("FlowCodeUnmixRcpp", quietly = TRUE ) &&
       "unmix.flowcode.cpp" %in% ls( getNamespace( "FlowCodeUnmixRcpp" ) ) ) {
    if ( verbose ) message( "Unmixing using FlowCodeUnmixRcpp" )
    tryCatch(
      unmixed.data <- FlowCodeUnmixRcpp::unmix.flowcode.cpp(
        raw.data = spectral.exprs,
        spectra = spectra,
        af.spectra = af.spectra,
        spectra.variants = spectra.variants,
        flowcode.spectra = flowcode.spectra,
        asp = asp,
        thresholds = flowcode.thresholds,
        k = k,
        parallel = parallel,
        threads = threads,
        verbose = verbose,
        optimize = optimize
      ),
      error = function( e ) {
        warning(
          "FlowCodeUnmixRcpp unmixing failed, falling back to standard FlowCodeUnmix: ",
          e$message,
          call. = FALSE
        )
        unmixed.data <- unmix.flowcode(
          raw.data = spectral.exprs,
          spectra = spectra,
          af.spectra = af.spectra,
          spectra.variants = spectra.variants,
          flowcode.spectra = flowcode.spectra,
          asp = asp,
          thresholds = flowcode.thresholds,
          k = k,
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
    unmixed.data <- unmix.flowcode(
      raw.data = spectral.exprs,
      spectra = spectra,
      af.spectra = af.spectra,
      spectra.variants = spectra.variants,
      flowcode.spectra = flowcode.spectra,
      asp = asp,
      thresholds = flowcode.thresholds,
      k = k,
      parallel = parallel,
      threads = threads,
      verbose = verbose,
      optimize = optimize
    )
  }
  rm( spectral.exprs )

  # combine with non-fluorescence data (scatter, time, etc.)
  final.colnames <- c( colnames( other.exprs ), colnames( unmixed.data ) )
  final.matrix <- matrix( 0, nrow = nrow( unmixed.data ), ncol = length( final.colnames ) )
  colnames( final.matrix ) <- final.colnames
  final.matrix[ , 1:ncol( other.exprs ) ] <- other.exprs
  final.matrix[ , ( ncol( final.matrix ) - ncol( unmixed.data ) + 1 ):ncol( final.matrix ) ] <- unmixed.data
  rm( unmixed.data, other.exprs )

  # fix any NA values (e.g., plate location with S8)
  if ( anyNA( final.matrix ) )
    final.matrix[ is.na( final.matrix ) ] <- 0

  # append "-A" to fluorophore and AF channel names
  colnames( final.matrix ) <-
    ifelse( colnames( final.matrix ) %in% c( rownames( spectra ), "AF" ),
            paste0( colnames( final.matrix ), "-A" ),
            colnames( final.matrix ) )

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
    file.name = file.name
  )

  # save file ---------
  if ( verbose ) message( paste( "Writing:", file.name ) )
  writeFCS( final.matrix, new.keywords, file.name, output.dir )
}



