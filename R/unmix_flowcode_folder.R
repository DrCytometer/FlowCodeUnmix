# unmix_flowcode_folder.r

#' @title Unmix All FCS Files in a Directory
#'
#' @description
#' This function unmixes all FCS files in a specified directory using the
#' provided spectra and method, and saves the unmixed FCS files to an output
#' directory of the user's choice.
#'
#' @param fcs.dir Directory containing FCS files to be unmixed.
#' @param spectra Matrix containing spectra information.
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
#' @param weights Optional numeric vector of weights: one per fluorescent
#' detector. Default is `NULL`, in which case weighting will be done by
#' channel means.
#' @param k Numeric, controls the number of variants tested for each fluorophore,
#' autofluorescence and FRET spectrum. Default is `10`. Values up to `10` provide
#' additional benefit in unmixing quality, `1` will be fastest.
#' @param output.dir Directory to save the unmixed FCS files
#' (default is `asp$unmixed.fcs.dir`).
#' @param file.suffix A character string to append to the output file name.
#' Default is `NULL`
#' @param include.imaging Logical indicating whether to include imaging data in
#' the written FCS file: relevant for S8 and A8. Default is `FALSE`
#' @param parallel Logical, default is `TRUE`, which enables parallel processing
#' for per-cell unmixing methods.
#' @param threads Numeric, defaults to all available cores if `parallel=TRUE`.
#' @param verbose Logical, controls messaging. Default is `TRUE`.
#'
#' @return None. Saves the unmixed FCS files to the specified output directory.
#'
#' @export

unmix.flowcode.folder <- function(
    fcs.dir,
    spectra,
    asp,
    flow.control,
    af.spectra,
    spectra.variants,
    flowcode.spectra,
    thresholds.file = NULL,
    weights = NULL,
    k = 10,
    output.dir = NULL,
    file.suffix = NULL,
    include.imaging = FALSE,
    parallel = TRUE,
    threads = if ( parallel ) 0 else 1,
    verbose = TRUE
) {

  # set up, create output folders where FCS files will go
  if ( is.null( output.dir ) )
    output.dir <- asp$unmixed.fcs.dir
  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

  if ( parallel & is.null( threads ) )
    threads <- asp$worker.process.n
  else if ( parallel & threads == 0 )
    threads <- parallelly::availableCores()

  # list all the FCS files
  files.to.unmix <- list.files( fcs.dir, pattern = ".fcs", full.names = TRUE )

  # unmix all files in list
  for ( fc in files.to.unmix ) {
    if ( verbose ) message( paste0( "Unmixing ", fc ) )

    unmix.flowcode.fcs(
      fcs.file = fc,
      spectra = spectra,
      asp = asp,
      flow.control = flow.control,
      af.spectra = af.spectra,
      spectra.variants,
      flowcode.spectra,
      thresholds.file = thresholds.file,
      weights = weights,
      k = k,
      output.dir = output.dir,
      file.suffix = file.suffix,
      include.imaging = include.imaging,
      parallel = parallel,
      threads = threads,
      verbose = verbose
    )
  }

  if ( verbose ) message( "Unmixing complete!" )
}
