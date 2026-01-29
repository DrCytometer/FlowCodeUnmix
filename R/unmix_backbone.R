# unmix_backbone.r

#' @title Unmix Backbone
#'
#' @description
#' Unmix the backbone sample, stained only with anti-FlowCode epitope antibodies,
#' for downstream identification of FRET/unmixing errors on a per-combo basis.
#'
#' @importFrom sp point.in.polygon
#'
#' @param flowcode.backbone.fcs File name and path for the FlowCode backbone
#' control FCS file. This should be a sample of cells, ideally the same cell
#' source as your single-stained control samples. These cells should be stained
#' with all the FlowCode epitope tag antibodies present in the fully stained
#' sample and nothing else. All tag combinations should be present and well
#' represented for best results. A minimum of 100 cells per combination is
#' recommended; 2000+ is ideal.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#' @param flow.control A list containing flow cytometry control parameters.
#' @param asp The AutoSpectral parameter list. Prepare using
#' `get.autospectral.param`.
#' @param flowcode.combo.file File name and path to the CSV file containing the
#' information describing your FlowCode library. Describes the valid combinations
#' of FlowCodes. Structure: One row per combination. Columns are `Id`,
#' `Procode.tag1`, `Procode.tag2` and `Procode.tag3`, describing the name (e.g.,
#'  CRISPR target), and three epitopes for the combination, respectively.
#' @param n.cells.gate Numeric, default `100`. The number of cells that will be
#' used to define the gate. A gate region will be defined using `define.gate()`
#' around the brightest `n.cells.gate` for each fluorophore used to identify
#' FlowCode epitope tags. Scatter coordinates of all brightest events will be
#' combined to determine the optimal gating region.
#' @param output.dir Destination folder for output plots and saved spectra.
#' Default is `./flowcode_spectra`.
#' @param filename Name of the output RDS file, default is `FlowCode_Backbone.rds`.
#' @param verbose Logical, whether to send messages to the console.
#' Default is `TRUE`.
#'
#' @return None. Saves RDS objects containing the processed data to disk.
#'
#' @export

unmix.backbone <- function(
    flowcode.backbone.fcs,
    spectra,
    af.spectra,
    flow.control,
    asp,
    flowcode.combo.file,
    n.cells.gate = 100,
    output.dir = "./flowcode_spectra",
    filename = "FlowCode_Backbone.rds",
    verbose = TRUE
  ) {

  if ( !requireNamespace( "AutoSpectral", quietly = TRUE ) ) {
    stop(
      "The AutoSpectral package is required but is not installed or not available.",
      call. = FALSE
    )
  }

  # set up
  if ( nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided as a matrix." )

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  fluorophores <- rownames( spectra )
  fluorophore.n <- length( fluorophores )
  spectral.channel <- colnames( spectra )
  detector.n <- length( spectral.channel )
  af.n <- nrow( af.spectra )

  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

  if ( verbose ) message( "Reading in data" )

  # read in combo file
  combo.df <- utils::read.csv( flowcode.combo.file )
  flowcode.tags <- unique( unlist( combo.df[ , -1 ] ) )

  # define flowcode channel-tag correspondence, case-independent
  # check against flow.control
  tag.lookup <- toupper( flowcode.tags )
  antigen.lookup <- toupper( flow.control$antigen )

  flowcode.fluors <- flow.control$fluorophore[ match( tag.lookup, antigen.lookup ) ]

  names( flowcode.fluors ) <- flowcode.tags

  # read in backbone
  backbone <- suppressWarnings(
    flowCore::exprs(
      flowCore::read.FCS(
        flowcode.backbone.fcs,
        transformation = NULL,
        truncate_max_range = FALSE,
        emptyValue = FALSE
      )
    )
  )
  raw.data <- backbone[ , spectral.channel ]
  scatter.data <- backbone[ , asp$default.scatter.parameter ]

  # initial unmixing without any AF
  if ( verbose ) message( "Initializing unmix for backbone" )

  # set WLS unmixing algorithm
  unmix <- AutoSpectral::unmix.wls.fast

  # define the weights
  weights <- abs( colMeans( raw.data ) )
  weights[ weights < 1e-6 ] <- 1e-6
  weights <- 1 / weights

  # unmix the data
  unmixed <- unmix( raw.data, spectra, weights )

  # before going any further,
  # define gate based on location of FlowCode-expressing cells
  gate <- define.gate(
    scatter.data = scatter.data,
    unmixed.data = unmixed,
    flowcode.fluors = flowcode.fluors,
    asp = asp,
    n.cells = n.cells.gate
  )

  inside <- sp::point.in.polygon(
    point.x = scatter.data[ , 1 ],
    point.y = scatter.data[ , 2 ],
    pol.x = gate$x,
    pol.y = gate$y
  )

  cells.in.gate <- which( inside > 0 )

  # restrict to gated events
  raw.data <- raw.data[ cells.in.gate, ]
  unmixed <- unmixed[ cells.in.gate, ]
  scatter.data <- scatter.data[ cells.in.gate, ]

  # calculate initial error (worst-case: incorrect fluorophore signal)
  error <- rowSums( abs( unmixed[ , fluorophores, drop = FALSE ] ) )

  # set up
  combined.spectra <- matrix(
    NA_real_,
    nrow = fluorophore.n + 1,
    ncol = detector.n
  )
  colnames( combined.spectra ) <- colnames( spectra )
  fluors.af <- c( fluorophores, "AF" )
  rownames( combined.spectra ) <- fluors.af
  combined.spectra[ 1:fluorophore.n, ] <- spectra
  cell.n <- nrow( raw.data )
  initial.af <- matrix(
    0,
    nrow = cell.n,
    ncol = 2
  )
  colnames( initial.af ) <- c( "AF", "AF Index" )
  unmixed <- cbind( unmixed, initial.af )
  fitted.af <- matrix(
    0,
    nrow = cell.n,
    ncol = detector.n
  )

  if ( verbose ) message( "Extracting AF cell-by-cell from backbone" )

  # unmix, removing AF
  for ( af in seq_len( af.n ) ) {

    # set this AF as the spectrum to use
    combined.spectra[ fluorophore.n + 1, ] <- af.spectra[ af, , drop = FALSE ]

    # unmix with this AF
    unmixed.af <- unmix( raw.data, combined.spectra, weights )

    error.af <- rowSums( abs( unmixed.af[ , fluorophores, drop = FALSE ] ) )
    improved <- which( error.af < error )

    # update improved cells
    if ( length( improved ) > 0 ) {
      error[ improved ] <- error.af[ improved ]
      unmixed[ improved, fluors.af ] <- unmixed.af[ improved, ]
      unmixed[ improved, "AF Index" ] <- af

      # update AF fitted values with improved cells
      fitted.af[ improved, ] <- unmixed.af[ improved, "AF", drop = FALSE ] %*%
        af.spectra[ af, , drop = FALSE ]
    }
  }

  # set remaining raw
  remaining.raw <- raw.data - fitted.af

  # construct RDS
  backbone.data <- list(
    Unmixed = cbind( unmixed, scatter.data ),
    Raw = remaining.raw,
    Flowcode.fluors = flowcode.fluors,
    Spectra = spectra,
    Combos = combo.df
  )

  # save all data (large) in an RDS file
  saveRDS(
    backbone.data,
    file.path( output.dir, filename )
  )

  # downsample
  if ( cell.n > 30000 ) {
    set.seed( 42 )
    idx <- sample( seq_len( cell.n ), 30000 )
  } else {
    idx <- seq_len( cell.n )
  }

  # construct RDS
  backbone.data <- list(
    Unmixed = cbind( unmixed[ idx, ], scatter.data[ idx, ] ),
    Flowcode.fluors = flowcode.fluors,
    Spectra = spectra,
    Combos = combo.df
  )

  # save the smaller version of the data in an RDS file
  filename <- paste( "Small", filename, sep = "_" )
  saveRDS(
    backbone.data,
    file.path( output.dir, filename )
  )

}
