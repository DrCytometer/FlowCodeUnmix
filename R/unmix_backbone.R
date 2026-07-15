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
#' @param spectra.variants Optional. The list returned by
#' `AutoSpectral::get.spectral.variants()` (must contain a `$variants`
#' element). When supplied, each combo fluorophore's per-cell best-fitting
#' spectral variant is substituted before the FRET residual is computed
#' (applied identically to the combo's own cells and to the untransduced
#' background pool), so normal single-colour spectral drift is not mistaken
#' for FRET. Also required, together with `af.spectra`, for
#' `plot.corrections = TRUE`. Default `NULL`, which skips the pre-clean step
#' and disables `plot.corrections`.
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
#' @param singlet.quantiles Numeric, default `c( 0.85, 0.975 )`. Quantile
#' cutoffs (FSC ratio, then SSC ratio) for the two-pass Area/Height
#' scatter-ratio doublet gate, mirroring the same gate in
#' `get_spectra_automated.R` (`flowstate::select_singlets`-style). The SSC
#' cutoff is computed only from events that already passed the FSC cutoff
#' (sequential, not independent, gating). Skipped with a message if Height
#' channels for `asp$default.scatter.parameter` aren't present in the file.
#' @param output.dir Destination folder for output plots and saved spectra.
#' Default is `./flowcode_spectra`.
#' @param filename Name of the output RDS file, default is `FlowCode_Backbone.rds`.
#' @param verbose Logical, whether to send messages to the console.
#' Default is `TRUE`.
#' @param parallel Logical, whether to use parallel processing for the per-cell
#' unmixing. Default is `FALSE`.
#' @param threads Numeric. Number of threads to use for parallel processing.
#' Defaults to `1` for sequential processing, or `0` (all cores) if `parallel=TRUE`.
#'
#' @return None. Saves RDS objects containing the processed data to disk.
#'
#' @export

unmix.backbone <- function(
    flowcode.backbone.fcs,
    spectra,
    af.spectra,
    spectra.variants = NULL,
    flow.control,
    asp,
    flowcode.combo.file,
    n.cells.gate = 100,
    singlet.quantiles = c( 0.85, 0.975 ),
    output.dir = "./flowcode_spectra",
    filename = "FlowCode_Backbone.rds",
    verbose = TRUE,
    parallel = TRUE,
    threads = if (parallel) 0 else 1
  ) {

  if ( !requireNamespace( "AutoSpectral", quietly = TRUE ) ) {
    stop(
      "The AutoSpectral package is required but is not installed or not available.",
      call. = FALSE
    )
  }

  # set number of threads to use
  if ( parallel ) {
    if ( is.null( threads ) ) threads <- asp$worker.process.n
    if ( threads == 0 ) threads <- parallelly::availableCores()
  } else {
    threads <- 1
  }
  # set up checks
  if ( nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided as a matrix." )

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  fluorophores <- rownames( spectra )
  fluorophore.n <- length( fluorophores )
  spectral.channel <- colnames( spectra )
  detector.n <- length( spectral.channel )
  af.n <- nrow( af.spectra )

  if ( !dir.exists( output.dir ) ) dir.create( output.dir )

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
  backbone <- AutoSpectral::readFCS( flowcode.backbone.fcs )
  raw.data <- backbone[ , spectral.channel ]
  scatter.data <- backbone[ , asp$default.scatter.parameter ]

  # exclude saturated events
  not.saturated <- rowSums( raw.data >= asp$expr.data.max ) == 0
  raw.data      <- raw.data[ not.saturated, , drop = FALSE ]
  scatter.data  <- scatter.data[ not.saturated, , drop = FALSE ]

  # exclude doublets/multiplets via a two-pass Area/Height scatter-ratio gate
  # (mirrors flowstate::select_singlets, as used in get_spectra_automated.R).
  # A doublet carrying two different barcodes will still debarcode to a
  # single Id downstream; without this, its unmixing residual -- driven by
  # the second, uncounted barcode's own spillover -- would look exactly like
  # FRET and get folded into the measured FRET spectra.
  fsc.a <- asp$default.scatter.parameter[ 1L ]
  ssc.a <- asp$default.scatter.parameter[ 2L ]
  fsc.h <- sub( "-A$", "-H", fsc.a )
  ssc.h <- sub( "-A$", "-H", ssc.a )

  if ( all( c( fsc.h, ssc.h ) %in% colnames( backbone ) ) ) {
    height.data <- backbone[ not.saturated, c( fsc.h, ssc.h ), drop = FALSE ]

    fsc.ratio   <- scatter.data[ , fsc.a ] / ( height.data[ , fsc.h ] + 1e-9 )
    singlet.idx <- fsc.ratio < stats::quantile( fsc.ratio, probs = singlet.quantiles[ 1L ] )

    # SSC cutoff computed only from events already passing the FSC cutoff --
    # sequential, not independent, gating (matches get_spectra_automated.R)
    ssc.ratio <- scatter.data[ singlet.idx, ssc.a ] / ( height.data[ singlet.idx, ssc.h ] + 1e-9 )
    singlet.idx[ singlet.idx ] <- ssc.ratio < stats::quantile( ssc.ratio, probs = singlet.quantiles[ 2L ] )

    raw.data     <- raw.data[ singlet.idx, , drop = FALSE ]
    scatter.data <- scatter.data[ singlet.idx, , drop = FALSE ]
  } else if ( verbose ) {
    message(
      paste0(
        "\033[33m",
        "Height channels for ", fsc.a, "/", ssc.a, " not found: skipping singlet gating.",
        "\033[0m"
      )
    )
  }

  # initial unmixing without any AF
  if ( verbose ) message( "Initializing unmix for backbone" )

  # unmix the data
  unmixed <- AutoSpectral::unmix.ols( raw.data, spectra )

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

  # unmix with AF extraction
  if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) ) {
    unmixed <- AutoSpectralRcpp::unmix.autospectral.rcpp(
      raw.data = raw.data,
      spectra = spectra,
      af.spectra = af.spectra,
      spectra.variants = spectra.variants,
      parallel = parallel,
      threads = threads,
      verbose = FALSE
    )
  } else {
    unmixed <- AutoSpectral::unmix.autospectral(
      raw.data = raw.data,
      spectra = spectra,
      af.spectra = af.spectra,
      spectra.variants = spectra.variants,
      asp = asp,
      parallel = parallel,
      threads = threads,
      verbose = FALSE
    )
  }

  if ( verbose ) message( "Saving data" )

  # construct RDS
  backbone.data <- list(
    Unmixed = cbind( unmixed, scatter.data ),
    Raw = raw.data,
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
  if ( length( cells.in.gate ) > 30000 ) {
    set.seed( asp$bird.seed )
    idx <- sample( seq_len( length( cells.in.gate ) ), 30000 )
  } else {
    idx <- seq_len( length( cells.in.gate ) )
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
