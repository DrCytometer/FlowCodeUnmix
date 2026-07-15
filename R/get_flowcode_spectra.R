# get_flowcode_spectra.r

#' @title Get FlowCode Spectra
#'
#' @description
#' Measure spectral unmixing errors per barcode in a sample stained only with
#' FlowCode tag antibodies (`flowcode.backbone.fcs`). Measures FRET error and
#' returns the calculated correction matrices and thresholds for debarcoding.
#'
#' @importFrom FlowSOM SOM
#'
#' @param backbone.rds File name and path to the full (large) rds file produced
#' by running `unmix.backbone`. For this to work well, all combinations present
#' in the experimental samples should be well represented in the
#' backbone sample used.
#' @param thresholds.file File name and path to the thresholds CSV file produced
#' using the threshold setting app.
#' @param asp The AutoSpectral parameter list. Prepare using
#' `get.autospectral.param`.
#' @param output.dir Destination folder for output plots and saved spectra.
#' Default is `./flowcode_spectra`.
#' @param filename Logical, whether to create and save plots of the data
#' correction. Default is `TRUE`.
#' @param figures Logical, controls whether graphs of the process are produced.
#' Default is `TRUE`.
#' @param plot.corrections Logical, default is `FALSE`. If `TRUE`, the full
#' FlowCode unmixing pipeline (`unmix.flowcode()` /
#' `FlowCodeUnmixRcpp::unmix.flowcode.cpp.joint()`) is run on the backbone and
#' before/after plots of the FRET-corrected data (FlowCode channels only) are
#' produced in `output.dir`. Requires both `spectra.variants` and `af.spectra`
#' to be supplied.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress
#' messaging.
#' @param background.pool.n Integer, default `5000`. Maximum number of
#' untransduced cells used as the scatter-matching / AF-PCA reference pool.
#' Downsampled from `untransduced.idx` when more are available.
#' @param spectra.variants Optional. The list returned by
#' `AutoSpectral::get.spectral.variants()` (must contain a `$variants`
#' element). When supplied, each combo fluorophore's per-cell best-fitting
#' spectral variant is substituted before the FRET residual is computed
#' (applied identically to the combo's own cells and to the untransduced
#' background pool), so normal single-colour spectral drift is not mistaken
#' for FRET. Also required, together with `af.spectra`, for
#' `plot.corrections = TRUE`. Default `NULL`, which skips the pre-clean step
#' and disables `plot.corrections`.
#' @param af.spectra Discrete autofluorescence spectral library (fluorophores
#' in rows are not expected here -- AF components in rows, detectors in
#' columns; at least 2 rows), as produced by `AutoSpectral::get.af.spectra()`
#' and already used to build `backbone$Raw` via `unmix.backbone()`. Required
#' only when `plot.corrections = TRUE`: used to reconstruct the
#' pre-AF-subtraction raw signal (so the full pipeline's own AF fitting isn't
#' applied on top of data that's already had AF removed) and to run that
#' pipeline. Default `NULL`.
#' @param parallel Logical, default is `TRUE`, which enables parallel processing
#' for per-cell unmixing methods.
#' @param threads Numeric, defaults to all available cores if `parallel=TRUE`.
#' @param n.cells.per.combo Integer, default `500`. Maximum number of
#' brightest events (ranked by summed signal in a combo's own tag channels,
#' using the existing `backbone$Unmixed` values) retained per combo before
#' FRET spectra are measured.
#' @param som.dim Integer, default `10`. Base SOM grid side length; shrunk
#' per-combo (in a loop-local variable, never mutating this default) for
#' combos with fewer than 500 retained cells.
#' @param fret.median.qc Logical, default `FALSE`. If `TRUE`, events are
#' screened before SOM clustering by cosine similarity to the combo's own
#' median FRET direction (`median.fret.spectrum`, computed from the full,
#' unfiltered event set so the reference isn't circular). This is a loose
#' sanity filter, not a same-shape filter: unlike a single fluorophore's
#' spectrum, FRET within a combo isn't guaranteed to have one true direction
#' (different tag pairings can plausibly produce more than one), so this is
#' only meant to drop events pointing away from the combo's overall trend
#' (see `fret.median.sim.threshold`), leaving SOM free to resolve genuine
#' multi-modal structure among the rest.
#' @param fret.median.sim.threshold Numeric in `[-1, 1]`, default `0`. Cosine
#' similarity cutoff used by `fret.median.qc`; `0` only excludes events
#' pointing into the opposite half-space from the combo's median direction.
#' Raise towards e.g. `0.3`-`0.5` for a stricter filter. Ignored if
#' `fret.median.qc = FALSE`. If fewer than 20 events would survive for a
#' given combo, this QC step is skipped for that combo (with a warning)
#' rather than starving SOM of input.
#' @param fret.background.qc Logical, default `FALSE`. If `TRUE`, candidate
#' FRET variants (SOM centroids, never the median row) are screened after
#' clustering against the same FRET-shaped signal computed from the
#' untransduced/background pool, run through the identical combo-specific
#' `clean.combo.data()` processing. Background cells cannot have real FRET
#' for this combo, so whatever FRET-shaped direction shows up in them after
#' identical processing is leftover contamination (residual AF, non-specific
#' staining, background-subtraction artifacts) rather than signal; variants
#' that closely resemble it (by absolute cosine similarity, since
#' contamination has no consistent sign) are dropped. Unlike
#' `fret.median.qc`, this doesn't assume a single true FRET direction, so it
#' complements rather than duplicates it.
#' @param fret.background.sim.threshold Numeric in `[0, 1]`, default `0.8`.
#' Absolute cosine similarity cutoff used by `fret.background.qc`; a
#' candidate variant is dropped if its similarity to the background pool's
#' own FRET-shaped direction is at or above this value. Ignored if
#' `fret.background.qc = FALSE`.
#'
#' @return A named list with two elements:
#' 1) Thresholds: numeric thresholds for debarcoding
#' 2) FRET: A list of matrices of FRET spectral variations per barcode combo
#' 3) Flowcode.fluors: Character vector of fluorophores linked to the FlowCodes
#' 4) Combos: The valid combinations, specified by the combination file.
#'
#' @export

get.flowcode.spectra <- function(
    backbone.rds,
    thresholds.file,
    asp,
    output.dir = "./flowcode_spectra",
    filename = "FlowCode_Spectra.rds",
    figures = TRUE,
    plot.corrections = TRUE,
    verbose = TRUE,
    background.pool.n = 5000,
    spectra.variants = NULL,
    af.spectra = NULL,
    parallel = TRUE,
    threads = if ( parallel ) 0 else 1,
    n.cells.per.combo = 500,
    som.dim = 10L,
    fret.median.qc = FALSE,
    fret.median.sim.threshold = 0,
    fret.background.qc = FALSE,
    fret.background.sim.threshold = 0.8
) {

  if ( !requireNamespace( "AutoSpectral", quietly = TRUE ) )
    stop( "The AutoSpectral package is required but is not installed or not available.", call. = FALSE )
  if ( !requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) ||
       !( "unmix.autospectral.rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) )
    stop( "This requires the `AutoSpectralRcpp` package (providing `unmix.autospectral.rcpp()`).", call. = FALSE )
  if ( is.null( af.spectra ) || nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided.", call. = FALSE )

  if ( verbose ) message( "Reading thresholds file" )
  flowcode.thresholds.df <- utils::read.csv( thresholds.file )
  flowcode.thresholds <- flowcode.thresholds.df$Threshold_Raw
  names( flowcode.thresholds ) <- flowcode.thresholds.df$Fluor

  if ( verbose ) message( "Loading RDS file" )
  backbone <- readRDS( backbone.rds )

  expected <- list(
    Raw             = list( class = c("matrix"), nrow_min = 100 ),
    Flowcode.fluors = list( class = c("character"), length_min = 4 ),
    Spectra         = list( class = c("matrix"), nrow_min = 2 ),
    Combos          = list( class = c("data.frame"), nrow_min = 10 )
  )
  check.rds.contents( backbone, required = expected )
  if ( verbose ) message( "RDS file has passed checks" )

  spectra <- backbone$Spectra
  combo.df <- backbone$Combos
  flowcode.fluors <- backbone$Flowcode.fluors

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  fluorophores <- rownames( spectra )
  spectral.channel <- colnames( spectra )
  detector.n <- length( spectral.channel )
  combo.n <- nrow( combo.df )
  non.fc.fluors <- setdiff( fluorophores, flowcode.fluors )

  if ( !dir.exists( output.dir ) ) dir.create( output.dir )

  # set number of threads to use
  if ( parallel ) {
    if ( is.null( threads ) || threads == 0 ) threads <- parallelly::availableCores()
  } else {
    threads <- 1
  }
  if ( is.null( noise.floor ) ) noise.floor <- 125

  # combo-logical matrix -- built up front now, since background computation
  # (below) needs it before the per-combo loop, not after
  flowcode.combo.logical <- matrix(
    0L, nrow = combo.n, ncol = length( flowcode.fluors ),
    dimnames = list( combo.df$Id, flowcode.fluors )
  )
  tag.cols <- c( "Procode.tag1", "Procode.tag2", "Procode.tag3" )
  for ( i in seq_len( combo.n ) ) {
    tags <- unlist( combo.df[ i, tag.cols ], use.names = FALSE )
    tags <- flowcode.fluors[ tags ]
    flowcode.combo.logical[ i, tags ] <- 1L
  }

  ### Stage 1: debarcoding-quality unmix -- same algorithm and tuning as
  ### unmix.flowcode.cpp.staged()'s Stage 1, so characterization and
  ### correction agree on debarcoding, background, and AF from the start.
  if ( verbose ) message( paste0( "\033[33m", "Stage 1: AutoSpectral unmix for debarcoding", "\033[0m" ) )

  step1 <- AutoSpectralRcpp::unmix.autospectral.rcpp(
    raw.data               = backbone$Raw,
    spectra                = spectra,
    af.spectra             = af.spectra,
    spectra.variants       = spectra.variants,
    pipeline               = "joint",
    verbose                = FALSE,
    parallel               = parallel,
    threads                = threads
  )

  if ( verbose ) message( paste0( "\033[33m", "Debarcoding the backbone", "\033[0m" ) )
  flowcode.ids <- debarcode(
    unmixed = step1,
    flowcode.fluors = flowcode.fluors,
    flowcode.thresholds = flowcode.thresholds,
    combo.df = combo.df
  )
  valid.idx <- which( flowcode.ids %in% seq_len( combo.n ) )

  if ( figures ) {
    tag.expression.plot(
      unmixed = step1[ valid.idx, ],
      flowcode.fluors = flowcode.fluors,
      flowcode.ids = flowcode.ids[ valid.idx ],
      combo.df = combo.df,
      output.dir = output.dir,
      asp = asp
    )
  }

  # resolved AF spectrum, exported below for reuse at correction time
  Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
  af.idx <- Mode( step1[ valid.idx, "AF Index" ] )
  af.mean.row <- as.numeric( af.spectra[ af.idx, ] )

  ### Shared background computation (same function unmix_flowcode_staged.R uses)
  if ( verbose ) message( paste0( "\033[33m", "Measuring pooled background", "\033[0m" ) )
  bg <- subtract.flowcode.background(
    raw.data               = backbone$Raw,
    step1.unmixed          = step1,
    flowcode.ids           = flowcode.ids,
    fc.fluor.names         = flowcode.fluors,
    flowcode.combo.logical = flowcode.combo.logical,
    spectra                = spectra,
    combo.n                = combo.n
  )
  raw.bg.corrected <- bg$raw.bg.corrected

  untransduced.idx <- which( flowcode.ids == 0 )
  set.seed( asp$bird.seed )
  background.pool.idx <- if ( length( untransduced.idx ) > background.pool.n ) {
    sample( untransduced.idx, background.pool.n )
  } else {
    untransduced.idx
  }

  ### Full-panel + fixed-AF basis, precomputed once
  basis.full <- rbind( spectra, AF = af.mean.row )

  # Unmix everything, clamp what we structurally know is spurious, take the
  # residual. `zero.other = TRUE` is valid ONLY here (backbone is FlowCode-
  # antibody-only, so every non-FlowCode fluor really is 0) -- correction
  # time (fit_flowcode_fret_cpp.cpp) never does this, since other markers
  # are genuinely expressed on real samples and must stay free-fit.
  clamped.residual <- function( raw.bg.sub, off.fluors, zero.other = TRUE ) {
    full.coefs <- AutoSpectral::unmix.ols( raw.bg.sub, basis.full )
    full.coefs[ , off.fluors ] <- 0
    if ( zero.other ) full.coefs[ , non.fc.fluors ] <- 0
    full.coefs[ full.coefs < 0 ] <- 0
    recon <- full.coefs %*% basis.full
    raw.bg.sub - recon
  }

  if ( verbose ) message( paste0( "\033[33m", "Calculating FRET errors per FlowCode combination", "\033[0m" ) )

  fret.spectra <- list()

  for ( id in seq_len( combo.n ) ) {
    combo.idx <- which( flowcode.ids == id )
    n.cells <- length( combo.idx )

    combo.tags <- unlist( combo.df[ id, tag.cols ] )
    combo.fluors <- flowcode.fluors[ match( combo.tags, names( flowcode.fluors ) ) ]
    off.fluors <- setdiff( flowcode.fluors, combo.fluors )

    if ( n.cells > 3 ) {

      # rank by Step 1's own coefficients now, not a separately-computed unmix
      combo.brightness <- step1[ combo.idx, combo.fluors, drop = FALSE ]
      if ( n.cells > n.cells.per.combo ) {
        brightest.events <- order( rowSums( combo.brightness ), decreasing = TRUE )[ 1:n.cells.per.combo ]
        combo.idx <- combo.idx[ brightest.events ]
        n.cells   <- n.cells.per.combo
      }

      input.data <- clamped.residual( raw.bg.corrected[ combo.idx, , drop = FALSE ], off.fluors )

      median.fret <- apply( input.data, 2, stats::median )
      median.fret.spectrum <- median.fret / max( abs( median.fret ) )

      if ( fret.background.qc ) {
        background.input.data <- clamped.residual(
          raw.bg.corrected[ background.pool.idx, , drop = FALSE ], flowcode.fluors  # everything's "off" for the untransduced pool
        )
        background.fret <- apply( background.input.data, 2, stats::median )
        background.fret.denom <- max( abs( background.fret ) )
        if ( background.fret.denom <= 0 ) background.fret.denom <- 1
        background.fret.spectrum <- background.fret / background.fret.denom
      }

      if ( fret.median.qc ) {
        ev.max  <- apply( abs( input.data ), 1, max )
        ev.max[ ev.max <= 0 ] <- 1
        ev.norm <- input.data / ev.max
        ev.cosine <- as.vector( ev.norm %*% median.fret.spectrum ) /
          ( sqrt( rowSums( ev.norm^2 ) ) * sqrt( sum( median.fret.spectrum^2 ) ) )
        median.qc.keep <- which( ev.cosine >= fret.median.sim.threshold )

        if ( length( median.qc.keep ) < 20 ) {
          warning( paste0(
            "Fewer than 20 events passed fret.median.qc for combo ",
            combo.df$Id[ id ], "; skipping this QC step for this combo."
          ) )
        } else {
          input.data <- input.data[ median.qc.keep, , drop = FALSE ]
        }
      }

      n.cells.som <- nrow( input.data )
      som.dim.combo <- if ( n.cells.som < 500L ) {
        max( 2L, floor( sqrt( n.cells.som / 3 ) ) )
      } else {
        som.dim
      }

      set.seed( asp$bird.seed )
      if ( requireNamespace( "EmbedSOM", quietly = TRUE ) ) {
        map <- EmbedSOM::SOM(
          input.data,
          xdim = som.dim.combo,
          ydim = som.dim.combo,
          batch = TRUE,
          parallel = TRUE
        )
      } else {
        map <- FlowSOM::SOM(
          input.data,
          xdim = som.dim.combo,
          ydim = som.dim.combo,
          silent = TRUE
        )
      }

      combo.spectra <- t( apply( map$codes[ , spectral.channel ], 1, function( x ) x / max( x ) ) )
      combo.spectra <- as.matrix( stats::na.omit( combo.spectra ) )
      combo.spectra <- rbind( median.fret.spectrum, combo.spectra )

      if ( fret.background.qc && nrow( combo.spectra ) > 1 ) {
        variant.rows <- 2:nrow( combo.spectra )
        variant.mat  <- combo.spectra[ variant.rows, , drop = FALSE ]
        bg.cosine <- as.vector( variant.mat %*% background.fret.spectrum ) /
          ( sqrt( rowSums( variant.mat^2 ) ) * sqrt( sum( background.fret.spectrum^2 ) ) )
        contaminated <- variant.rows[ abs( bg.cosine ) >= fret.background.sim.threshold ]
        if ( length( contaminated ) > 0 )
          combo.spectra <- combo.spectra[ -contaminated, , drop = FALSE ]
      }

      rownames( combo.spectra ) <- paste0( id, 1:nrow( combo.spectra ) )

      if ( figures ) {
        combo.label <- combo.df$Id[ id ]
        plot.title <- paste( "FRET_errors_for", combo.label, paste( combo.fluors, collapse = "_" ), sep = "_" )
        AutoSpectral::spectral.variant.plot.dens(
          spectra.variants = combo.spectra,
          median.spectrum = combo.spectra[ 1, ],
          title = plot.title,
          save = TRUE,
          plot.dir = output.dir
        )
      }

      fret.spectra[[ id ]] <- combo.spectra

    } else {
      combo.label <- combo.df$Id[ id ]
      warning( paste( "Not enough cells (fewer than 3) were present for", combo.label, ":", combo.fluors ) )
      fret.spectra[[ id ]] <- matrix( 0, nrow = 1, ncol = detector.n )
    }
  }

  names( fret.spectra ) <- combo.df$Id

  fret.delta.list <- lapply( names( fret.spectra ), function( fc ) {
    fret.spectra[[ fc ]] - matrix(
      fret.spectra[[ fc ]][ 1, ], nrow = nrow( fret.spectra[[ fc ]] ), ncol = detector.n, byrow = TRUE
    )
  } )
  names( fret.delta.list ) <- names( fret.spectra )
  fret.delta.norms <- lapply( fret.delta.list, function( d ) sqrt( rowSums( d^2 ) ) )

  if ( figures && plot.corrections ) {
    if ( verbose ) message( paste0( "\033[33m", "Re-unmixing the backbone with the production FlowCode pipeline", "\033[0m" ) )

    flowcode.spectra.preview <- list(
      Thresholds = flowcode.thresholds,
      FRET = fret.spectra,
      Flowcode.fluors = flowcode.fluors,
      Combos = combo.df,
      Delta = fret.delta.list,
      Delta.norms = fret.delta.norms,
      Logical.combo = flowcode.combo.logical,
      AF.Spectrum = af.mean.row
    )

    corrected <- unmix.flowcode.cpp.staged(
      raw.data = backbone$Raw,
      spectra = spectra,
      af.spectra = af.spectra,
      spectra.variants = spectra.variants,
      flowcode.spectra = flowcode.spectra.preview,
      asp = asp,
      verbose = FALSE
    )

    tag.expression.plot(
      unmixed = corrected[ valid.idx, ],
      flowcode.fluors = flowcode.fluors,
      flowcode.ids = flowcode.ids[ valid.idx ],
      combo.df = combo.df,
      output.dir = output.dir,
      asp = asp,
      title = "Corrected"
    )
    flowcode.combo.plot(
      unmixed = step1,
      corrected = corrected,
      flowcode.fluors = flowcode.fluors,
      asp = asp,
      output.dir = output.dir
    )
  }

  if ( verbose ) message( paste0( "\033[33m", "Spectral variation computed!", "\033[0m" ) )

  flowcode.spectra <- list(
    Thresholds = flowcode.thresholds,
    FRET = fret.spectra,
    Flowcode.fluors = flowcode.fluors,
    Combos = combo.df,
    Delta = fret.delta.list,
    Delta.norms = fret.delta.norms,
    Logical.combo = flowcode.combo.logical,
    AF.Spectrum = af.mean.row
  )

  saveRDS( flowcode.spectra, file = file.path( output.dir, filename ) )
  return( flowcode.spectra )
}
