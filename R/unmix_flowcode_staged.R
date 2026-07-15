# unmix_flowcode_cpp_staged.r

#' @title Unmix FlowCode (staged C++ pipeline)
#'
#' @description
#' Unmix FlowCode samples using the staged/sequenced pipeline design: a
#' first generic AutoSpectral unmix for debarcoding, a pooled cross-combo
#' background subtraction, a
#' small joint FRET fit against the background-corrected raw data, and a
#' second generic AutoSpectral unmix for the final reported values. This
#' replaces the single fused C++ pass used by `unmix.flowcode.cpp.joint()`
#' with two calls into the generic, FlowCode-agnostic
#' `AutoSpectralRcpp::unmix.autospectral.rcpp()` machinery, bracketing a
#' thin FlowCode-specific stage that only handles debarcoding, background
#' subtraction, and FRET fitting.
#'
#' @details
#' Pipeline stages:
#' \enumerate{
#'   \item \strong{Initial unmix} - `AutoSpectralRcpp::unmix.autospectral.rcpp()`
#'     (`pipeline = "joint"`), unmodified, generic. Off-target precision
#'     doesn't matter yet; this is only used to get fluor coefficients clean
#'     enough for debarcoding.
#'   \item \strong{Debarcode} - `debarcode()` on stage 1's output, unchanged.
#'   \item \strong{Background subtraction} - one population-wide, per-FlowCode-
#'     fluorophore pooled median background is estimated from off-target
#'     observations across every valid combo, then subtracted uniformly from
#'     every cell using the same full-panel reconstruction vector -- the
#'     correction does not vary by combo identity or debarcode id.
#'   \item \strong{FRET fit} - `fit_flowcode_fret_cpp()`. For valid-ID cells
#'     only, every candidate FRET row for the cell's combo is fit jointly
#'     with the combo's active tag fluorophores and the mean AF spectrum
#'     (all free coefficients, off-target coefficients not part of the
#'     basis at all) against the background-corrected raw data; the
#'     candidate with the lowest L1 residual is selected and its fitted
#'     FRET component subtracted.
#'   \item \strong{Final unmix} - `AutoSpectralRcpp::unmix.autospectral.rcpp()`
#'     again, on the background-and-FRET-corrected raw data (valid-ID
#'     cells) / background-corrected-only raw data (id 0 / unrecognized
#'     cells). This is the reported result.
#' }
#'
#' See `flowcode_fret_af_redesign_summary.md` for the full diagnosis this
#' design addresses (AF/FRET competition from a shared error metric, and
#' non-specific-staining over-correction from a single free scalar having to
#' explain both baseline and true FRET at once).
#'
#' This pipeline requires the compiled `AutoSpectralRcpp` package (providing
#' `unmix.autospectral.rcpp()`); there is currently no pure-R fallback for
#' this staged design.
#'
#' @importFrom parallelly availableCores
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in `spectra`.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Row
#' 1 (the mean AF component, as produced by `AutoSpectral::get.af.spectra()`)
#' is used as the fixed-shape AF basis vector in the FRET fit (stage 4); its
#' coefficient there is free, not fixed.
#' @param spectra.variants Named list (names are fluorophores) carrying
#' matrices of spectral signature variations for each fluorophore. Prepare
#' using `get.spectral.variants()`. Ignored (AF-only unmixing at stages 1
#' and 5) if `optimize = FALSE`.
#' @param flowcode.spectra Structured output from `get.flowcode.spectra()`.
#' @param asp The AutoSpectral parameter list. Prepare using
#' `get.autospectral.param()`. Not used directly by this pipeline yet -
#' retained for interface parity with `unmix.flowcode.cpp.joint()`.
#' @param thresholds Optional named numeric vector of positivity thresholds
#' for the FlowCode fluorophores, overriding `flowcode.spectra$Thresholds`.
#' Default `NULL`.
#' @param parallel Logical, whether to use parallel processing. Default `TRUE`.
#' @param threads Numeric, number of threads. Default `0` (all cores) if
#' `parallel = TRUE`, else `1`.
#' @param verbose Logical, default `TRUE`.
#' @param optimize Logical, default `TRUE`. If `FALSE`, `spectra.variants` is
#' not passed through to either AutoSpectral call (AF-only unmixing at
#' stages 1 and 5); debarcoding, background subtraction, and FRET fitting
#' still run regardless.
#' @param cell.weight Logical, default `FALSE`. Per-cell IRLS-style detector
#' weighting, passed through to both `unmix.autospectral.rcpp()` calls. Not
#' applied to the FRET fit itself (stage 4), which is always unweighted in
#' this design.
#' @param noise.floor Numeric, default `NULL` (falls back to `125`). Passed
#' through to `unmix.autospectral.rcpp()`.
#' @param alpha Numeric in \[0, 1\], default `0.5`. Passed through to
#' `unmix.autospectral.rcpp()`.
#' @param collinear.thresh Numeric in \[0, 1\], default `0.5`. Passed through
#' to `unmix.autospectral.rcpp()` as `collinear.threshold`.
#' @param joint.pair.resolution Logical, default `TRUE`. Passed through to
#' `unmix.autospectral.rcpp()`.
#' @param n.passes Numeric, default `1`. Passed through to
#' `unmix.autospectral.rcpp()`.
#' @param n.af.passes Numeric, default `1`. Passed through to
#' `unmix.autospectral.rcpp()`.
#' @param refine.af.quantile Numeric in \[0, 1\], default `0.5`. Passed
#' through to `unmix.autospectral.rcpp()`.
#' @param fret.alpha Numeric in \[0, 1\], default `0.5`. Passed through to
#' `fit_flowcode_fret_cpp()` if used (currently not in use).
#' @param fret.median.only Logical, default `FALSE`. If `TRUE`, skips the
#' per-cell FRET variant scan and always fits/subtracts the first (median)
#' row of the combo's FRET spectra matrix.
#' @param return.diagnostics Logical, default `FALSE`. If `TRUE`, returns a
#' list (`$unmixed`, `$diagnostics`) instead of just the unmixed matrix.
#' `$diagnostics` contains `median.background` (named numeric vector, one
#' entry per FlowCode fluorophore - see stage 3), `flowcode.ids`, `fret.k`,
#' and `fret.index` (the latter two `0` for cells with no FRET fit
#' attempted).
#'
#' @return Unmixed data with cells in rows and fluorophores in columns (plus
#' `AF`, `AF Index`, `FlowCode`, and one column per valid combination), or a
#' list including diagnostics if `return.diagnostics = TRUE`.
#'
#' @export

unmix.flowcode.cpp.staged <- function(
    raw.data,
    spectra,
    af.spectra,
    spectra.variants,
    flowcode.spectra,
    asp,
    thresholds            = NULL,
    parallel              = TRUE,
    threads               = if ( parallel ) 0 else 1,
    verbose               = TRUE,
    optimize              = TRUE,
    cell.weight           = FALSE,
    noise.floor           = NULL,
    alpha                 = 0.5,
    collinear.thresh      = 0.5,
    joint.pair.resolution = TRUE,
    n.passes              = 1,
    n.af.passes           = 1,
    refine.af.quantile    = 0.5,
    fret.alpha            = 0.5,
    fret.median.only      = FALSE,
    return.diagnostics    = FALSE
) {

  ### Setup -----------

  if ( !requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) ||
       !( "unmix.autospectral.rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) )
    stop(
      "This staged FlowCode pipeline requires the `AutoSpectralRcpp` package ",
      "(providing `unmix.autospectral.rcpp()`); there is currently no pure-R ",
      "fallback for this pipeline design.",
      call. = FALSE
    )

  # check for AF in spectra, remove if present
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  fluorophores <- rownames( spectra )

  if ( is.null( af.spectra ) || nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided.", call. = FALSE )

  # unpack FlowCode FRET, thresholds, and combo definitions
  combo.df <- flowcode.spectra$Combos
  combo.n  <- nrow( combo.df )

  flowcode.thresholds <- if ( !is.null( thresholds ) ) thresholds else flowcode.spectra$Thresholds
  flowcode.fluors     <- flowcode.spectra$Flowcode.fluors

  if ( is.null( flowcode.thresholds ) )
    stop( "Check that FlowCode thresholds have been supplied correctly (use get.flowcode.spectra).", call. = FALSE )
  if ( !all( names( flowcode.thresholds ) %in% fluorophores ) )
    stop( "FlowCode thresholds don't match the fluorophore names supplied.", call. = FALSE )

  # Defensively reorder FRET spectra and the combo-logical matrix by
  # combo.df$Id -- flowcode.spectra$FRET is a named list, but downstream we
  # need list position `id` (from debarcode()) to match combo.df's row
  # order exactly, so we fix that correspondence explicitly here rather than
  # relying on construction order.
  fret.spectra <- flowcode.spectra$FRET[ as.character( combo.df$Id ) ]
  flowcode.combo.logical <- flowcode.spectra$Logical.combo[ as.character( combo.df$Id ), , drop = FALSE ]
  fc.fluor.names <- colnames( flowcode.combo.logical )

  bad.fret <- which( vapply(
    fret.spectra,
    function( fm ) !is.matrix( fm ) || ncol( fm ) != ncol( spectra ) || nrow( fm ) < 1,
    logical( 1 )
  ) )
  if ( length( bad.fret ) > 0 )
    stop(
      "`flowcode.spectra$FRET` has malformed entries for combo(s): ",
      paste( names( fret.spectra )[ bad.fret ], collapse = ", " ),
      ". Each element must be a matrix with ", ncol( spectra ), " columns.",
      call. = FALSE
    )

  # set number of threads to use
  if ( parallel ) {
    if ( is.null( threads ) || threads == 0 ) threads <- parallelly::availableCores()
  } else {
    threads <- 1
  }

  if ( is.null( noise.floor ) ) noise.floor <- 125

  sv.pass <- if ( isTRUE( optimize ) ) spectra.variants else NULL

  # shared arguments for both AutoSpectral calls (stages 1 and 5); only
  # `raw.data` differs between them
  common.args <- list(
    spectra               = spectra,
    af.spectra             = af.spectra,
    spectra.variants       = sv.pass,
    pipeline               = "joint",
    verbose                = FALSE,
    parallel               = parallel,
    threads                = threads,
    n.af.passes            = n.af.passes,
    n.passes               = n.passes,
    cell.weight            = cell.weight,
    noise.floor            = noise.floor,
    alpha                  = alpha,
    collinear.threshold    = collinear.thresh,
    joint.pair.resolution  = joint.pair.resolution,
    refine.af.quantile     = refine.af.quantile
  )


  ### Stage 1: initial generic AutoSpectral unmix (for debarcoding) -----------

  if ( verbose ) message( "Stage 1/5: initial AutoSpectral unmix..." )

  step1 <- do.call(
    AutoSpectralRcpp::unmix.autospectral.rcpp,
    c( list( raw.data = raw.data ), common.args )
  )


  ### Stage 2: debarcode -----------

  if ( verbose ) message( "Stage 2/5: debarcoding..." )

  flowcode.ids <- debarcode(
    unmixed             = step1,
    flowcode.fluors     = flowcode.fluors,
    flowcode.thresholds = flowcode.thresholds,
    combo.df            = combo.df
  )


  ### Stage 3: pooled cross-combo background subtraction -----------

  if ( verbose ) message( "Stage 3/5: pooled background subtraction..." )

  bg <- .subtract.flowcode.background(
    raw.data               = raw.data,
    step1.unmixed          = step1,
    flowcode.ids           = flowcode.ids,
    fc.fluor.names         = fc.fluor.names,
    flowcode.combo.logical = flowcode.combo.logical,
    spectra                = spectra,
    combo.n                = combo.n
  )

  raw.bg.corrected  <- bg$raw.bg.corrected
  median.background <- bg$median.background


  ### Stage 4: FRET fit (valid-ID cells only) -----------

  if ( verbose ) message( "Stage 4/5: fitting FRET..." )

  active.fc.idx.list <- lapply( seq_len( combo.n ), function( id ) {
    active.fluors <- fc.fluor.names[ flowcode.combo.logical[ id, ] == 1 ]
    idx <- match( active.fluors, fluorophores )
    if ( anyNA( idx ) )
      stop(
        "FlowCode fluorophore(s) not found in `spectra` for combo '",
        combo.df$Id[ id ], "': ",
        paste( active.fluors[ is.na( idx ) ], collapse = ", " ),
        call. = FALSE
      )
    idx - 1L   # 0-based for C++
  } )

  # determine most used AF spectrum
  # 0-based indices of every FlowCode fluor within `spectra` -- used by the
  # per-cell FRET fit to know which coefficients are eligible to be clamped
  # (off-combo) vs. must stay free-fit ("other", non-FlowCode markers).
  fc.all.idx <- match( fc.fluor.names, fluorophores ) - 1L
  if ( anyNA( fc.all.idx ) )
    stop( "FlowCode fluorophore(s) in `fc.fluor.names` not found in `spectra`.", call. = FALSE )

  # FRET-designated AF spectrum: reuse whatever was resolved at
  # characterization time (get.flowcode.spectra()) for continuity between
  # characterization and correction. Fall back to the empirical mode of this
  # sample's own Step 1 AF Index only if it wasn't exported (older
  # flowcode.spectra objects).
  if ( !is.null( flowcode.spectra$AF.Spectrum ) ) {
    af.mean.row <- as.numeric( flowcode.spectra$AF.Spectrum )
  } else {
    Mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    valid.idx <- which( flowcode.ids %in% seq_len( combo.n ) )
    af.idx <- Mode( step1[ valid.idx, "AF Index" ] )
    af.mean.row <- as.numeric( af.spectra[ af.idx, ] )
  }
  if ( length( af.mean.row ) != ncol( spectra ) )
    stop( "`flowcode.spectra$AF.Spectrum` length must equal ncol(spectra).", call. = FALSE )



  fret.result <- FlowCodeUnmixRcpp::fit_flowcode_fret_lowd_cpp(
    raw_bg_corrected_in = as.matrix( raw.bg.corrected ),
    spectra             = as.matrix( spectra ),
    af_mean_row         = af.mean.row,
    flowcode_ids        = as.integer( flowcode.ids ),
    active_fc_idx_list  = active.fc.idx.list,
    #fc_all_idx          = as.integer( fc.all.idx ),
    fret_spectra_list   = fret.spectra,
    #fret_alpha          = as.numeric( fret.alpha ),
    n_threads           = as.integer( threads ),
    fret_median_only    = isTRUE( fret.median.only )
  )

  raw.corrected <- fret.result$raw.corrected
  colnames( raw.corrected ) <- colnames( raw.bg.corrected )


  ### Stage 5: final generic AutoSpectral unmix -----------

  if ( verbose ) message( "Stage 5/5: final AutoSpectral unmix..." )

  final.unmixed <- do.call(
    AutoSpectralRcpp::unmix.autospectral.rcpp,
    c( list( raw.data = raw.corrected ), common.args )
  )


  ### Build "FlowCode" summary + per-combo channels -----------

  cell.n <- nrow( final.unmixed )
  FlowCode <- rep( 0, cell.n )
  fc.data <- matrix(
    0, nrow = cell.n, ncol = combo.n,
    dimnames = list( NULL, combo.df$Id )
  )

  for ( id in seq_len( combo.n ) ) {
    idx <- which( flowcode.ids == id )
    if ( length( idx ) == 0 ) next

    fc.fluors.id <- fluorophores[ active.fc.idx.list[[ id ]] + 1L ]
    vals <- final.unmixed[ idx, fc.fluors.id, drop = FALSE ]
    med  <- if ( length( fc.fluors.id ) > 1 ) apply( vals, 1, stats::median ) else vals[ , 1 ]

    FlowCode[ idx ]    <- med
    fc.data[ idx, id ] <- med
  }

  final.matrix <- cbind( final.unmixed, FlowCode = FlowCode, fc.data )
  colnames( final.matrix ) <- c( fluorophores, "AF", "AF Index", "FlowCode", combo.df$Id )

  if ( verbose ) message( "Unmixing complete." )

  if ( !return.diagnostics ) return( final.matrix )

  list(
    unmixed = final.matrix,
    diagnostics = list(
      median.background = median.background,
      flowcode.ids       = flowcode.ids,
      fret.k             = as.numeric( fret.result$fret.k ),
      fret.index         = as.integer( fret.result$fret.index ),
      resid.ratio        = as.numeric( fret.result$resid.ratio ),
      leakage.ratio      = as.numeric( fret.result$leakage.ratio )
    )
  )
}


### Internal helpers ---------------------------------------------------------

# Pooled, cross-combo median background per FlowCode fluorophore. For each
# FlowCode fluorophore `f`, pools stage-1 unmixed values in channel `f` across
# every valid-ID cell whose combo does not have `f` active, then takes one
# median. Pooling across combos (rather than per-combo) is what separates
# true non-specific background from FRET: any one combo's systematic FRET in
# `f` is diluted by every other combo where `f` is not elevated.
.pool.flowcode.background <- function(
    step1.unmixed,
    flowcode.ids,
    fc.fluor.names,
    flowcode.combo.logical,
    combo.n
) {

  valid.idx <- which( flowcode.ids %in% seq_len( combo.n ) )

  median.background <- stats::setNames( numeric( length( fc.fluor.names ) ), fc.fluor.names )

  for ( f in fc.fluor.names ) {
    active.combos.f <- which( flowcode.combo.logical[ , f ] == 1 )
    pool.idx <- valid.idx[ !( flowcode.ids[ valid.idx ] %in% active.combos.f ) ]

    if ( length( pool.idx ) == 0 ) {
      warning(
        "No off-target observations available to estimate background for '",
        f, "'; using 0.", call. = FALSE
      )
      next
    }

    median.background[ f ] <- stats::median( step1.unmixed[ pool.idx, f ] )
  }

  median.background
}

# Subtracts the pooled background reconstruction (in raw/spectral space)
.subtract.flowcode.background <- function(
    raw.data,
    step1.unmixed,
    flowcode.ids,
    fc.fluor.names,
    flowcode.combo.logical,
    spectra,
    combo.n
) {

  median.background <- .pool.flowcode.background(
    step1.unmixed, flowcode.ids, fc.fluor.names, flowcode.combo.logical, combo.n
  )

  # reconstruct raw background
  bg.recon <- as.numeric( median.background[ fc.fluor.names ] %*% spectra[ fc.fluor.names, , drop = FALSE ] )

  raw.bg.corrected <- if ( all( bg.recon == 0 ) ) {
    raw.data
  } else {
    sweep( raw.data, 2, bg.recon, "-" )
  }

  list( raw.bg.corrected = raw.bg.corrected, median.background = median.background )
}
