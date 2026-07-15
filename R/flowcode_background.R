# flowcode_background.r
#
# Shared by unmix_flowcode_staged.R and get_flowcode_spectra.R so the
# background-subtraction methodology used at correction time and at
# characterization time can never silently drift apart.

# Pooled, cross-combo median background per FlowCode fluorophore. For each
# FlowCode fluorophore `f`, pools stage-1 unmixed values in channel `f` across
# every valid-ID cell whose combo does not have `f` active, then takes one
# median. Pooling across combos (rather than per-combo) is what separates
# true non-specific background from FRET: any one combo's systematic FRET in
# `f` is diluted by every other combo where `f` is not elevated.
pool.flowcode.background <- function(
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

# Subtracts the pooled background reconstruction (in raw/spectral space),
# uniformly, from every cell -- including this combo's own active tags, and
# including id-0/unrecognized cells (non-specific staining is a property of
# the fluorophore/reagent, not of transduction status).
subtract.flowcode.background <- function(
    raw.data,
    step1.unmixed,
    flowcode.ids,
    fc.fluor.names,
    flowcode.combo.logical,
    spectra,
    combo.n
) {

  median.background <- pool.flowcode.background(
    step1.unmixed, flowcode.ids, fc.fluor.names, flowcode.combo.logical, combo.n
  )

  bg.recon <- as.numeric( median.background[ fc.fluor.names ] %*% spectra[ fc.fluor.names, , drop = FALSE ] )

  raw.bg.corrected <- if ( all( bg.recon == 0 ) ) {
    raw.data
  } else {
    sweep( raw.data, 2, bg.recon, "-" )
  }

  list( raw.bg.corrected = raw.bg.corrected, median.background = median.background )
}
