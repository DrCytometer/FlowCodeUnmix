# fit_af_joint.r

#' @title Fit Autofluorescence Joint (Proportional Dual-Criterion)
#'
#' @description
#' Selects the best-fitting autofluorescence spectrum per cell using a combined
#' proportional score that balances two complementary error signals:
#'
#' \itemize{
#'   \item \strong{Fluorophore spillover} (`p_fluor`): L1 norm of unmixed
#'     fluorophore values after subtracting the AF shadow, expressed as a
#'     proportion of the baseline spillover with no AF correction at all.
#'   \item \strong{Raw-space residual} (`p_resid`): L2 norm of the detector-space
#'     residual after subtracting the fitted AF signal, expressed as a proportion
#'     of the baseline residual with no AF correction.
#' }
#'
#' The combined score for variant \eqn{j} is
#' \deqn{\text{score}_j = \alpha \cdot p_{\text{resid},j} + (1-\alpha) \cdot p_{\text{fluor},j}}
#' and the variant that minimises this score is selected. Values below 1 indicate
#' an improvement over the no-AF baseline; values above 1 indicate the variant
#' made things worse.
#'
#' The `alpha` parameter controls the trade-off: `alpha = 0` is pure fluorophore
#' minimisation (equivalent to the original `fit.af()` behaviour modulo the
#' proportional scaling); `alpha = 1` is pure raw-residual minimisation.
#' `alpha = 0.5` balances both equally.
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
#'   detectors in columns. Columns should be fluorescent data only and must
#'   match the columns in `spectra`.
#' @param unmixed Unmixed data (unmixed using fluorophore spectra only, no AF).
#'   Cells in rows and fluorophores in columns.
#' @param unmixing.matrix The OLS unmixing matrix for the fluorophore spectra
#'   (`solve(spectra %*% t(spectra)) %*% spectra`).
#' @param spectra Spectral signatures of fluorophores, normalised between 0
#'   and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalised
#'   between 0 and 1, with autofluorescence variants in rows and detectors in
#'   columns. Prepare using `get.af.spectra`.
#' @param alpha Numeric in \[0, 1\], default `0.5`. Weight given to the
#'   raw-detector residual term relative to the fluorophore spillover term.
#'   See Description for details.
#'
#' @return A list with four elements:
#' \describe{
#'   \item{`unmixed`}{Matrix (cells × fluorophores) of AF-corrected unmixed
#'     values.}
#'   \item{`AF`}{Numeric vector (length `nrow(raw.data)`) of per-cell AF
#'     intensities (the fitted scalar \eqn{k}).}
#'   \item{`af.idx`}{Integer vector of the AF variant index selected per cell
#'     (1-based).}
#'   \item{`fitted.af`}{Matrix (cells × detectors) of the AF signal projected
#'     back into raw detector space, for subtraction from raw data downstream.}
#' }
#'
#' @seealso [fit.af()] for the original fluorophore-only criterion.
#'
#' @export

fit.af.joint <- function(
    raw.data,
    unmixed,
    unmixing.matrix,
    spectra,
    af.spectra,
    alpha = 0.5
) {

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  cell.n       <- nrow( raw.data )
  detector.n   <- ncol( raw.data )
  af.n         <- nrow( af.spectra )
  fluorophore.n <- nrow( spectra )
  S            <- t( spectra )   # detectors × fluorophores

  # ---- Pre-compute AF projection library (shared with original fit.af) -------

  # v.library[f, j] = how much AF variant j leaks into fluorophore f
  v.library <- unmixing.matrix %*% t( af.spectra )   # fluorophores × af.n

  # r.library[d, j] = the part of AF variant j that fluorophore spectra cannot
  #                   explain (orthogonal complement in detector space)
  r.library <- t( af.spectra ) - ( S %*% v.library ) # detectors × af.n

  # Solve for optimal AF intensity k[i, j] for every cell i and variant j
  # k = (raw · r_j) / (r_j · r_j)   [unconstrained OLS for the AF scalar]
  numerator   <- raw.data %*% r.library              # cell.n × af.n
  denominator <- colSums( r.library^2 )
  denominator[ denominator == 0 ] <- 1e-10
  k.matrix <- sweep( numerator, 2, denominator, "/" )
  k.matrix[ k.matrix < 0 ] <- 0   # physical constraint: AF intensity >= 0

  # ---- Baseline errors (denominators for proportional scores) ----------------

  # Fluorophore baseline: spillover when no AF is extracted at all.
  # Use non-negative-clamped unmixed values so that the baseline fitted signal
  # matches what would actually be projected back (negative unmixed values do
  # not contribute real spectral signal).
  unmixed.nonneg    <- unmixed
  unmixed.nonneg[ unmixed.nonneg < 0 ] <- 0
  resids.initial    <- raw.data - ( unmixed.nonneg %*% spectra ) # cells × detectors

  base.e.fluor <- rowSums( abs( unmixed ) ) + 1e-6  # L1 fluorophore, unclamped
  base.e.resid <- sqrt( rowSums( resids.initial^2 ) ) + 1e-6  # L2 raw residual

  # ---- Per-variant proportional scoring (vectorised over cells) --------------

  prop.score.matrix <- matrix( 0, nrow = cell.n, ncol = af.n )

  for ( j in seq_len( af.n ) ) {

    k_j <- k.matrix[ , j ]          # length cell.n

    # Fluorophore spillover after AF correction:
    # unmixed_corrected = unmixed - k_j * v_j
    v_j        <- v.library[ , j ]                       # length fluorophore.n
    adj.fluors <- unmixed - tcrossprod( k_j, v_j )       # cell.n × fluorophore.n
    e.fluor.j  <- rowSums( abs( adj.fluors ) )

    # Raw-space residual after subtracting the fitted AF signal:
    # resid_corrected = resids_initial - k_j * r_j
    r_j        <- r.library[ , j ]                       # length detector.n
    adj.resids <- resids.initial - tcrossprod( k_j, r_j ) # cell.n × detectors
    e.resid.j  <- sqrt( rowSums( adj.resids^2 ) )

    # Proportional improvement (< 1 = better than no AF, > 1 = worse)
    p.fluor <- e.fluor.j / base.e.fluor
    p.resid <- e.resid.j / base.e.resid

    prop.score.matrix[ , j ] <- alpha * p.resid + ( 1 - alpha ) * p.fluor
  }

  # Select the variant that minimises the combined proportional score per cell
  af.idx <- max.col( -prop.score.matrix, ties.method = "first" )

  # ---- Final reconstruction (identical to original fit.af) -------------------

  final.fluors      <- matrix( 0, nrow = cell.n, ncol = fluorophore.n )
  fitted.af         <- matrix( 0, nrow = cell.n, ncol = detector.n )
  best.af.intensity <- numeric( cell.n )

  for ( j in seq_len( af.n ) ) {
    idx <- which( af.idx == j )
    if ( length( idx ) > 0 ) {
      ki_final   <- k.matrix[ idx, j, drop = FALSE ]         # n_idx × 1
      vi_final   <- v.library[ , j, drop = FALSE ]           # fluorophore.n × 1
      spec_final <- af.spectra[ j, , drop = FALSE ]          # 1 × detector.n

      final.fluors[ idx, ]      <- unmixed[ idx, ] - tcrossprod( ki_final, vi_final )
      fitted.af[ idx, ]         <- ki_final %*% spec_final
      best.af.intensity[ idx ]  <- ki_final
    }
  }

  colnames( final.fluors ) <- rownames( spectra )
  colnames( fitted.af )    <- colnames( raw.data )

  return( list(
    unmixed   = final.fluors,
    AF        = best.af.intensity,
    af.idx    = af.idx,
    fitted.af = fitted.af
  ) )
}
