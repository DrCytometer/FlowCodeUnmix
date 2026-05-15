# fit_fret_joint.r

#' @title Fit FRET Signals — Joint Proportional Dual-Criterion
#'
#' @description
#' Select and extract the best-fitting FRET correction spectrum per cell using a
#' combined proportional score that balances two complementary error signals,
#' analogous to \code{fit.af.joint()}:
#'
#' \itemize{
#'   \item \strong{Off-channel fluorophore spillover} (\code{p_fluor}): L1 norm
#'     of unmixed values in the fluorophores that are \emph{absent} from the
#'     cell's FlowCode barcode (the "OFF channels"), after subtracting the FRET
#'     shadow.  Expressed as a proportion of the baseline OFF-channel spillover
#'     with no FRET correction at all.
#'   \item \strong{Raw-space residual} (\code{p_resid}): L2 norm of the
#'     detector-space residual after subtracting the fitted FRET signal,
#'     expressed as a proportion of the baseline residual with no FRET
#'     correction.
#' }
#'
#' The combined score for variant \eqn{j} is
#' \deqn{\text{score}_j = \alpha \cdot p_{\text{resid},j}
#'       + (1-\alpha) \cdot p_{\text{fluor},j}}
#' and the variant that minimises this score is selected.  Values below 1
#' indicate an improvement over the no-FRET baseline; values above 1 indicate
#' the variant made things worse.
#'
#' The \code{alpha} parameter controls the trade-off:
#' \code{alpha = 0} is pure off-channel-fluorophore minimisation (equivalent to
#' the original \code{fit.fret()} behaviour modulo the proportional scaling);
#' \code{alpha = 1} is pure raw-residual minimisation;
#' \code{alpha = 0.5} balances both equally (default).
#'
#' @param raw.data_subset Expression data from raw FCS files.  Cells in rows and
#'   detectors in columns.  Columns should be fluorescent data only and must
#'   match the columns in \code{spectra}.  Provide only the subset of data for
#'   cells that carry FlowCodes.
#' @param unmixed_subset Unmixed data (unmixed using fluorophore spectra only,
#'   no FRET correction).  Cells in rows and fluorophores in columns.  Provide
#'   only the FlowCode-carrying subset.
#' @param unmixing.matrix The OLS unmixing matrix for the fluorophore spectra
#'   (\code{solve(spectra \%*\% t(spectra)) \%*\% spectra}).
#' @param spectra Spectral signatures of fluorophores, normalised between 0 and
#'   1, with fluorophores in rows and detectors in columns.
#' @param flowcode.ids_subset The FlowCode barcode ID for each cell in the
#'   subset (integer, matching rows of \code{flowcode.combo.logical}).
#' @param fret.spectra A named list (one element per barcode ID) of matrices
#'   containing the FRET correction spectra for that barcode.  Each matrix has
#'   FRET variants in rows and detectors in columns.
#' @param flowcode.fluors Named character vector mapping ProCode tag names to
#'   fluorophore names.
#' @param flowcode.combo.logical Integer (1/0) matrix with one row per barcode
#'   combo and one column per FlowCode fluorophore, indicating which fluorophores
#'   are active in each combo.
#' @param alpha Numeric in \[0, 1\], default \code{0.5}.  Weight given to the
#'   raw-detector residual term relative to the OFF-channel fluorophore spillover
#'   term.  See Description for details.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{unmixed}}{Matrix (cells × fluorophores) of FRET-corrected
#'     unmixed values.}
#'   \item{\code{fitted.fret}}{Matrix (cells × detectors) of the FRET signal
#'     projected back into raw detector space.}
#' }
#'
#' @seealso [fit.fret()] for the original OFF-channel-only criterion.
#' @seealso [fit.af.joint()] for the analogous AF joint scorer.
#'
#' @export

fit.fret.joint <- function(
    raw.data_subset,
    unmixed_subset,
    unmixing.matrix,
    spectra,
    flowcode.ids_subset,
    fret.spectra,
    flowcode.fluors,
    flowcode.combo.logical,
    alpha = 0.5
) {

  # ---- Setup ------------------------------------------------------------------

  # Strip AF (or any other non-fluorophore column) from spectra and from
  # unmixed_subset before any dimensional arithmetic.
  #
  # backbone$Unmixed is produced by the full unmixing pipeline, which appends
  # an AF column.  spectra has already had AF removed upstream in
  # get.flowcode.spectra.  If the columns of unmixed_subset do not match
  # rownames(spectra) exactly, keep only the shared fluorophore columns so
  # that unmixed %*% spectra is conformable and off-channel indexing is
  # correct.  This mirrors the AF-row guard in fit.af.joint().
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  fluor.names <- rownames( spectra )

  if ( !identical( colnames( unmixed_subset ), fluor.names ) ) {
    shared <- intersect( fluor.names, colnames( unmixed_subset ) )
    if ( length( shared ) < length( fluor.names ) )
      stop(
        "fit.fret.joint: unmixed_subset is missing fluorophore columns: ",
        paste( setdiff( fluor.names, colnames( unmixed_subset ) ), collapse = ", " )
      )
    unmixed_subset <- unmixed_subset[ , fluor.names, drop = FALSE ]
  }

  S           <- t( spectra )                        # detectors × fluorophores
  cell.n      <- nrow( raw.data_subset )
  detector.n  <- ncol( raw.data_subset )
  fluor.n     <- length( fluor.names )
  fc.fluor.idx <- match( flowcode.fluors, fluor.names )

  # Storage — pre-filled with baseline (no-correction) unmixed values so that
  # cells that belong to a placeholder combo are returned as-is.
  final_unmixed <- unmixed_subset   # already fluor.names columns only
  colnames( final_unmixed ) <- fluor.names
  fitted.fret <- matrix(
    0,
    nrow = cell.n,
    ncol = detector.n,
    dimnames = list( NULL, colnames( raw.data_subset ) )
  )

  # ---- Per-barcode loop -------------------------------------------------------

  unique_ids <- unique( flowcode.ids_subset )

  for ( id in unique_ids ) {

    cell_mask <- which( flowcode.ids_subset == id )
    n_cells_id <- length( cell_mask )
    id_raw    <- raw.data_subset[ cell_mask, , drop = FALSE ]
    id_unmixed <- unmixed_subset[ cell_mask, , drop = FALSE ]

    # Identify OFF-channel fluorophore indices for this barcode
    active_fc_mask   <- flowcode.combo.logical[ id, ] == 1
    active_fc_indices <- fc.fluor.idx[ active_fc_mask ]
    off_indices      <- setdiff( seq_len( fluor.n ), active_fc_indices )

    fret_library <- fret.spectra[[ id ]]

    # Skip placeholder combos (single all-zero row written by get.flowcode.spectra)
    if ( nrow( fret_library ) == 1 && all( fret_library == 0 ) ) next

    fret_n <- nrow( fret_library )

    # ---- Pre-compute FRET projection library (shared with fit.fret) -----------
    #
    #   v.library[f, j] = how much FRET variant j leaks into fluorophore f
    #   r.library[d, j] = detector-space residual of variant j that fluorophore
    #                     spectra cannot explain (orthogonal complement)
    #
    v.library <- unmixing.matrix %*% t( fret_library )   # fluor.n × fret_n
    r.library <- t( fret_library ) - ( S %*% v.library )  # detector.n × fret_n

    # Optimal FRET intensity k[i, j] = (raw_i · r_j) / (r_j · r_j), clamp >= 0
    numerator   <- id_raw %*% r.library                   # n_cells_id × fret_n
    denominator <- colSums( r.library^2 )
    denominator[ denominator == 0 ] <- 1e-10
    k.matrix <- sweep( numerator, 2, denominator, "/" )
    k.matrix[ k.matrix < 0 ] <- 0

    # ---- Baseline errors (denominators for proportional scores) ---------------
    #
    # Fluorophore baseline: OFF-channel spillover with no FRET correction.
    # Use unclamped unmixed values (negative values are real unmixing artifacts
    # and should be counted in the L1 penalty).
    base.e.fluor <- rowSums( abs( id_unmixed[ , off_indices, drop = FALSE ] ) ) + 1e-6

    # Raw-space baseline: residual of the non-negative-clamped fluorophore fit.
    unmixed.nonneg <- id_unmixed
    unmixed.nonneg[ unmixed.nonneg < 0 ] <- 0
    resids.initial <- id_raw - ( unmixed.nonneg %*% spectra )  # n_cells_id × detector.n
    base.e.resid <- sqrt( rowSums( resids.initial^2 ) ) + 1e-6

    # ---- Per-variant proportional scoring (vectorised over cells) -------------

    prop.score.matrix <- matrix( 0, nrow = n_cells_id, ncol = fret_n )

    for ( j in seq_len( fret_n ) ) {

      k_j <- k.matrix[ , j ]             # length n_cells_id

      # OFF-channel fluorophore spillover after FRET correction
      v_j_off    <- v.library[ off_indices, j ]            # length |off_indices|
      adj.off    <- id_unmixed[ , off_indices, drop = FALSE ] -
        tcrossprod( k_j, v_j_off )           # n_cells_id × |off|
      e.fluor.j  <- rowSums( abs( adj.off ) )

      # Raw-space residual after subtracting the fitted FRET signal
      r_j        <- r.library[ , j ]                       # length detector.n
      adj.resids <- resids.initial - tcrossprod( k_j, r_j ) # n_cells_id × detector.n
      e.resid.j  <- sqrt( rowSums( adj.resids^2 ) )

      # Proportional improvement (< 1 = better than no FRET, > 1 = worse)
      p.fluor <- e.fluor.j / base.e.fluor
      p.resid <- e.resid.j / base.e.resid

      prop.score.matrix[ , j ] <- alpha * p.resid + ( 1 - alpha ) * p.fluor
    }

    # Select the variant that minimises the combined proportional score per cell
    best_fret_indices <- max.col( -prop.score.matrix, ties.method = "first" )

    # ---- Final reconstruction (analogous to fit.fret) -------------------------

    for ( j in seq_len( fret_n ) ) {
      sub_mask <- which( best_fret_indices == j )
      if ( length( sub_mask ) == 0 ) next

      g_idx       <- cell_mask[ sub_mask ]
      ki_final    <- k.matrix[ sub_mask, j, drop = FALSE ]  # n_sub × 1
      vi_final    <- v.library[ , j, drop = FALSE ]          # fluor.n × 1
      spec_final  <- fret_library[ j, , drop = FALSE ]       # 1 × detector.n

      # Fluorophore-space correction
      final_unmixed[ g_idx, ] <- unmixed_subset[ g_idx, ] -
        tcrossprod( ki_final, vi_final )

      # Raw-space FRET signal reconstruction
      fitted.fret[ g_idx, ] <- ki_final %*% spec_final
    }
  }

  return( list(
    unmixed     = final_unmixed,
    fitted.fret = fitted.fret
  ) )
}
