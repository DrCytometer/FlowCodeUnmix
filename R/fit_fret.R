# fit_fret.r

#' @title Fit FRET Signals
#'
#' @description
#' Select and extract the best fitting FRET correction spectrum per cell. Uses
#' minimization of fluorophore spillover in the channels that do not belong to
#' that cell's FlowCode barcode as the criterion.
#'
#' @param raw.data_subset Expression data from raw FCS files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra. Only provide the subset of data for the cells
#' that have FlowCodes.
#' @param unmixed_subset Unmixed data (unmixed using fluorophore spectra only).
#' Cells in rows and fluorophores in columns. Only provide the subset of data
#' for the cells that have FlowCodes.
#' @param unmixing.matrix The OLS unmixing matrix for the fluorophore spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param flowcode.ids_subset The FlowCode barcode IDs for each cell.
#' @param fret.spectra A list of matrices containing the FRET correction spectra
#' per barcode.
#' @param flowcode.fluors The names of the fluorophores used to identify the
#' FlowCodes
#' @param flowcode.combo.logical Integer (1/0) matrix describing the combinations
#' of fluorophores in each FlowCode barcode.
#'
#' @return A list containing the adjusted unmixed data (`unmixed`) and the raw
#' space projection for the FRET (`fitted.fret`).
#'
#' @export

fit.fret <- function(
    raw.data_subset,
    unmixed_subset,
    unmixing.matrix,
    spectra,
    flowcode.ids_subset,
    fret.spectra,
    flowcode.fluors,
    flowcode.combo.logical
  ) {

  # transpose spectra & define quantities
  S <- t(spectra)
  cell.n <- nrow(raw.data_subset)
  detector.n <- ncol(raw.data_subset)
  fluor.names <- rownames(spectra)
  fluor.n <- length(fluor.names)
  fc.fluor.idx <- match(flowcode.fluors, fluor.names)

  # Storage
  final_unmixed <- matrix(0, nrow = cell.n, ncol = fluor.n, dimnames = list(NULL, fluor.names))
  fitted.fret <- matrix(0, nrow = cell.n, ncol = detector.n, dimnames = list(NULL, colnames(raw.data_subset)))

  # cycle through each FlowCode ID to correct FRET per barcode
  unique_ids <- unique(flowcode.ids_subset)

  for (id in unique_ids) {
    cell_mask <- which(flowcode.ids_subset == id)
    id_raw <- raw.data_subset[cell_mask, , drop = FALSE]

    # Identify OFF-channels
    active_fc_mask <- flowcode.combo.logical[id, ] == 1
    active_fc_indices <- fc.fluor.idx[active_fc_mask]
    off_indices <- setdiff(seq_len(fluor.n), active_fc_indices)

    fret_library <- fret.spectra[[id]]
    fret_n <- nrow(fret_library)

    # Project FRET library
    v.library <- unmixing.matrix %*% t(fret_library)
    r.library <- t(fret_library) - (S %*% v.library)

    # Calculate k (FRET intensity)
    numerator <- id_raw %*% r.library
    denominator <- colSums(r.library^2)
    denominator[denominator == 0] <- 1e-10
    k.matrix <- scale(numerator, center = FALSE, scale = denominator)
    k.matrix[k.matrix < 0] <- 0

    # Optimization (L1 on OFF-channels)
    best_errors <- rep(Inf, length(cell_mask))
    best_fret_indices <- integer(length(cell_mask))
    id_unmixed_base_off <- unmixed_subset[cell_mask, off_indices, drop = FALSE]

    for (i in seq_len(fret_n)) {
      ki <- k.matrix[, i]
      vi_off <- v.library[off_indices, i]
      current_errors <- rowSums(abs(id_unmixed_base_off - (ki %*% t(vi_off))))

      better <- current_errors < best_errors
      best_errors[better] <- current_errors[better]
      best_fret_indices[better] <- i
    }

    # Final Extraction of Unmixed and Raw-Space FRET
    for (i in unique(best_fret_indices)) {
      sub_mask <- which(best_fret_indices == i)
      g_idx <- cell_mask[sub_mask]

      ki_final <- k.matrix[sub_mask, i, drop = FALSE]
      vi_final <- v.library[, i, drop = FALSE]
      spec_final <- fret_library[i, , drop = FALSE] # The raw detector shape

      # Fluorophore space correction
      final_unmixed[g_idx, ] <- unmixed_subset[g_idx, ] - (ki_final %*% t(vi_final))

      # Raw space FRET signal reconstruction
      fitted.fret[g_idx, ] <- ki_final %*% spec_final
    }
  }

  return(list(
    unmixed = final_unmixed,
    fitted.fret = fitted.fret
  ))
}
