# fit_af.r

#' @title Fit Autofluorescence
#'
#' @description
#' Select and extract the best fitting autofluorescence spectrum per cell. Uses
#' minimization of fluorophore spillover as the criterion.
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param unmixed Unmixed data (unmixed using fluorophore spectra only). Cells
#' in rows and fluorophores in columns.
#' @param unmixing.matrix The OLS unmixing matrix for the fluorophore spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#'
#' @return A list containing the adjusted unmixed data (`unmixed`), the abundances
#' of autofluorescence per cell (`AF`), the index corresponding to the spectrum
#' selected per cell (`af.idx`) and the raw space projection for the
#' autofluorescence (`fitted.af`).
#'
#' @export

fit.af <- function(
    raw.data,
    unmixed,
    unmixing.matrix,
    spectra,
    af.spectra
) {

  if ("AF" %in% rownames(spectra))
    spectra <- spectra[rownames(spectra) != "AF", , drop = FALSE]

  cell.n <- nrow(raw.data)
  detector.n <- ncol(raw.data)
  af.n <- nrow(af.spectra)
  fluorophore.n <- nrow(spectra)
  S <- t(spectra)

  # How much each AF variant 'looks like' each fluorophore
  v.library <- unmixing.matrix %*% t(af.spectra)

  # Calculate the 'residual AF' (the part of AF fluors can't explain)
  r.library <- t(af.spectra) - (S %*% v.library)

  # Predicted AF intensity (k)
  numerator <- raw.data %*% r.library
  denominator <- colSums(r.library^2)
  denominator[denominator == 0] <- 1e-10
  k.matrix <- sweep(numerator, 2, denominator, "/")
  k.matrix[k.matrix < 0] <- 0 # Physical constraint: positive AF intensity

  # Optimization: Pick AF variant that minimizes signal in fluorophore channels
  best.errors <- rep(Inf, cell.n)
  af.idx <- integer(cell.n)

  for (i in seq_len(af.n)) {
    ki <- k.matrix[, i]
    vi <- v.library[, i]

    # Current unmixed signal minus the AF "shadow" cast into fluorophore space
    # (ki %*% t(vi)) creates the [cell.n x fluor.n] correction matrix
    current.errors <- rowSums(abs(unmixed - (ki %*% t(vi))))

    improved <- current.errors < best.errors
    best.errors[improved] <- current.errors[improved]
    af.idx[improved] <- i
  }

  # Final Reconstruction
  final.fluors <- matrix(0, nrow = cell.n, ncol = fluorophore.n)
  fitted.af <- matrix(0, nrow = cell.n, ncol = detector.n)
  best.af.intensity <- numeric(cell.n)

  for (i in seq_len(af.n)) {
    idx <- which(af.idx == i)
    if (length(idx) > 0) {
      ki_final <- k.matrix[idx, i, drop = FALSE]
      vi_final <- v.library[, i, drop = FALSE]
      spec_final <- af.spectra[i, , drop = FALSE]

      # Correct fluorophore values
      final.fluors[idx, ] <- unmixed[idx, ] - (ki_final %*% t(vi_final))

      # Reconstruct raw detector space AF signal
      fitted.af[idx, ] <- ki_final %*% spec_final

      # Store the intensity (k)
      best.af.intensity[idx] <- ki_final
    }
  }

  colnames(final.fluors) <- rownames(spectra)
  colnames(fitted.af) <- colnames(raw.data)

  results <- list(
    unmixed = final.fluors,
    AF = best.af.intensity,
    af.idx = af.idx,
    fitted.af = fitted.af
  )

  return(results)
}
