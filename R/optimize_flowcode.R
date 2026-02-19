# optimize_flowcode.r

#' @title Optimize FlowCode spectral unmixing
#'
#' @description
#' Parallel backend for FlowCode-based spectral optimization in R.
#'
#' @importFrom AutoSpectral create.parallel.lapply parallel.backend
#' @importFrom AutoSpectral unmix.ols.fast
#' @importFrom parallelly availableCores
#' @importFrom RhpcBLASctl blas_set_num_threads omp_set_num_threads
#'
#' @param raw.data Numeric matrix (n_cells x n_detectors)
#' @param unmixed Numeric matrix (n_cells x n_fluors)
#' @param spectra Numeric matrix (n_fluors x n_detectors)
#' @param pos.thresholds Numeric vector (n_fluors)
#' @param flowcode.ids Integer vector (n_cells), 0 = no flowcode
#' @param has.flowcode Integer/logical vector (n_cells)
#' @param flowcode.combo.logical Integer matrix (n_combos x n_flowcode.fluors)
#' @param flowcode.fluors Character vector
#' @param optimize.fluors Character vector
#' @param variants List of variant matrices per fluorophore
#' @param delta.list List of delta matrices per fluorophore
#' @param delta.norms List of delta norms per fluorophore
#' @param asp The AutoSpectral parameter list.
#' @param k Integer, number of variants to test
#' @param threads Integer, number of threads
#' @param parallel Logical, whether to use parallel processing
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

optimize.flowcode <- function(
    raw.data,
    unmixed,
    spectra,
    pos.thresholds,
    flowcode.ids,
    has.flowcode,
    flowcode.combo.logical,
    flowcode.fluors,
    optimize.fluors,
    variants,
    delta.list,
    delta.norms,
    asp,
    k = 10L,
    threads = 1L,
    parallel = TRUE
) {

  cell.n <- nrow( raw.data )

  # set up parallel backend
  if ( parallel ) {

    result <- AutoSpectral::create.parallel.lapply(
      asp,
      # modify exports as needed
      exports = c(
        "raw.data", "unmixed", "spectra", "pos.thresholds","flowcode.ids",
        "has.flowcode", "flowcode.combo.logical", "flowcode.fluors",
        "optimize.fluors", "variants", "delta.list", "delta.norms", "k"
      ),
      parallel = TRUE,
      threads = threads,
      export.env = environment()
    )
    lapply.function <- result$lapply
  } else {
    lapply.function <- lapply
    result <- list( cleanup = NULL )
  }

  # loop over each cell, optimizing fluorophore spectra
  unmixed.opt <- tryCatch(
    expr = {
      lapply.function( seq_len( cell.n ), function( cell ) {

        # get cell's data
        cell.raw <- raw.data[ cell, , drop = FALSE ]
        cell.unmixed <- unmixed[ cell, , drop = FALSE ]

        # check whether this cell has any fluorophores present
        pos.fluors <- as.numeric( cell.unmixed ) >= pos.thresholds

        if ( !any( pos.fluors ) ) return( cell.unmixed )

        if ( has.flowcode[ cell ] ) {
          # which FlowCode combo is this?
          id <- flowcode.ids[ cell ]
          # which FlowCode fluors are expected in this cell?
          allowed.idx <- flowcode.fluors[ flowcode.combo.logical[ id, ] == 1 ]
          pos.fluors[ flowcode.fluors ] <- FALSE
          pos.fluors[ allowed.idx ] <- TRUE
        }

        # set baseline spectra
        cell.spectra.final <- spectra

        # remove absent fluorophores for optimization, unmix
        cell.spectra.curr <- cell.spectra.final[ pos.fluors, , drop = FALSE ]
        cell.unmixed <- AutoSpectral::unmix.ols.fast( cell.raw, cell.spectra.curr )

        # set baseline unmixed and residuals
        resid <- cell.raw - ( cell.unmixed %*% cell.spectra.curr )
        error.final <- sum( abs( resid ) )


        ########################################
        ### Fluorophore Optimization Section ###
        ########################################

        # restrict optimization to fluorophores we have variants for
        ### TBD: use indexing rather than %in%
        fluors.to.sort <- optimize.fluors[
          optimize.fluors %in% names( pos.fluors )[ pos.fluors ] ]

        if ( length( fluors.to.sort ) > 0 ) {
          # sort by abundance to optimize brightest fluors first (error is proportional to signal)
          fluor.order <- sort( cell.unmixed[ , fluors.to.sort ], decreasing = TRUE )

          for ( fl in names( fluor.order ) ) {
            fl.variants <- variants[[ fl ]]
            delta.fl <- delta.list[[ fl ]]
            delta.norm  <- delta.norms[[ fl ]]

            # score variants
            joint.score <- as.numeric( delta.fl %*% t( resid ) ) * cell.unmixed[ , fl ]
            joint.score <- joint.score / delta.norm / sqrt( sum( resid^2 ) )

            # select k variants up to the max we have available
            k.eff <- min( k, length( joint.score ) )
            topK <- order( joint.score, decreasing = TRUE )[ seq_len( k.eff ) ]

            # test the top k scoring variants
            for ( var in topK ) {
              # supplant the base spectrum with this variant
              cell.spectra.curr[ fl, ] <- fl.variants[ var, ]

              # reunmix with this variant
              trial.unmix <- AutoSpectral::unmix.ols.fast( cell.raw, cell.spectra.curr )

              # assess the residual error with this variant
              trial.resid <- cell.raw - ( trial.unmix %*% cell.spectra.curr )
              trial.error <- sum( abs( trial.resid ) )

              # accept change if residual is lower
              if ( trial.error < error.final ) {
                error.final <- trial.error
                cell.spectra.final[ fl, ] <- cell.spectra.curr[ fl, ]
                resid <- trial.resid
              } else {
                # reject if not
                cell.spectra.curr[ fl, ] <- cell.spectra.final[ fl, ]
              }
            }
          }
        }


        ##############################################################
        ### Final Unmix Using Optimized Spectra (All Fluorophores) ###
        ##############################################################

        cell.unmixed <- AutoSpectral::unmix.ols.fast( cell.raw, cell.spectra.final )

        return( cell.unmixed )
      } )
    },
    finally = {
      if ( !is.null( result$cleanup ) ) result$cleanup()
    }
  )

  # combine data into a matrix
  unmixed.opt <- do.call( rbind, unmixed.opt )

  return( unmixed.opt )
}
