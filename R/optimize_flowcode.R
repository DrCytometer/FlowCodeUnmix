# optimize_flowcode.r

#' @title Optimize FlowCode spectral unmixing
#'
#' @description
#' Parallel backend for FlowCode-based spectral optimization in R.
#'
#' @param raw_data Numeric matrix (n_cells x n_detectors)
#' @param unmixed Numeric matrix (n_cells x n_fluors)
#' @param unmix_fun Unmixing function (WLS)
#' @param combined_spectra Numeric matrix (n_fluors x n_detectors)
#' @param weights Numeric vector (n_detectors)
#' @param pos_thresholds Numeric vector (n_fluors)
#' @param af_idx Integer vector (n_cells), 1-based AF indices
#' @param af_spectra Numeric matrix (n_af_variants x n_detectors)
#' @param flowcode_ids Integer vector (n_cells), 0 = no flowcode
#' @param has_flowcode Integer/logical vector (n_cells)
#' @param combo_fret List of FRET spectra matrices
#' @param fret_delta_list List of delta matrices
#' @param fret_delta_norms List of norm vectors
#' @param flowcode_combo_logical Integer matrix (n_combos x n_flowcode_fluors)
#' @param flowcode_fluors Character vector
#' @param optimize_fluors Character vector
#' @param variants List of variant matrices per fluorophore
#' @param delta_list List of delta matrices per fluorophore
#' @param delta_norms List of delta norms per fluorophore
#' @param fluorophores Character vector of fluorophore names
#' @param af_idx_in_spectra Integer, 0-based index of AF in spectra
#' @param asp The AutoSpectral parameter list.
#' @param k Integer, number of variants to test
#' @param cell_weighting Logical, use cell-specific weights
#' @param cell_weight_regularize Logical, regularize cell weights
#' @param nthreads Integer, number of threads
#' @param parallel Logical, whether to use parallel processing
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

optimize.flowcode <- function(
    raw_data,
    unmixed,
    unmix_fun,
    combined_spectra,
    weights,
    pos_thresholds,
    af_idx,
    af_spectra,
    flowcode_ids,
    has_flowcode,
    combo_fret,
    fret_delta_list,
    fret_delta_norms,
    flowcode_combo_logical,
    flowcode_fluors,
    optimize_fluors,
    variants,
    delta_list,
    delta_norms,
    fluorophores,
    af_idx_in_spectra,
    asp,
    k = 10L,
    cell_weighting = FALSE,
    cell_weight_regularize = FALSE,
    nthreads = 1L,
    parallel = TRUE
) {

  # this would probably be faster if we moved AF to position 1

  cell.n <- nrow( raw_data )

  # set up parallel backend
  if ( parallel ) {

    result <- create.parallel.lapply( # call from AutoSpectral
      asp,
      # modify exports as needed
      exports = c(
        "raw_data", "unmixed", "unmix_fun",
        "combined_spectra", "weights", "pos_thresholds",
        "af_idx", "af_spectra", "flowcode_ids",
        "has_flowcode", "combo_fret", "fret_delta_list", "fret_delta_norms",
        "flowcode_combo_logical", "flowcode_fluors", "optimize_fluors",
        "variants", "delta_list", "delta_norms", "fluorophores",
        "af_idx_in_spectra", "k", "cell_weighting", "cell_weight_regularize"
      ),
      parallel = TRUE,
      threads = nthreads,
      export.env = environment()
    )
    lapply.function <- result$lapply
  } else {
    lapply.function <- lapply
    result <- list( cleanup = NULL )
  }

  # loop over each cell, optimizing fluorophore spectra
  # also extracting FRET errors if the cell has a FlowCode combination
  unmixed.opt <- tryCatch(
    expr = {
      lapply.function( seq_len( cell.n ), function( cell ) {

        # get cell's data
        cell.raw <- raw_data[ cell, , drop = FALSE ]
        cell.unmixed <- unmixed[ cell, , drop = FALSE ]

        # set weights in a cell-specific manner
        cell.weights <- weights

        if ( cell_weighting ) {
          # use cell-specific weighting (Poisson-like)
          w <- abs( cell.raw )
          w[ w < 1e-6 ] <- 1e-6
          cell.weights <- 1 / w

          if ( cell_weight_regularize ) {
            # regularize weight towards weights for full data
            cell.weights <- ( cell.weights + weights ) * 0.5
          }
        }

        # determine which AF has been selected
        cell.af.idx <- af_idx[ cell ]
        cell.af <- af_spectra[ cell.af.idx, ]

        # calculate delta matrix and norm for this AF spectrum
        af.delta <- sweep( af_spectra, 2, cell.af, "-" )
        af.delta.norm <- sqrt( rowSums( af.delta^2 ) )

        # set baseline spectra
        cell.spectra.final <- combined_spectra
        cell.spectra.final[ af_idx_in_spectra, ] <- cell.af

        # check whether this cell has any fluorophores present
        pos.fluors <- as.numeric( cell.unmixed ) >= pos_thresholds

        # remove absent fluorophores for optimization, unmix
        cell.spectra.curr <- cell.spectra.final[ which( pos.fluors ), , drop = FALSE ]
        cell.unmixed <- unmix_fun( cell.raw, cell.spectra.curr, cell.weights )

        # set baseline unmixed and residuals
        resid <- cell.raw - ( cell.unmixed %*% cell.spectra.curr )
        error.final <- sum( abs( resid ) )

        ###############################
        ### Early Exit (AF Only)    ###
        ###############################

        # if nothing is positive, return early--change to reset AF here
        if ( !any( pos.fluors[ fluorophores ] ) ) {

          # score AF variants for alignment with the residual
          joint.score <- as.numeric( af.delta %*% t( resid ) ) * cell.unmixed[ , "AF" ]
          joint.score <- joint.score / af.delta.norm / sqrt( sum( resid^2 ) )

          # pick top k candidates to test
          k.eff <- min( k, length( joint.score ) )
          topK <- order( joint.score, decreasing = TRUE )[ seq_len( k.eff ) ]

          # test the top k scoring variants
          for ( var in topK ) {
            # supplant the base spectrum with this variant
            cell.spectra.curr[ "AF", ] <- af_spectra[ var, ]

            # reunmix with this variant
            trial.unmix <- unmix_fun( cell.raw, cell.spectra.curr, cell.weights )

            # assess the residual error with this variant
            trial.resid <- cell.raw - ( trial.unmix %*% cell.spectra.curr )
            trial.error <- sum( abs( trial.resid ) )

            # accept change if residual is lower
            if ( trial.error < error.final ) {
              error.final <- trial.error
              cell.spectra.final[ af_idx_in_spectra, ] <- af_spectra[ var, ]
              resid <- trial.resid
            }
          }

          cell.unmixed <- unmix_fun( cell.raw, cell.spectra.final, cell.weights )

          return( cell.unmixed )
        }


        ###############################
        ### FRET Correction Section ###
        ###############################

        ### check whether this cell has a valid FlowCode combo
        # if it does, remove any non-combo FlowCode fluors (incorrect) for the
        # optimization steps.
        if ( has_flowcode[ cell ] ) {
          # which FlowCode combo is this?
          id <- flowcode_ids[ cell ]

          # use alignment between delta and residuals to select best FRET spectrum
          variants.fr <- combo_fret[[ id ]]
          delta.fr <- fret_delta_list[[ id ]]
          delta.norm <- fret_delta_norms[[ id ]]

          # remove non-combo FlowCode fluors from those considered present
          disallowed.idx <- flowcode_fluors[ flowcode_combo_logical[id, ] == 0 ]
          pos.fluors[ disallowed.idx ] <- FALSE
          # drop non-combo FlowCode fluors from spectra for testing
          cell.spectra.curr <- cell.spectra.final[ which( pos.fluors ), , drop = FALSE ]

          # unmix with median FRET spectrum
          cell.spectra.fret <- rbind(
            variants.fr[ 1, ], # first row is median
            cell.spectra.curr
          )
          trial.unmix <-  unmix_fun( cell.raw, cell.spectra.fret, cell.weights )

          # check if this lowers the residual
          trial.resid <- cell.raw - ( trial.unmix %*% cell.spectra.fret )
          trial.error <- sum( abs( trial.resid ) )

          # set this as the target to beat
          resid <- trial.resid
          error.final <- trial.error
          fitted.fret <- trial.unmix[ , 1, drop = FALSE ] %*%
            cell.spectra.fret[ 1, , drop = FALSE ]

          # score FRET variants based on alignment to the residual
          joint.score <- as.numeric( delta.fr %*% t( resid ) ) * trial.unmix[ , 1 ]
          joint.score <- joint.score / delta.norm / sqrt( sum( resid^2 ) )

          # select number of high-scoring variants to test
          k.eff <- min( k, length( joint.score ) )
          topK <- order( joint.score, decreasing = TRUE )[ seq_len( k.eff ) ]

          # test the top k scoring variants sequentially
          for ( var in topK ) {
            # supplant the base spectrum with this variant
            cell.spectra.fret[ 1, ] <- variants.fr[ var, ]

            # reunmix with this variant
            trial.unmix <- unmix_fun( cell.raw, cell.spectra.fret, cell.weights )

            # assess the residual error with this variant
            trial.resid <- cell.raw - ( trial.unmix %*% cell.spectra.fret )
            trial.error <- sum( abs( trial.resid ) )

            # accept change if residual is lower
            if ( trial.error < error.final ) {
              error.final <- trial.error
              resid <- trial.resid
              fitted.fret <- trial.unmix[ , 1, drop = FALSE ] %*%
                cell.spectra.fret[ 1, , drop = FALSE ]
            }
          }

          # subtract FRET from cell.raw for fluorophore optimization
          cell.raw <- cell.raw - fitted.fret

          # recompute baseline unmix and residuals after FRET correction
          cell.unmixed <- unmix_fun( cell.raw, cell.spectra.curr, cell.weights )
          resid <- cell.raw - ( cell.unmixed %*% cell.spectra.curr )
          error.final <- sum( abs( resid ) )
        }


        #############################################
        ### Autofluorescence Optimization Section ###
        #############################################

        # score AF spectra for alignment with residual
        joint.score <- as.numeric( af.delta %*% t( resid ) ) * cell.unmixed[ , "AF" ]
        joint.score <- joint.score / af.delta.norm / sqrt( sum( resid^2 ) )

        # pick top k candidates to test
        k.eff <- min( k, length( joint.score ) )
        topK <- order( joint.score, decreasing = TRUE )[ seq_len( k.eff ) ]

        # test the top k scoring variants
        for ( var in topK ) {
          # supplant the base spectrum with this variant
          cell.spectra.curr[ "AF", ] <- af_spectra[ var, ]

          # reunmix with this variant
          trial.unmix <- unmix_fun( cell.raw, cell.spectra.curr, cell.weights )

          # assess the residual error with this variant
          trial.resid <- cell.raw - ( trial.unmix %*% cell.spectra.curr )
          trial.error <- sum( abs( trial.resid ) )

          # accept change if residual is lower
          if ( trial.error < error.final ) {
            error.final <- trial.error
            cell.spectra.final[ af_idx_in_spectra, ] <- af_spectra[ var, ]
            resid <- trial.resid
          } else {
            # reject if not
            cell.spectra.curr[ "AF", ] <- cell.spectra.final[ af_idx_in_spectra, ]
          }
        }

        # recompute baseline unmix after AF correction
        cell.unmixed <- unmix_fun( cell.raw, cell.spectra.final, cell.weights )

        # reassess fluorophore positivity--this may not be really necessary
        pos.fluors <- as.numeric( cell.unmixed ) >= pos_thresholds

        # if it is necessary, we can again drop out cells with only AF remaining
        if ( !any( pos.fluors[ fluorophores ] ) ) return( cell.unmixed )

        if ( has_flowcode[ cell ] ) {
          # which FlowCode fluors are expected in this cell?
          allowed.idx <- flowcode_fluors[ flowcode_combo_logical[ id, ] == 1 ]
          pos.fluors[ flowcode_fluors ] <- FALSE
          pos.fluors[ allowed.idx ] <- TRUE
        }

        # reset baseline spectra
        cell.spectra.curr <- cell.spectra.final[ which( pos.fluors ), , drop = FALSE ]

        # reset lower dimensional unmixing and targets
        cell.unmixed <- unmix_fun( cell.raw, cell.spectra.curr, cell.weights )
        resid <- cell.raw - ( cell.unmixed %*% cell.spectra.curr )
        error.final <- sum( abs( resid ) )


        ########################################
        ### Fluorophore Optimization Section ###
        ########################################

        # restrict optimization to fluorophores we have variants for
        ### TBD: use indexing rather than %in%
        fluors.to.sort <- optimize_fluors[ optimize_fluors %in%
                                             names( pos.fluors )[ pos.fluors ] ]

        if ( length( fluors.to.sort ) > 0 ) {
          # sort by abundance to optimize brightest fluors first (error is proportional to signal)
          fluor.order <- sort( cell.unmixed[ , fluors.to.sort ], decreasing = TRUE )

          for ( fl in names( fluor.order ) ) {
            fl.variants <- variants[[ fl ]]
            delta.fl <- delta_list[[ fl ]]
            delta.norm  <- delta_norms[[ fl ]]

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
              trial.unmix <- unmix_fun( cell.raw, cell.spectra.curr, cell.weights )

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

        cell.unmixed <- unmix_fun( cell.raw, cell.spectra.final, cell.weights )

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
