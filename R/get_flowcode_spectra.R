# get_flowcode_spectra.r

#' @title Get FlowCode Spectra
#'
#' @description
#' Measure spectral unmixing errors per barcode in a sample stained only with
#' FlowCode tag antibodies (`flowcode.backbone.fcs`). Measures FRET error and
#' returns the calculated correction matrices and thresholds for debarcoding.
#'
#' @importFrom flowCore read.FCS exprs
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
#' @param plot.corrections Logical, default is `FALSE`. If `TRUE`, the FRET
#' unmixing will be applied to `flowcode.backbone.fcs` and plots of the
#' corrected data (FlowCode channels only) will be produced in `output.dir`.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress
#' messaging.
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
    verbose = TRUE
  ) {

  if ( !requireNamespace( "AutoSpectral", quietly = TRUE ) ) {
    stop(
      "The AutoSpectral package is required but is not installed or not available.",
      call. = FALSE
    )
  }

  # set up
  if ( verbose ) message( "Reading thresholds file" )
  flowcode.thresholds.df <- utils::read.csv( thresholds.file )
  flowcode.thresholds <- flowcode.thresholds.df$Threshold_Raw
  names( flowcode.thresholds ) <- flowcode.thresholds.df$Fluor

  # read in the data and run checks
  if ( verbose ) message( "Loading RDS file" )
  backbone <- readRDS( backbone.rds )

  expected <- list(
    Unmixed = list(
      class = c("matrix"),
      nrow_min = 100
    ),
    Raw = list(
      class = c("matrix"),
      nrow_min = 100
    ),
    Flowcode.fluors = list(
      class = c("character"),
      length_min = 4
    ),
    Spectra = list(
      class = c("matrix"),
      nrow_min = 2
    ),
    Combos = list(
      class = c("data.frame"),
      nrow_min = 10
    )
  )

  # check that the backbone file contains the right stuff
  check.rds.contents( backbone, required = expected )

  # if it didn't fail, we're good to go
  if ( verbose ) message( "RDS file has passed checks" )

  # unpack small stuff
  spectra <- backbone$Spectra
  combo.df <- backbone$Combos
  flowcode.fluors <- backbone$Flowcode.fluors

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  fluorophores <- rownames( spectra )
  fluorophore.n <- length( fluorophores )
  spectral.channel <- colnames( spectra )
  detector.n <- length( spectral.channel )
  combo.n <- nrow( combo.df )

  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

  weights <- abs( colMeans( backbone$Raw ) )
  weights[ weights < 1e-6 ] <- 1e-6
  weights <- 1 / weights

  # debarcode backbone
  if ( verbose ) {
    message( paste0( "\033[33m", "Debarcoding the backbone", "\033[0m" ) )
  }

  flowcode.ids <- debarcode(
    unmixed = backbone$Unmixed,
    flowcode.fluors = flowcode.fluors,
    flowcode.thresholds = flowcode.thresholds,
    combo.df = combo.df
  )

  if ( figures ) {
    # QC plot showing density histograms of flowcode tags by combo
    valid.combos <- which( flowcode.ids %in% seq_len( combo.n ) )

    tag.expression.plot(
      unmixed = backbone$Unmixed[ valid.combos, ],
      flowcode.fluors = flowcode.fluors,
      flowcode.ids = flowcode.ids[ valid.combos ],
      combo.df = combo.df,
      output.dir = output.dir,
      asp = asp
    )
  }


  ### Measuring the per-combination error for each FlowCode
  if ( verbose ) message( paste0( "\033[33m", "Measuring background", "\033[0m" ) )

  # which cells lack any FlowCode expression?
  untransduced.idx <- which( flowcode.ids == 0 )

  # get a sample of 500 of these cells for downstream processing of background
  set.seed( 42 )
  background.idx <- sample( untransduced.idx, 500 )

  if ( verbose ) {
    message(
      paste0(
        "\033[33m",
        "Calculating FRET errors per FlowCode combination",
        "\033[0m"
      )
    )
  }

  # set up collection
  fret.spectra <- list()

  # loop through each combo
  # measuring the spillover/unmixing error in non-combo channels
  # this can be parallelized using parLapply/mclapply later
  for ( id in seq_len( combo.n ) ) {
    # get cells with this combo
    combo.idx <- which( flowcode.ids == id )
    n.cells <- length( combo.idx )

    # only proceed with FRET assessment if we have at least three cells
    if ( n.cells > 3 ) {
      # define combo and non-combo fluors
      combo.tags <- unlist( combo.df[  id, -1 ] )
      combo.fluors <- flowcode.fluors[
        match( combo.tags, names( flowcode.fluors ) ) ]
      non.combo.fluors <- rownames( spectra )[
        !( rownames( spectra ) %in% combo.fluors ) ]



      ### Option 1: directly subtract untransduced raw from combo+ cells
      ## this can be done outside the loop in a single step
      #differential.raw <- sweep( backbone$Raw[ idx, ], 2, untransduced.median, "-" )
      #diff.unmix <- unmix.wls( differential.raw, spectra, weights )
      #resid <- differential.raw - diff.unmix %*% spectra
      #non.combo.fitted <- diff.unmix[ , non.combo.fluors, drop = FALSE ] %*%
      #  spectra[ non.combo.fluors, , drop = FALSE ]

      ### Option 2: unmix both positive and background
      combo.unmixed <- AutoSpectral::unmix.wls(
        backbone$Raw[ combo.idx, ],
        spectra,
        weights
      )
      background.unmixed <- AutoSpectral::unmix.wls(
        backbone$Raw[ background.idx, ],
        spectra,
        weights
      )

      # select brightest events only in the FlowCode channels for this combo
      if ( n.cells > 500 ) {
        brightest.events <- order(
          rowSums(
            combo.unmixed[ , combo.fluors, drop = FALSE ]
          ),
          decreasing = TRUE
        )[ 1:500 ]
        combo.idx <- combo.idx[ brightest.events ]

        combo.unmixed <- combo.unmixed[ combo.idx, , drop = FALSE ]
        n.cells <- 500
      }

      # project back non-combo fluors
      combo.non.combo.fitted <- combo.unmixed[ , non.combo.fluors, drop = FALSE ] %*%
        spectra[ non.combo.fluors, , drop = FALSE ]
      background.non.combo.fitted <- background.unmixed[ , non.combo.fluors, drop = FALSE ] %*%
        spectra[ non.combo.fluors, , drop = FALSE ]

      # get the median background
      median.background <- apply(
        background.non.combo.fitted,
        2,
        stats::median
      )

      # subtract background
      non.combo.fitted <- sweep(
        combo.non.combo.fitted,
        2,
        median.background,
        FUN = "-"
      )

      # calculate residuals for combo+ cells
      resid <- backbone$Raw[ combo.idx, ] - combo.unmixed %*% spectra

      # scale som.dim to n.cells, with limits between 2 and 10
      som.dim <- if ( n.cells >= 500 ) {
        10
      } else {
        min( 10, max( ceiling( sqrt( n.cells ) ) - 1, 2 ) )
      }

      # cluster on incorrect data + residuals to get variation
      input.data <- resid + non.combo.fitted
      set.seed( 42 )
      map <- FlowSOM::SOM(
        input.data,
        xdim = som.dim,
        ydim = som.dim,
        silent = TRUE
      )

      # L-Inf (peak) normalize to get FRET spectra
      combo.spectra <- t(
        apply( map$codes[ , spectral.channel ], 1, function( x ) {
          max.x <- ifelse( max( abs( x ) ) > max( x ), min( x ), max( x ) )
          x / max.x } )
      )
      combo.spectra <- as.matrix( stats::na.omit( combo.spectra ) )

      # get the median of the incorrect data to use as the starting point for correction
      median.fret <- apply(
        resid + non.combo.fitted,
        2,
        stats::median
      )

      # L-Inf (peak) normalize
      median.fret.spectrum <- median.fret / max( abs( median.fret ) )

      # add combo median as row 1 (first FRET variant spectrum)
      combo.spectra <- rbind( median.fret.spectrum, combo.spectra )

      # add rownames for tracking if desired (not used, but doesn't hurt)
      rownames( combo.spectra ) <- paste0( id, 1:nrow( combo.spectra ) )

      if ( figures ) {
        combo.label <- combo.df$Id[ id ]
        if ( verbose ) {
          message(
            paste0(
              "\033[33m",
              "Plotting FRET variation for ",
              combo.label,
              "\033[0m"
            )
          )
        }

        plot.title <- paste(
          combo.fluors,
          collapse = "_"
        )
        plot.title <- paste(
          "FRET_errors_for",
          combo.label,
          plot.title,
          sep = "_"
        )

        # plot what the FRET variation looks like--very messy
        AutoSpectral::spectral.variant.plot(
          spectra.variants = combo.spectra,
          median.spectrum = combo.spectra[ 1, ], # median spectrum
          title = plot.title,
          save = TRUE,
          plot.dir = output.dir
        )
      }

      fret.spectra[[ id ]] <- combo.spectra

    } else {
      # warn and store a placeholder zero spectrum to preserve structure
      combo.label <- combo.df$Id[ id ]
      warning(
        paste(
          "Not enough cells (fewer than 3) were present for",
          combo.label, ":", combo.fluors
        )
      )
      fret.spectra[[ id ]] <- matrix( 0, nrow = 1, ncol = detector.n )
    }
  }

  names( fret.spectra ) <- combo.df$Id


  ### calculate additional data required for unmixing

  # calculate delta (differences between FRET variants) for each combo
  fret.delta.list <- lapply( names( fret.spectra ), function( fc ) {
    fret.spectra[[ fc ]] - matrix(
      fret.spectra[[ fc ]][ 1, ], # median spectrum
      nrow = nrow( fret.spectra[[ fc ]] ),
      ncol = detector.n,
      byrow = TRUE
    )
  } )
  names( fret.delta.list ) <- names( fret.spectra )

  # calculate delta norms (penalty distance for each variant) for each combo
  fret.delta.norms <- lapply( fret.delta.list, function( d ) {
    sqrt( rowSums( d^2 ) )
  } )

  # assign logical matrix describing which fluors are present and absent in each combo
  flowcode.combo.logical <- matrix(
    0L,
    nrow = combo.n,
    ncol = length( flowcode.fluors ),
    dimnames = list( combo.df$Id, flowcode.fluors )
  )

  # columns in combo.df that contain fluorophore tags
  tag.cols <- c( "Procode.tag1", "Procode.tag2", "Procode.tag3" )

  # fill with 1s where a fluorophore is present on that combo
  for ( i in seq_len( combo.n ) ) {
    tags <- unlist( combo.df[ i, tag.cols ], use.names = FALSE )
    tags <- flowcode.fluors[ tags ]
    flowcode.combo.logical[ i, tags ] <- 1L
  }

  ### Plotting and corrections--this needs to be updated following unmix.flowcode()
  # calculate and plot corrections for the backbone, if desired
  if ( figures & plot.corrections ) {
    if ( verbose ) {
      message(
        paste0(
          "\033[33m",
          "Calculating per-cell FRET corrections",
          "\033[0m"
        )
      )
    }

    corrected <- backbone$Unmixed

    # reunmix backbone, plot before/after for flowcode tags
    for ( id in seq_len( combo.n ) ) {

      # get cells with this combo
      idx <- which( flowcode.ids == id )

      if ( length( idx ) == 0 )
        next

      # get FRET spectra for this combo
      combo.spectra <- fret.spectra[[ id ]]

      # check for all zero matrix, if so, next
      if ( nrow( combo.spectra ) == 1 && all( combo.spectra == 0 ) )
        next

      # this needs to be a per-cell loop internally here
      cell.raw <- backbone$Raw[ idx, , drop = FALSE ]
      weights <- abs( colMeans( cell.raw ) )
      weights[ weights < 1e-6 ] <- 1e-6
      weights <- 1 / weights

      combo.tags <- unlist( combo.df[  id, -1 ] )
      combo.fluors <- flowcode.fluors[ match( combo.tags, names( flowcode.fluors ) ) ]

      cell.spectra.curr <- spectra[ combo.fluors, ]
      cell.unmixed <- AutoSpectral::unmix.wls(
        cell.raw,
        cell.spectra.curr,
        weights
      )
      resid  <- cell.raw - cell.unmixed %*% cell.spectra.curr
      error.final <- rowSums( abs( resid ) )
      best.fret <- rep( 0, length( idx ) )

      # unmix all these cells with each fret combo plus combo fluors
      # determine best FRET match
      ## change to scoring
      for ( fret in seq_len( nrow( combo.spectra ) ) ) {

        combined.spectra <- rbind( cell.spectra.curr, combo.spectra[ fret, ] )
        unmixed.fret <- AutoSpectral::unmix.wls(
          cell.raw,
          combined.spectra,
          weights
        )
        error.fret <- rowSums( abs( cell.raw - ( unmixed.fret %*% combined.spectra ) ) )
        improved <- which( error.fret < error.final )

        if ( length( improved ) != 0 ) {
          error.final[ improved ] <- error.fret[ improved ]
          best.fret[ improved ] <- fret
        }
      }
      # re-unmix with full fluorophore spectra + FRET
      for ( fret in unique( best.fret ) ) {
        combined.spectra <- rbind( spectra, combo.spectra[ fret, ] )
        fret.idx <- which( best.fret == fret )
        new.unmix <- AutoSpectral::unmix.wls(
          cell.raw[ fret.idx, , drop = FALSE ],
          combined.spectra,
          weights
        )
        corrected[ idx[ fret.idx ], fluorophores ] <- new.unmix[ , fluorophores ]
      }
    }

    if ( verbose ) {
      message(
        paste0(
          "\033[33m",
          "Plotting FlowCode data with and without FRET correction",
          "\033[0m"
        )
      )
    }

    # plot the corrections as ridgeline (density) plots
    tag.expression.plot(
      unmixed = corrected[ valid.combos, ],
      flowcode.fluors = flowcode.fluors,
      flowcode.ids = flowcode.ids[ valid.combos ],
      combo.df = combo.df,
      output.dir = output.dir,
      asp = asp,
      title = "Corrected"
    )

    # create side-by-side pseudocolor biplots
    flowcode.combo.plot(
      unmixed = backbone$Unmixed,
      corrected = corrected,
      flowcode.fluors = flowcode.fluors,
      asp = asp,
      output.dir = output.dir
    )
  }

  if ( verbose )
    message( paste0( "\033[33m", "Spectral variation computed!", "\033[0m" ) )

  # construct list of items to return
  flowcode.spectra <- list(
    Thresholds = flowcode.thresholds,
    FRET = fret.spectra,
    Flowcode.fluors = flowcode.fluors,
    Combos = combo.df,
    Delta = fret.delta.list,
    Delta.norms = fret.delta.norms,
    Logical.combo = flowcode.combo.logical
  )

  # save it for later in case subsequent processing is needed
  saveRDS( flowcode.spectra, file = file.path( output.dir, filename ) )

  return( flowcode.spectra )
}
