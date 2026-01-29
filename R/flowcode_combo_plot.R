# flowcode_combo_plot.r

#' @title FlowCode Combo Plot
#'
#' @description
#' Creates side-by-side pseudocolor density dot plots (standard flow cytometry)
#' for the FlowCode markers. Plots the unmixed data on a biexponential scale
#' using `AutoSpectral::create.biplot()`.
#'
#' @importFrom AutoSpectral create.biplot
#' @importFrom ggplot2 ggsave
#' @importFrom cowplot plot_grid
#'
#' @param unmixed Matrix of unmixed flow cytometry expression data. Cells in
#' rows and fluorophores in columns. Fluorophores must include `flowcode.fluors`.
#' @param corrected Corrected matrix of unmixed flow cytometry expression data
#' after processing for FRET errors.
#' @param flowcode.fluors Named character vector of fluorophores used to identify
#' the FlowCodes. Names should be FlowCode epitope tags, values should be the
#' corresponding fluorophores.
#' @param asp The AutoSpectral parameter list. Prepare using
#' `get.autospectral.param`.
#' @param output.dir File path determining where the plot will be saved. Default
#' is `./flowcode_spectra`.
#' @param filename Character, name of the saved plot. Default is
#' `Corrected_biplots.jpg`.
#' @param width Numeric, default `10`. Width (in inches) of the saved plot.
#' @param height Numeric, default `5`. Height (in inches) of the saved plot.
#'
#' @return None. Saves the plot to the specified file path in `output.dir`.
#'
#' @export

flowcode.combo.plot <- function(
    unmixed,
    corrected,
    flowcode.fluors,
    asp,
    output.dir = "./flowcode_spectra",
    filename = "Corrected_biplots.jpg",
    width = 5,
    height = 10
  ) {

  if ( !requireNamespace( "AutoSpectral", quietly = TRUE ) ) {
    stop(
      "The AutoSpectral package is required but is not installed or not available.",
      call. = FALSE
    )
  }

  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

  n.fluors <- length( flowcode.fluors )

  if ( n.fluors %% 2 != 0 ) {
    flowcode.fluors <- c( flowcode.fluors, flowcode.fluors[ 1 ] )
    n.fluors <- length( flowcode.fluors )
  }

  all.plots <- list()

  pair.index <- 1

  for ( i in seq( 1, n.fluors, by = 2 ) ) {

    x.dim <- flowcode.fluors[ i ]
    y.dim <- flowcode.fluors[ i + 1 ]

    title.orig <- sprintf( "%s vs %s : original", x.dim, y.dim )
    title.corr <- sprintf( "%s vs %s : corrected", x.dim, y.dim )

    p1 <- AutoSpectral::create.biplot(
      plot.data = unmixed,
      x.dim = x.dim,
      y.dim = y.dim,
      x.min = -5000, x.max = asp$expr.data.max,
      y.min = -5000, y.max = asp$expr.data.max,
      x.width.basis = -1000,
      y.width.basis = -1000,
      max.points = 3e5,
      asp = asp,
      title = title.orig,
      save = FALSE
    )

    p2 <- AutoSpectral::create.biplot(
      plot.data = corrected,
      x.dim = x.dim,
      y.dim = y.dim,
      x.min = -5000, x.max = asp$expr.data.max,
      y.min = -5000, y.max = asp$expr.data.max,
      x.width.basis = -1000,
      y.width.basis = -1000,
      max.points = 3e5,
      asp = asp,
      title = title.corr,
      save = FALSE
    )

    # Store pair as two consecutive entries
    all.plots[[ pair.index     ]] <- p1
    all.plots[[ pair.index + 1 ]] <- p2

    pair.index <- pair.index + 2
  }

  # Number of rows = number of pairs
  n.pairs <- n.fluors / 2

  combined <- cowplot::plot_grid(
    plotlist = all.plots,
    ncol = 2,
    nrow = n.pairs,
    labels = rep( c( "original", "corrected" ), n.pairs ),
    label_size = 12,
    label_x = 0.01,
    label_y = 0.99,
    hjust = 0,
    vjust = 1
  )

  ggsave(
    filename = file.path( output.dir, filename ),
    plot = combined,
    width = width,
    height = height * n.pairs
  )
}

