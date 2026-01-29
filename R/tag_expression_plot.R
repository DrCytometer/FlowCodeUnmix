# tag_expression_plot.r

#' @title Plot Tag Expression
#'
#' @description
#' Creates ridgeline (density) plots of the FlowCode fluorophore expression
#' levels from unmixed data.
#'
#' @importFrom ggplot2 ggplot aes facet_wrap scale_x_continuous theme_classic
#' @importFrom ggplot2 labs ggsave
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom ggridges geom_density_ridges
#' @importFrom dplyr mutate group_by sample_n ungroup filter n
#'
#' @param unmixed Matrix of unmixed flow cytometry expression data. Cells in
#' rows and fluorophores in columns. Fluorophores must include `flowcode.fluors`.
#' @param flowcode.fluors Named character vector of fluorophores used to identify
#' the FlowCodes. Names should be FlowCode epitope tags, values should be the
#' corresponding fluorophores.
#' @param flowcode.ids Numeric vector containing the FlowCode combination
#' assignments per cell (result of `debarcode()`).
#' @param combo.df Data frame describing the valid combinations of FlowCodes.
#' Structure: One row per combination. Columns are `Id`, `Procode.tag1`,
#' `Procode.tag2` and `Procode.tag3`, describing the name (e.g., CRISPR target),
#' and three epitopes for the combination, respectively.
#' @param asp The AutoSpectral parameter list.
#' @param output.dir File path determining where the plot will be saved. Default
#' is `./flowcode_spectra`.
#' @param title Character, name to add to the saved plots. Default is `Uncorrected`.
#' @param max.per.group Numeric, maximum number of cells to plot per group
#' (FlowCode combination). Larger numbers will take longer. Default is `2000`.
#' @param height Numeric, default `20`. Width (in inches) of the saved plot.
#' @param width Numeric, default `20`. Height (in inches) of the saved plot.
#'
#' @return None. Saves the plot to the specified file path in `output.dir`.
#'
#' @export

tag.expression.plot <- function(
    unmixed,
    flowcode.fluors,
    flowcode.ids,
    combo.df,
    asp,
    output.dir = "./flowcode_spectra",
    title = "Uncorrected",
    max.per.group = 2000,
    height = 20,
    width = 20
) {

  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

  # set title & filename
  title <- paste( title, "debarcoded_tag_expression", sep = "_" )
  filename <- paste0( title, ".jpg" )

  # axes and transforms
  x.breaks <- asp$ribbon.breaks[ asp$ribbon.breaks < asp$expr.data.max ]
  x.axis.labels <- sapply( x.breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0("10^", log10( abs( x ) ) ) )
  } )
  x.limits <- c(
    asp$default.transformation.param$width * 5,
    asp$expr.data.max
  )

  # FlowJo-like biexp transform
  biexp.transform <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = asp$default.transformation.param$max.range,
    pos          = asp$default.transformation.param$pos,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = asp$default.transformation.param$width,
    inverse      = FALSE
  )

  # build long-form dataset
  plot.list <- lapply( names( flowcode.fluors ), function( tag ) {
    channel <- flowcode.fluors[[ tag ]]
    stopifnot( channel %in% colnames( unmixed ) )

    expr <- unmixed[ , channel ]

    data.frame(
      Fluor      = tag,
      Expression = expr,
      GroupID    = flowcode.ids,
      Label      = combo.df$Id[ flowcode.ids ],
      stringsAsFactors = FALSE
    )
  } )

  plot.df <- do.call( rbind, plot.list )

  # remove any missing info
  plot.df <- plot.df |>
    filter( GroupID > 0, !is.na( Label ) )

  # downsample if needed
  plot.df <- plot.df |>
    group_by( Fluor, Label ) |>
    sample_n( size = min( n(), max.per.group ) ) |>
    ungroup()

  plot.df <- plot.df |>
    mutate( Expr.trans = biexp.transform( Expression ) )

  p <- suppressMessages(
    ggplot( plot.df, aes( x = Expr.trans, y = Label, group = Label ) ) +
      geom_density_ridges(
        rel_min_height = 0.01,
        fill = "steelblue",
        alpha = 0.45,
        scale = 2,
        na.rm = TRUE
      ) +
      facet_wrap( ~ Fluor, scales = "free_x" ) +
      scale_x_continuous(
        name   = "Expression",
        breaks = biexp.transform( x.breaks ),
        labels = x.axis.labels,
        limits = biexp.transform( x.limits )
      ) +
      theme_classic( base_size = 14 ) +
      labs(
        title = title,
        y = "Combo ID (Gene Label)"
      )
  )

  suppressMessages(
    ggsave(
      filename = file.path( output.dir, filename ),
      plot     = p,
      width    = width,
      height   = height
    )
  )
}


