# launch_threshold_app.r

#' @title Launch the Thresholding Shiny App
#'
#' @description
#' A simple wrapper function to launch the manual threshold selection app.
#'
#' @importFrom shiny runApp shinyOptions
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarLayout sidebarPanel
#' @importFrom shiny mainPanel fileInput selectInput actionButton
#' @importFrom shiny plotOutput tableOutput verbatimTextOutput
#' @importFrom shiny reactiveValues observeEvent renderPlot renderTable
#' @importFrom shiny updateSelectInput req showNotification brushOpts
#' @importFrom shiny tags
#' @importFrom shinyFiles shinyDirButton shinyDirChoose parseDirPath
#' @importFrom ggplot2 ggplot aes geom_vline theme_classic theme labs
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous expansion
#' @importFrom scattermore geom_scattermore
#' @importFrom dplyr filter
#' @importFrom AutoSpectral create.biplot
#' @importFrom sp point.in.polygon
#'
#' @param backbone.rds Optional, default is `NULL`. Path to a backbone RDS
#' file produced by `unmix.backbone()` (e.g. `"FlowCode_Backbone.rds"`, or its
#' downsampled `"Small_..."` counterpart). When given, the app loads this file
#' automatically on startup instead of requiring you to browse for it via the
#' file picker.
#' @param flowcode.fluors Optional, default is `NULL`. Named character vector
#' of fluorophores used to identify FlowCode epitope tags (names = tags,
#' values = fluorophore/detector names). If supplied, takes precedence over
#' whatever is embedded in `backbone.rds` or derived from
#' `flowcode.combo.file`/`flow.control`. Provide this directly, OR provide
#' `flowcode.combo.file` + `flow.control`, OR let the app read it from
#' `backbone.rds` — one of the three is needed to see fluorophore labels when
#' working from a raw FCS file.
#' @param spectra Optional, default is `NULL`. Spectral signatures of
#' fluorophores, normalized between 0 and 1, with fluorophores in rows and
#' detectors in columns. If supplied, takes precedence over whatever is
#' embedded in `backbone.rds`.
#' @param asp Optional, default is `NULL`. The AutoSpectral parameter list
#' (see `get.autospectral.param()`), used to seed default transform values
#' and plot styling. If omitted, the app falls back to built-in defaults.
#' @param flowcode.combo.file Optional, default is `NULL`. File name and path
#' to the CSV file describing your FlowCode library, used together with
#' `flow.control` to derive `flowcode.fluors` when `flowcode.fluors` isn't
#' supplied directly. Structure: one row per combination, columns `Id`,
#' `Procode.tag1`, `Procode.tag2`, `Procode.tag3`.
#' @param flow.control Optional, default is `NULL`. A list containing flow
#' cytometry control parameters, used together with `flowcode.combo.file`.
#' @param output.dir Optional, default is `NULL`. Initial output folder for
#' saved thresholds/plots. Defaults to `./flowcode_spectra` under the current
#' working directory if that exists, otherwise the working directory itself.
#'
#' @export

launch.threshold.app <- function(
    backbone.rds = NULL,
    flowcode.fluors = NULL,
    spectra = NULL,
    asp = NULL,
    flowcode.combo.file = NULL,
    flow.control = NULL,
    output.dir = NULL
) {

  if ( !requireNamespace( "AutoSpectral", quietly = TRUE ) ) {
    stop(
      "The AutoSpectral package is required but is not installed or not available.",
      call. = FALSE
    )
  }

  # `flowcode.combo.file` and `flow.control` must be given together, or not
  # at all — a partial pair used to fail confusingly downstream instead of
  # here.
  combo.given   <- !is.null( flowcode.combo.file )
  control.given <- !is.null( flow.control )
  if ( xor( combo.given, control.given ) ) {
    stop(
      "Provide both `flowcode.combo.file` and `flow.control` together, or neither.",
      call. = FALSE
    )
  }

  if ( is.null( flowcode.fluors ) && combo.given && control.given ) {
    combo.df <- utils::read.csv( flowcode.combo.file )
    flowcode.tags <- unique( unlist( combo.df[ , -1 ] ) )

    # define flowcode channel-tag correspondence, case-independent
    # check against flow.control
    tag.lookup <- toupper( flowcode.tags )
    antigen.lookup <- toupper( flow.control$antigen )

    flowcode.fluors <- flow.control$fluorophore[ match( tag.lookup, antigen.lookup ) ]
    names( flowcode.fluors ) <- flowcode.tags
  }

  if ( !is.null( backbone.rds ) && !file.exists( backbone.rds ) )
    stop( "`backbone.rds` not found: ", backbone.rds, call. = FALSE )

  if ( is.null( output.dir ) ) {
    candidate  <- file.path( getwd(), "flowcode_spectra" )
    output.dir <- if ( dir.exists( candidate ) ) candidate else getwd()
  }

  # Pass data to the app via shinyOptions() — scoped to shiny (unlike base
  # options(), which would leak into the wider R session), and readable
  # inside the app without ever reaching into .GlobalEnv.
  shiny::shinyOptions(
    flowcode_fluors_list = flowcode.fluors,
    spectra              = spectra,
    asp                  = asp,
    backbone_rds         = backbone.rds,
    threshold_app_wd     = getwd(),
    threshold_app_outdir = output.dir
  )

  # check for the app
  app.dir <- system.file( "shiny", "ThresholdApp", package = "FlowCodeUnmix" )

  if ( app.dir == "" )
    stop( "App not found. Make sure the package `FlowCodeUnmix` installed correctly." )

  # launch the app
  shiny::runApp( appDir = app.dir, display.mode = "normal" )
}
