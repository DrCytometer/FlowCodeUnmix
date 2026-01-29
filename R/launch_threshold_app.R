# launch_threshold_app.r

#' @title Launch the Thresholding Shiny App
#'
#' @description
#' A simple wrapper function to launch the manual threshold selection app.
#'
#' @importFrom shiny runApp
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
#'
#' @param flowcode.combo.file Optional, default is `NULL`. Provide this when you
#' wish to establish new thresholds using the app and will be working with an
#' FCS file rather than loading the RDS object. When given: File name and path
#' to the CSV file containing the information describing your FlowCode library.
#' Describes the valid combinations of FlowCodes. Structure: One row per
#' combination. Columns are `Id`, `Procode.tag1`, `Procode.tag2` and
#' `Procode.tag3`, describing the name (e.g., CRISPR target), and three epitopes
#' for the combination, respectively.
#' @param flow.control Optional, default is `NULL`. Provide this when you
#' wish to establish new thresholds using the app and will be working with an
#' FCS file rather than loading the RDS object. A list containing flow cytometry
#' control parameters.
#' @param spectra Optional, default is `NULL`. Provide this when you wish to
#' establish new thresholds using the app and will be working with an FCS file
#' rather than loading the RDS object.Spectral signatures of fluorophores,
#' normalized between 0 and 1, with fluorophores in rows and detectors in columns.
#'
#' @export

launch.threshold.app <- function(
    flowcode.combo.file = NULL,
    flow.control = NULL,
    spectra = NULL
  ) {

  if ( !requireNamespace( "AutoSpectral", quietly = TRUE ) ) {
    stop(
      "The AutoSpectral package is required but is not installed or not available.",
      call. = FALSE
    )
  }

  if ( !( is.null( flowcode.combo.file ) & is.null( flow.control ) ) ) {
    # read in combo file
    combo.df <- utils::read.csv( flowcode.combo.file )
    flowcode.tags <- unique( unlist( combo.df[ , -1 ] ) )

    # define flowcode channel-tag correspondence, case-independent
    # check against flow.control
    tag.lookup <- toupper( flowcode.tags )
    antigen.lookup <- toupper( flow.control$antigen )

    flowcode.fluors <- flow.control$fluorophore[ match( tag.lookup, antigen.lookup ) ]

    names( flowcode.fluors ) <- flowcode.tags
  }

  # store flowcode fluors and spectra for access in the app
  options(
    flowcode_fluors_list = flowcode.fluors %||% NULL,
    Spectra = spectra,
    threshold_app_wd = getwd()
  )

  # check for the app
  app.dir <- system.file( "shiny", "ThresholdApp", package = "FlowCodeUnmix" )

  if ( app.dir == "" )
    stop( "App not found. Make sure the package `FlowCodeUnmix` installed correctly." )

  # launch the app
  shiny::runApp( appDir = app.dir, display.mode = "normal" )
}
