# app.R  —  FlowCodeUnmix Interactive Thresholding App

# Allow large uploads (e.g., 2 GB)
options(shiny.maxRequestSize = 2 * 1024^3)

# ==============================================================================
# Libraries
# ==============================================================================
library(shiny)
library(bslib)
library(shinyjs)
library(shinyFiles)
library(ggplot2)
library(dplyr)

# ==============================================================================
# Helpers
# ==============================================================================

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# Build a biexponential transform function using flowWorkspace
make_biexp <- function(length, max.range, pos, neg, width, inverse = FALSE) {
  tryCatch(
    flowWorkspace::flowjo_biexp(
      channelRange = length,
      maxValue     = max.range,
      pos          = pos,
      neg          = neg,
      widthBasis   = width,
      inverse      = inverse
    ),
    error = function(e) NULL
  )
}

# Apply biexp transform to a numeric vector, returning transformed values.
# Returns NULL on failure.
apply_biexp <- function(vals, length, max.range, pos, neg, width) {
  xf <- make_biexp(length, max.range, pos, neg, width, inverse = FALSE)
  if (is.null(xf)) return(NULL)
  xf(vals)
}

# Check whether a point is inside a polygon
points_in_polygon <- function(px, py, poly_x, poly_y) {
  sp::point.in.polygon(px, py, poly_x, poly_y) > 0
}

# ==============================================================================
# UI
# ==============================================================================

ui <- page_navbar(
  title  = tags$span(
    tags$img(src = "https://www.r-project.org/logo/Rlogo.svg",
             height = "22px", style = "margin-right:8px; vertical-align:middle;"),
    "FlowCodeUnmix — Thresholding"
  ),
  theme  = bs_theme(
    bootswatch = "flatly",
    primary    = "#2C6FAC",
    base_font  = font_google("Inter")
  ),
  header = useShinyjs(),

  # ============================================================
  # TAB 1 — Data & Gate
  # ============================================================
  nav_panel(
    "Gate",
    icon = icon("crop"),

    layout_sidebar(
      sidebar = sidebar(
        width = 300,

        # ---- File loading ----
        card(
          card_header(icon("folder-open"), " Load Data"),
          card_body(
            shinyFilesButton(
              "file", "Load FCS / RDS",
              "Please select a file",
              multiple = FALSE,
              class    = "btn-primary btn-sm w-100"
            ),
            br(), br(),
            verbatimTextOutput("file_status", placeholder = TRUE)
          )
        ),

        # ---- Gate mode ----
        card(
          card_header(icon("draw-polygon"), " Gate Mode"),
          card_body(
            radioButtons(
              "gate_mode", NULL,
              choices  = c("Rectangular brush" = "rect",
                           "Free-form polygon" = "poly"),
              selected = "rect",
              inline   = FALSE
            ),
            conditionalPanel(
              condition = "input.gate_mode == 'poly'",
              helpText(
                icon("info-circle"),
                "Click to add vertices. Double-click the plot to close & apply gate."
              ),
              actionButton("clear_polygon", "Clear polygon",
                           class = "btn-outline-danger btn-sm w-100",
                           icon  = icon("trash"))
            ),
            hr(),
            actionButton("refresh_gate_plot", "Refresh gate plot",
                         class = "btn-outline-secondary btn-sm w-100",
                         icon  = icon("sync"))
          )
        ),

        # ---- Gate stats ----
        card(
          card_header(icon("filter"), " Gated cells"),
          card_body(
            verbatimTextOutput("gate_stats"),
            uiOutput("gate_msg_ui")
          )
        )
      ),

      # ---- Main area ----
      layout_columns(
        col_widths = c(12),
        card(
          card_header("FSC-A vs SSC-A  — draw a gate"),
          card_body(
            plotOutput(
              "gate_plot",
              height   = "420px",
              brush    = brushOpts(id = "gate_brush", resetOnNew = FALSE),
              click    = "gate_click",
              dblclick = "gate_dblclick"
            )
          )
        )
      )
    )
  ),

  # ============================================================
  # TAB 2 — Transform & Threshold
  # ============================================================
  nav_panel(
    "Threshold",
    icon = icon("sliders"),

    layout_sidebar(
      sidebar = sidebar(
        width = 300,

        # ---- Channel selector ----
        card(
          card_header("Channels"),
          card_body(
            tags$b("X-axis (threshold) channel:"),
            div(
              id    = "fluor_container",
              style = paste(
                "max-height:240px; overflow-y:auto;",
                "background:#f8f9fa; padding:8px;",
                "border:1px solid #dee2e6; border-radius:6px;"
              ),
              uiOutput("fluor_selector_ui")
            ),
            hr(),
            selectInput("y_axis_select", "Y-axis channel:", choices = NULL)
          )
        ),

        # ---- Transform tuning ----
        card(
          card_header(icon("wand-magic-sparkles"), " Transform"),
          card_body(
            checkboxInput("link_axes", "Link X & Y transforms", value = TRUE),
            hr(),
            h6("X-axis"),
            sliderInput("xform_width_x", "Width basis",
                        min = -1000, max = 0, value = -1000, step = 25),
            numericInput("xform_pos_x", "Positive log decades",
                         value = 4.418, step = 0.1, min = 1, max = 8),
            numericInput("xform_neg_x", "Negative log decades",
                         value = 0, step = 0.1, min = 0, max = 4),
            conditionalPanel(
              condition = "!input.link_axes",
              hr(),
              h6("Y-axis"),
              sliderInput("xform_width_y", "Width basis",
                          min = -1000, max = 0, value = -1000, step = 25),
              numericInput("xform_pos_y", "Positive log decades",
                           value = 4.418, step = 0.1, min = 1, max = 8),
              numericInput("xform_neg_y", "Negative log decades",
                           value = 0, step = 0.1, min = 0, max = 4)
            ),
            hr(),
            actionButton("refresh_threshold_plot",
                         "Refresh plot",
                         class = "btn-primary btn-sm w-100",
                         icon  = icon("sync")),
            helpText(
              icon("info-circle"),
              "Click the plot to place a threshold line."
            )
          )
        ),

        # ---- Thresholds table ----
        card(
          card_header(icon("table"), " Thresholds set"),
          card_body(tableOutput("threshold_table"))
        )
      ),

      # ---- Main area ----
      card(
        card_header("Biexponential scatter plot — click to set threshold"),
        card_body(
          plotOutput(
            "threshold_plot",
            height = "480px",
            click  = "threshold_click"
          )
        )
      )
    )
  ),

  # ============================================================
  # TAB 3 — Export
  # ============================================================
  nav_panel(
    "Export",
    icon = icon("download"),

    layout_columns(
      col_widths = c(6, 6),

      card(
        card_header(icon("floppy-disk"), " Save Thresholds"),
        card_body(
          shinyDirButton(
            "output_dir", "Change output folder",
            "Select output folder",
            class = "btn-outline-secondary btn-sm"
          ),
          br(), br(),
          verbatimTextOutput("output_dir_display"),
          hr(),
          textInput("out_filename", "Output CSV filename:", value = "thresholds.csv"),
          actionButton("save_thr", "Save thresholds & plots",
                       class = "btn-success",
                       icon  = icon("floppy-disk"))
        )
      ),

      card(
        card_header(icon("circle-info"), " Help"),
        card_body(
          tags$ul(
            tags$li("Load an FCS or RDS file on the ", tags$b("Gate"), " tab."),
            tags$li("Draw a gate using brush or polygon, then switch to ", tags$b("Threshold"), "."),
            tags$li("Tune the biexponential transform until populations separate cleanly."),
            tags$li("Click the scatter plot to place a threshold line."),
            tags$li("Save thresholds and per-fluorophore plots on the ", tags$b("Export"), " tab.")
          )
        )
      )
    )
  )
)

# ==============================================================================
# Server
# ==============================================================================

server <- function(input, output, session) {

  # --------------------------------------------------------------------------
  # Reactive state
  # --------------------------------------------------------------------------

  rv <- shiny::reactiveValues(
    data            = NULL,
    gated_data      = NULL,
    thresholds      = list(),
    current_fluor   = NULL,
    all_fluors      = NULL,
    flowcode_fluors = NULL,
    fluor_tags      = NULL,
    spectra         = NULL,
    cosine_sim      = NULL,
    output_path     = NULL,
    gate_msg        = NULL,
    # polygon gating
    poly_vertices_x = numeric(0),
    poly_vertices_y = numeric(0),
    poly_closed     = FALSE
  )

  # --------------------------------------------------------------------------
  # Path / volume setup
  # --------------------------------------------------------------------------

  initial_wd <- getShinyOption("threshold_app_wd", getwd())

  default_subfolder <- getShinyOption("threshold_app_outdir",
                                      file.path(initial_wd, "flowcode_spectra"))
  if (!dir.exists(default_subfolder)) {
    tryCatch(dir.create(default_subfolder, recursive = TRUE),
             error = function(e) NULL)
  }
  # Use a plain variable here — rv$ cannot be read/written at top-level server scope
  default_path <- if (dir.exists(default_subfolder)) default_subfolder else initial_wd

  # Initialise rv$output_path inside observe() so it runs in a reactive context
  observe({
    rv$output_path <- default_path
  })

  # Build volumes — only include C: on Windows
  volumes <- c(
    FlowCodeSpectra = default_path,
    WorkingDir      = initial_wd
  )
  if (.Platform$OS.type == "windows") {
    volumes <- c(volumes, `C:` = "C:/")
  }

  # --------------------------------------------------------------------------
  # shinyFiles wiring
  # --------------------------------------------------------------------------

  shinyFileChoose(input, "file",
                  roots    = volumes,
                  session  = session,
                  filetypes = c("fcs", "rds"))

  shinyDirChoose(input, "output_dir",
                 roots   = volumes,
                 session = session)

  observeEvent(input$output_dir, {
    chosen <- parseDirPath(volumes, input$output_dir)
    if (length(chosen) > 0 && nchar(chosen) > 0)
      rv$output_path <- chosen
  })

  output$output_dir_display <- renderText({
    if (is.null(rv$output_path)) "No folder selected"
    else paste("Saving to:", rv$output_path)
  })

  # --------------------------------------------------------------------------
  # File loading
  # --------------------------------------------------------------------------

  output$file_status <- renderText({
    if (is.null(rv$data)) "No file loaded." else {
      nr <- nrow(rv$data)
      nf <- length(rv$all_fluors)
      paste0(format(nr, big.mark = ","), " events  |  ", nf, " channels loaded.")
    }
  })

  load_data_file <- function(file_path, file_name) {

    withProgress(message = "Loading data…", value = 0, {

      # Reset state
      rv$data            <- NULL
      rv$gated_data      <- NULL
      rv$thresholds      <- list()
      rv$cosine_sim      <- NULL
      rv$poly_vertices_x <- numeric(0)
      rv$poly_vertices_y <- numeric(0)
      rv$poly_closed     <- FALSE

      loaded_raw_data <- NULL
      loaded_fluors   <- NULL
      loaded_tags     <- NULL

      incProgress(0.2, detail = "Reading file…")

      tryCatch({
        if (grepl("\\.fcs$", file_name, ignore.case = TRUE)) {
          ff <- flowCore::read.FCS(file_path,
                                   transformation       = FALSE,
                                   truncate_max_range   = FALSE,
                                   emptyValue           = FALSE)
          loaded_raw_data <- base::as.data.frame(flowCore::exprs(ff),
                                                 check.names = FALSE)
          opt_fluors <- getShinyOption("flowcode_fluors_list")
          if (!is.null(opt_fluors)) {
            loaded_fluors <- as.character(opt_fluors)
            loaded_tags   <- names(opt_fluors) %||% loaded_fluors
          } else {
            loaded_fluors <- colnames(loaded_raw_data)
            loaded_tags   <- loaded_fluors
          }
          rv$spectra <- getShinyOption("spectra") %||% NULL

        } else if (grepl("\\.rds$", file_name, ignore.case = TRUE)) {
          rds_obj         <- base::readRDS(file_path)
          loaded_raw_data <- base::as.data.frame(rds_obj$Unmixed, check.names = FALSE)

          opt_fluors <- getShinyOption("flowcode_fluors_list")
          if (!is.null(opt_fluors)) {
            loaded_fluors <- as.character(opt_fluors)
            loaded_tags   <- names(opt_fluors) %||% loaded_fluors
          } else {
            loaded_fluors <- as.character(rds_obj$Flowcode.fluors)
            loaded_tags   <- names(rds_obj$Flowcode.fluors) %||% loaded_fluors
          }
          rv$spectra <- getShinyOption("spectra") %||% rds_obj$Spectra %||% NULL
        } else {
          showNotification("Unsupported file type.", type = "error")
          return()
        }
      }, error = function(e) {
        showNotification(paste("Error loading file:", conditionMessage(e)), type = "error")
        return()
      })

      if (is.null(loaded_raw_data)) return()

      # Remove scatter channels
      fcs_scatter <- c("FSC-A", "SSC-A", "FSC-H", "SSC-H", "FSC-W", "SSC-W", "Time")
      clean_list  <- base::setdiff(loaded_fluors, fcs_scatter)

      rv$data            <- loaded_raw_data
      rv$flowcode_fluors <- loaded_fluors
      rv$all_fluors      <- clean_list
      rv$fluor_tags      <- stats::setNames(loaded_tags, loaded_fluors)
      rv$current_fluor   <- if (length(clean_list) > 0) clean_list[[1]] else NULL

      incProgress(0.5, detail = "Computing spectral similarity…")

      if (!is.null(rv$spectra) &&
          length(clean_list) > 0 &&
          all(clean_list %in% rownames(rv$spectra))) {
        tryCatch(
          rv$cosine_sim <- AutoSpectral::cosine.similarity(
            rv$spectra[clean_list, , drop = FALSE]
          ),
          error = function(e) rv$cosine_sim <- NULL
        )
      }

      updateSelectInput(session, "y_axis_select", choices = clean_list)
      incProgress(1, detail = "Done.")
    })
  }

  observeEvent(input$file, {
    req(input$file)
    file_info <- parseFilePaths(volumes, input$file)
    req(nrow(file_info) > 0)
    load_data_file(as.character(file_info$datapath), as.character(file_info$name))
  })

  # Auto-load a backbone RDS on startup if one was passed to
  # launch.threshold.app(backbone.rds = ...).
  observe({
    backbone_path <- getShinyOption("backbone_rds")
    if (!is.null(backbone_path) && file.exists(backbone_path))
      load_data_file(backbone_path, basename(backbone_path))
  }) |> bindEvent(TRUE, once = TRUE)

  # --------------------------------------------------------------------------
  # Fluorophore selector UI
  # --------------------------------------------------------------------------

  output$fluor_selector_ui <- renderUI({
    req(rv$all_fluors)
    if (length(rv$all_fluors) == 0) return(tags$em("No channels found."))
    radioButtons(
      "current_fluor_select",
      label    = NULL,
      choices  = rv$all_fluors,
      selected = rv$current_fluor %||% rv$all_fluors[[1]]
    )
  })

  observeEvent(input$current_fluor_select, {
    req(input$current_fluor_select)
    rv$current_fluor <- input$current_fluor_select

    # Auto-select the most spectrally dissimilar channel as Y
    if (!is.null(rv$cosine_sim)) {
      fl <- rv$current_fluor
      if (fl %in% rownames(rv$cosine_sim)) {
        sims   <- rv$cosine_sim[fl, ]
        best_y <- names(sort(sims, decreasing = TRUE))[2]
        if (!is.null(best_y))
          updateSelectInput(session, "y_axis_select", selected = best_y)
      }
    }
  })

  # --------------------------------------------------------------------------
  # Retrieve asp (set by launch.threshold.app())
  # --------------------------------------------------------------------------

  asp <- getShinyOption("asp")

  # Seed transform-tuning sliders from asp defaults (if asp is available)
  observe({
    p <- asp$default.transformation.param
    if (!is.null(p)) {
      updateSliderInput(session, "xform_width_x", value = max(p$width %||% -1000, -1000))
      updateNumericInput(session, "xform_pos_x",  value = p$pos   %||%  4.418)
      updateNumericInput(session, "xform_neg_x",  value = p$neg   %||%  0)
      updateSliderInput(session, "xform_width_y", value = max(p$width %||% -1000, -1000))
      updateNumericInput(session, "xform_pos_y",  value = p$pos   %||%  4.418)
      updateNumericInput(session, "xform_neg_y",  value = p$neg   %||%  0)
    }
  }) |> bindEvent(TRUE, once = TRUE)   # run exactly once at session start

  # --------------------------------------------------------------------------
  # Transform parameter helpers (X and Y)
  # --------------------------------------------------------------------------

  get_xform_x <- reactive({
    req(rv$data)
    p    <- asp$default.transformation.param
    len  <- p$length    %||% 256
    maxr <- p$max.range %||% 4194304
    list(
      length    = len,
      max.range = maxr,
      pos       = input$xform_pos_x,
      neg       = input$xform_neg_x,
      width     = input$xform_width_x
    )
  })

  get_xform_y <- reactive({
    req(rv$data)
    if (isTRUE(input$link_axes)) return(get_xform_x())
    p    <- asp$default.transformation.param
    len  <- p$length    %||% 256
    maxr <- p$max.range %||% 4194304
    list(
      length    = len,
      max.range = maxr,
      pos       = input$xform_pos_y,
      neg       = input$xform_neg_y,
      width     = input$xform_width_y
    )
  })

  # --------------------------------------------------------------------------
  # Gate plot  (FSC-A vs SSC-A)
  # --------------------------------------------------------------------------

  gate_plot_trigger <- reactiveVal(0L)
  observeEvent(input$refresh_gate_plot, {
    gate_plot_trigger(gate_plot_trigger() + 1L)
  })
  # Also trigger automatically when data is first loaded
  observeEvent(rv$data, {
    gate_plot_trigger(gate_plot_trigger() + 1L)
  })

  output$gate_plot <- renderPlot({
    gate_plot_trigger()   # take dependency on trigger
    req(rv$data)

    df <- rv$data
    if (!all(c("FSC-A", "SSC-A") %in% colnames(df))) {
      return(ggplot() + theme_void() +
               labs(title = "FSC-A / SSC-A channels not found in data."))
    }

    n_show <- min(nrow(df), 100000L)
    scatter_df <- df[seq_len(n_show), c("FSC-A", "SSC-A")]
    colnames(scatter_df) <- c("x", "y")

    # Use asp scatter limits when available, otherwise fall back to data quantiles
    if (!is.null(asp$scatter.data.min.x)) {
      xlim <- c(asp$scatter.data.min.x, asp$scatter.data.max.x)
      ylim <- c(asp$scatter.data.min.y, asp$scatter.data.max.y)
    } else {
      xlim <- quantile(scatter_df$x, c(0.001, 0.999), na.rm = TRUE)
      ylim <- quantile(scatter_df$y, c(0.001, 0.999), na.rm = TRUE)
    }

    p <- ggplot(scatter_df, aes(x, y)) +
      geom_hex(bins = 80, na.rm = TRUE) +
      scale_fill_viridis_c(option = "turbo") +
      scale_x_continuous(name = "FSC-A",
                         limits = xlim,
                         labels = scales::label_number(scale_cut = scales::cut_si(""))) +
      scale_y_continuous(name = "SSC-A",
                         limits = ylim,
                         labels = scales::label_number(scale_cut = scales::cut_si(""))) +
      theme_bw(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            legend.position  = "none",
            aspect.ratio     = 1)

    gate_title <- if (input$gate_mode == "poly") {
      if (rv$poly_closed)
        "Polygon gate applied — switch tab or refresh to redefine"
      else if (length(rv$poly_vertices_x) > 0)
        paste0("Polygon gate: ", length(rv$poly_vertices_x),
               " vertices — double-click to close")
      else
        "Click to add polygon vertices — double-click to close"
    } else {
      "Draw a rectangular brush gate"
    }

    p <- p + labs(title = gate_title)

    # Draw polygon vertices / outline if in poly mode
    if (input$gate_mode == "poly" && length(rv$poly_vertices_x) > 0) {
      vx <- rv$poly_vertices_x
      vy <- rv$poly_vertices_y
      poly_df <- data.frame(x = vx, y = vy)

      if (rv$poly_closed) {
        # Closed filled polygon
        p <- p +
          geom_polygon(data  = poly_df, aes(x, y),
                       fill  = "#2C6FAC", alpha = 0.20,
                       color = "#2C6FAC", linewidth = 1.2, inherit.aes = FALSE)
      } else {
        # Point markers always; path only when ≥ 2 vertices exist (avoids
        # "each group consists of only one observation" warning from geom_path)
        if (nrow(poly_df) >= 2)
          p <- p +
            geom_path(data     = poly_df, aes(x, y, group = 1),
                      color    = "#e74c3c", linewidth = 1,
                      linetype = "dashed", na.rm = TRUE, inherit.aes = FALSE)
        p <- p +
          geom_point(data  = poly_df, aes(x, y),
                     color = "#e74c3c", size = 3, na.rm = TRUE, inherit.aes = FALSE)
      }
    }

    p
  })

  # On-page gate confirmation message
  output$gate_msg_ui <- renderUI({
    msg <- rv$gate_msg
    if (is.null(msg)) return(NULL)
    tags$div(
      style = paste(
        "margin-top:8px; padding:8px 10px;",
        "background:#d4edda; border:1px solid #c3e6cb;",
        "border-radius:4px; color:#155724; font-size:0.9em;"
      ),
      icon("circle-check"), " ", msg
    )
  })

  # --------------------------------------------------------------------------
  # Gate extraction — rectangular
  # --------------------------------------------------------------------------

  observeEvent(input$gate_brush, {
    if (input$gate_mode != "rect") return()
    req(rv$data)
    brush <- input$gate_brush
    if (is.null(brush)) return()
    gated <- rv$data |>
      dplyr::filter(`FSC-A` >= brush$xmin, `FSC-A` <= brush$xmax,
                    `SSC-A` >= brush$ymin, `SSC-A` <= brush$ymax)
    rv$gated_data <- gated
    n <- nrow(gated)
    tot <- nrow(rv$data)
    rv$gate_msg <- paste0(
      "Rectangular gate applied: ",
      format(n, big.mark = ","), " / ", format(tot, big.mark = ","),
      " cells (", round(100 * n / tot, 1), "%)"
    )
  })

  # --------------------------------------------------------------------------
  # Gate extraction — polygon clicks
  # --------------------------------------------------------------------------

  observeEvent(input$gate_click, {
    if (input$gate_mode != "poly") return()
    if (rv$poly_closed) return()   # don't add after closing
    req(input$gate_click)
    rv$poly_vertices_x <- c(rv$poly_vertices_x, input$gate_click$x)
    rv$poly_vertices_y <- c(rv$poly_vertices_y, input$gate_click$y)
  })

  observeEvent(input$gate_dblclick, {
    if (input$gate_mode != "poly") return()
    req(rv$data)
    vx <- rv$poly_vertices_x
    vy <- rv$poly_vertices_y
    if (length(vx) < 3) {
      showNotification("Need at least 3 vertices to close a polygon gate.",
                       type = "warning")
      return()
    }
    rv$poly_closed <- TRUE
    inside <- points_in_polygon(rv$data[["FSC-A"]], rv$data[["SSC-A"]], vx, vy)
    gated <- rv$data[inside, ]
    rv$gated_data <- gated
    n   <- nrow(gated)
    tot <- nrow(rv$data)
    rv$gate_msg <- paste0(
      "Polygon gate applied (", length(vx), " vertices): ",
      format(n, big.mark = ","), " / ", format(tot, big.mark = ","),
      " cells (", round(100 * n / tot, 1), "%)"
    )
  })

  observeEvent(input$clear_polygon, {
    rv$poly_vertices_x <- numeric(0)
    rv$poly_vertices_y <- numeric(0)
    rv$poly_closed     <- FALSE
    rv$gated_data      <- NULL
    rv$gate_msg        <- NULL
  })

  # Gate stats
  output$gate_stats <- renderText({
    if (is.null(rv$data)) return("No data loaded.")
    if (is.null(rv$gated_data)) return("No gate applied.")
    total  <- nrow(rv$data)
    gated  <- nrow(rv$gated_data)
    pct    <- round(100 * gated / total, 1)
    paste0(format(gated, big.mark = ","), " / ",
           format(total, big.mark = ","), " cells  (", pct, "%)")
  })

  # --------------------------------------------------------------------------
  # Threshold plot  (triggered by button only to avoid rerender on slider drag)
  # --------------------------------------------------------------------------

  thresh_plot_trigger <- reactiveVal(0L)
  observeEvent(input$refresh_threshold_plot, {
    thresh_plot_trigger(thresh_plot_trigger() + 1L)
  })
  # Also re-trigger when the selected channel changes (UX convenience)
  observeEvent(rv$current_fluor, {
    thresh_plot_trigger(thresh_plot_trigger() + 1L)
  })
  observeEvent(rv$gated_data, {
    thresh_plot_trigger(thresh_plot_trigger() + 1L)
  })

  output$threshold_plot <- renderPlot({
    thresh_plot_trigger()   # dependency
    req(rv$gated_data, rv$current_fluor, input$y_axis_select)

    x_fluor <- rv$current_fluor
    y_fluor <- input$y_axis_select

    df <- rv$gated_data

    if (!all(c(x_fluor, y_fluor) %in% colnames(df))) {
      return(ggplot() + theme_void() +
               labs(title = paste("Channels not found:", x_fluor, "or", y_fluor)))
    }

    px <- get_xform_x()
    py <- get_xform_y()

    xf_x <- make_biexp(px$length, px$max.range, px$pos, px$neg, px$width)
    xf_y <- make_biexp(py$length, py$max.range, py$pos, py$neg, py$width)

    if (is.null(xf_x) || is.null(xf_y)) {
      return(ggplot() + theme_void() +
               labs(title = "Invalid transform parameters — adjust sliders."))
    }

    n_show  <- min(nrow(df), 2e5L)
    idx     <- seq_len(n_show)
    plot_df <- data.frame(
      x = xf_x(df[[x_fluor]][idx]),
      y = xf_y(df[[y_fluor]][idx])
    )

    # ---- Axis labels from asp$ribbon.breaks (matching create.biplot) ----------
    # Labels use scientific notation: "0" for zero, "10^N" for powers of 10.
    # Upper limit: asp$expr.data.max (the cytometer's fluorescence ceiling).
    # Lower limit: -1000 by default; pulled further negative when the user has
    #   set neg > 0, using the same inverse-transform approach as create.biplot.
    expr_max <- asp$expr.data.max %||% px$max.range
    ribbon_breaks <- asp$ribbon.breaks   # NULL-safe

    make_axis <- function(xf, neg_decades) {
      raw_max <- expr_max
      # Derive lower display limit: -1000 baseline, extended by neg transform
      # neg_decades is additive on top of the -1000 baseline:
      #   0 -> -1000,  1 -> -10000,  2 -> -100000, ...
      raw_min <- -(10^neg_decades * 1000)
      if (!is.null(ribbon_breaks)) {
        raw_breaks <- ribbon_breaks[ribbon_breaks < raw_max]
      } else {
        raw_breaks <- c(-1e4, -1e3, 0, 1e3, 1e4, 1e5, 1e6, 1e7)
        raw_breaks <- raw_breaks[raw_breaks < raw_max]
      }
      axis_labels <- sapply(raw_breaks, function(v) {
        if (v == 0) "0" else parse(text = paste0("10^", log10(abs(v))))
      })
      list(
        breaks = tryCatch(xf(raw_breaks), error = function(e) NULL),
        labels = axis_labels,
        limits = tryCatch(sort(xf(c(raw_min, raw_max))), error = function(e) NULL)
      )
    }

    ax_x <- make_axis(xf_x, px$neg)
    ax_y <- make_axis(xf_y, py$neg)

    x_tag <- rv$fluor_tags[[x_fluor]] %||% x_fluor
    y_tag <- rv$fluor_tags[[y_fluor]] %||% y_fluor

    # ---- Build plot (matching create.biplot layer order) ----------------------
    p <- ggplot(plot_df, aes(x, y)) +
      scattermore::geom_scattermore(
        pointsize = asp$figure.gate.point.size %||% 0.5,
        color     = "black",
        alpha     = 1,
        na.rm     = TRUE
      ) +
      stat_density_2d(
        aes(fill = after_stat(level)),
        geom  = "polygon",
        na.rm = TRUE
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.margin      = margin(
          asp$figure.margin %||% 5, asp$figure.margin %||% 5,
          asp$figure.margin %||% 5, asp$figure.margin %||% 5
        ),
        legend.position  = "none",
        aspect.ratio     = 1,
        axis.ticks       = element_line(
          linewidth = asp$figure.panel.line.size %||% 0.5),
        axis.text        = element_text(
          size = asp$figure.axis.text.size %||% 11),
        axis.title       = element_text(
          size = asp$figure.axis.title.size %||% 12),
        panel.border     = element_rect(
          fill = NA, linewidth = asp$figure.panel.line.size %||% 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = paste("Threshold:", x_tag, "(", x_fluor, ")"),
        x     = paste0(x_tag, "  [", x_fluor, "]"),
        y     = paste0(y_tag, "  [", y_fluor, "]")
      )

    # ---- Colour fill (matching create.biplot palette logic) -------------------
    viridis_opts <- c("magma","inferno","plasma","viridis",
                      "cividis","rocket","mako","turbo")
    color_palette <- asp$color.palette %||% "rainbow"
    if (color_palette %in% viridis_opts) {
      p <- p + scale_fill_viridis_c(option = color_palette)
    } else {
      p <- p + scale_fill_gradientn(
        colors = asp$density.palette.base.color %||%
          c("#00007F","blue","#007FFF","cyan","#7FFF7F",
            "yellow","#FF7F00","red","#7F0000")
      )
    }

    if (!is.null(ax_x$breaks))
      p <- p + scale_x_continuous(breaks = ax_x$breaks,
                                  labels = ax_x$labels,
                                  limits = ax_x$limits)
    if (!is.null(ax_y$breaks))
      p <- p + scale_y_continuous(breaks = ax_y$breaks,
                                  labels = ax_y$labels,
                                  limits = ax_y$limits)

    # Threshold vline
    thr <- rv$thresholds[[x_fluor]]
    if (!is.null(thr))
      p <- p + geom_vline(xintercept = thr, color = "#e74c3c",
                          linewidth = 1.2, linetype = "dashed")

    p
  })

  # --------------------------------------------------------------------------
  # Record threshold on click AND immediately save the plot
  # (captures the exact transform state; overwrites if user re-clicks)
  # --------------------------------------------------------------------------
  observeEvent(input$threshold_click, {
    req(rv$current_fluor, rv$gated_data)
    cx <- input$threshold_click$x
    if (is.null(cx)) return()

    fl    <- rv$current_fluor
    rv$thresholds[[fl]] <- cx

    # Only attempt to save if an output path exists
    if (is.null(rv$output_path) || !nzchar(rv$output_path)) return()

    x_tag <- rv$fluor_tags[[fl]] %||% fl

    # Choose Y channel (most spectrally dissimilar, or first available)
    y_fl <- if (!is.null(rv$cosine_sim) && fl %in% rownames(rv$cosine_sim)) {
      names(sort(rv$cosine_sim[fl, ], decreasing = TRUE))[2]
    } else {
      setdiff(rv$all_fluors, fl)[1]
    }
    y_tag <- rv$fluor_tags[[y_fl]] %||% y_fl

    tryCatch({
      px   <- get_xform_x()
      py   <- get_xform_y()
      xf_x <- make_biexp(px$length, px$max.range, px$pos, px$neg, px$width)
      xf_y <- make_biexp(py$length, py$max.range, py$pos, py$neg, py$width)
      if (is.null(xf_x) || is.null(xf_y)) return()

      # ---- Transform data ------------------------------------------------
      n_show  <- min(nrow(rv$gated_data), 2e5L)
      plot_df <- data.frame(
        x = xf_x(rv$gated_data[[fl  ]][seq_len(n_show)]),
        y = xf_y(rv$gated_data[[y_fl]][seq_len(n_show)])
      )

      # ---- Axis limits & breaks (same logic as threshold_plot) -----------
      expr_max      <- asp$expr.data.max %||% px$max.range
      ribbon_breaks <- asp$ribbon.breaks
      make_save_axis <- function(xf, neg_decades, max_raw) {
        raw_min <- -(10^neg_decades * 1000)
        if (!is.null(ribbon_breaks)) {
          raw_br <- ribbon_breaks[ribbon_breaks < max_raw]
        } else {
          raw_br <- c(-1e4, -1e3, 0, 1e3, 1e4, 1e5, 1e6, 1e7)
          raw_br <- raw_br[raw_br < max_raw]
        }
        axis_labels <- sapply(raw_br, function(v) {
          if (v == 0) "0" else parse(text = paste0("10^", log10(abs(v))))
        })
        list(
          breaks = tryCatch(xf(raw_br), error = function(e) NULL),
          labels = axis_labels,
          limits = tryCatch(sort(xf(c(raw_min, max_raw))), error = function(e) NULL)
        )
      }
      ax_x <- make_save_axis(xf_x, px$neg, expr_max)
      ax_y <- make_save_axis(xf_y, py$neg, expr_max)

      # ---- Build plot (mirrors threshold_plot exactly) -------------------
      p <- ggplot(plot_df, aes(x, y)) +
        scattermore::geom_scattermore(
          pointsize = asp$figure.gate.point.size %||% 0.5,
          color     = "black",
          alpha     = 1,
          na.rm     = TRUE
        ) +
        stat_density_2d(
          aes(fill = after_stat(level)),
          geom  = "polygon",
          na.rm = TRUE
        ) +
        geom_vline(xintercept = cx, color = "#e74c3c",
                   linewidth = 1.2, linetype = "dashed") +
        theme_bw(base_size = 13) +
        theme(
          plot.margin      = margin(
            asp$figure.margin %||% 5, asp$figure.margin %||% 5,
            asp$figure.margin %||% 5, asp$figure.margin %||% 5
          ),
          legend.position  = "none",
          aspect.ratio     = 1,
          axis.ticks       = element_line(
            linewidth = asp$figure.panel.line.size %||% 0.5),
          axis.text        = element_text(
            size = asp$figure.axis.text.size %||% 11),
          axis.title       = element_text(
            size = asp$figure.axis.title.size %||% 12),
          panel.border     = element_rect(
            fill = NA, linewidth = asp$figure.panel.line.size %||% 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        labs(
          title = paste("Threshold:", x_tag, "(", fl, ")"),
          x     = paste0(x_tag, "  [", fl, "]"),
          y     = paste0(y_tag, "  [", y_fl, "]")
        )

      # ---- Colour fill ---------------------------------------------------
      viridis_opts  <- c("magma","inferno","plasma","viridis",
                         "cividis","rocket","mako","turbo")
      color_palette <- asp$color.palette %||% "rainbow"
      if (color_palette %in% viridis_opts) {
        p <- p + scale_fill_viridis_c(option = color_palette)
      } else {
        p <- p + scale_fill_gradientn(
          colors = asp$density.palette.base.color %||%
            c("#00007F","blue","#007FFF","cyan","#7FFF7F",
              "yellow","#FF7F00","red","#7F0000")
        )
      }

      if (!is.null(ax_x$breaks))
        p <- p + scale_x_continuous(breaks = ax_x$breaks,
                                    labels = ax_x$labels,
                                    limits = ax_x$limits)
      if (!is.null(ax_y$breaks))
        p <- p + scale_y_continuous(breaks = ax_y$breaks,
                                    labels = ax_y$labels,
                                    limits = ax_y$limits)

      # ---- Save (overwrites previous plot for this channel) ---------------
      outname <- file.path(rv$output_path,
                           paste0(x_tag, "_", fl, "_threshold.png"))
      ggplot2::ggsave(outname, p, width = 6, height = 6, dpi = 200)

    }, error = function(e) {
      message("Threshold plot save failed for ", fl, ": ", conditionMessage(e))
    })
  })

  # Threshold table
  output$threshold_table <- renderTable({
    if (length(rv$thresholds) == 0) return(NULL)
    data.frame(
      Channel           = names(rv$thresholds),
      Transform_thresh  = round(unlist(rv$thresholds), 4),
      stringsAsFactors  = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # --------------------------------------------------------------------------
  # Save thresholds + plots
  # --------------------------------------------------------------------------

  observeEvent(input$save_thr, {
    if (length(rv$thresholds) == 0) {
      showNotification("No thresholds to save.", type = "warning")
      return()
    }
    if (is.null(rv$output_path) || !nzchar(rv$output_path)) {
      showNotification("Please choose an output folder first.", type = "error")
      return()
    }
    if (is.null(rv$gated_data)) {
      showNotification("No gated data — apply a gate first.", type = "error")
      return()
    }

    fname <- input$out_filename
    if (!grepl("\\.csv$", fname, ignore.case = TRUE))
      fname <- paste0(fname, ".csv")

    # Build inverse transform to recover raw threshold values
    px <- get_xform_x()
    biexp_inverse_fn <- make_biexp(px$length, px$max.range,
                                   px$pos, px$neg, px$width,
                                   inverse = TRUE)

    thresholds_raw <- if (!is.null(biexp_inverse_fn)) {
      sapply(rv$thresholds, function(x) as.numeric(biexp_inverse_fn(x)))
    } else {
      rep(NA_real_, length(rv$thresholds))
    }

    df_out <- data.frame(
      Fluor                 = names(rv$thresholds),
      Tag                   = rv$fluor_tags[names(rv$thresholds)] %||% names(rv$thresholds),
      Threshold_Raw         = thresholds_raw,
      Threshold_Transformed = unlist(rv$thresholds),
      stringsAsFactors      = FALSE
    )

    outfile  <- file.path(rv$output_path, fname)
    write_ok <- tryCatch({
      utils::write.csv(df_out, outfile, row.names = FALSE)
      TRUE
    }, error = function(e) {
      showNotification(paste("Could not write CSV:", conditionMessage(e)),
                       type = "error")
      FALSE
    })
    if (!write_ok) return()

    showNotification(
      paste0("Saved: ", fname, "  (plots saved automatically at each threshold click)"),
      type     = "message",
      duration = 8
    )
  })
}

# ==============================================================================
# Run
# ==============================================================================
shiny::shinyApp(ui, server)
