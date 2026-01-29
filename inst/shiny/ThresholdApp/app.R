# app.R

# Allow large uploads (e.g., 2 GB)
options(shiny.maxRequestSize = 2 * 1024^3)

# ---- UI ----
ui <- shiny::fluidPage(
  shiny::tags$head(
    shiny::tags$style(HTML("
      #fluor_container { max-height: 300px; overflow-y: auto; background: #f9f9f9;
                         padding: 10px; border: 1px solid #ccc; border-radius: 4px; }
      .shiny-notification { position: fixed; top: 10%; right: 10%; }
    "))
  ),

  shiny::titlePanel("Interactive Thresholding App"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shinyFiles::shinyFilesButton(
        "file",
        "Load FCS or RDS",
        "Please select a file",
        multiple = FALSE
      ),
      shiny::tags$hr(),

      shiny::tags$b("Select fluor to threshold (X-axis):"),
      shiny::tags$div(id = "fluor_container",
                      shiny::uiOutput("fluor_selector_ui")
      ),

      # Y-Axis Override
      shiny::selectInput("y_axis_select", "Y-axis Fluorophore:", choices = NULL),

      shiny::tags$hr(),
      shinyFiles::shinyDirButton(
        "output_dir",
        "Change output folder",
        "Select folder"
      ),
      shiny::verbatimTextOutput("output_dir_display"),

      # User Filename Input
      shiny::textInput("out_filename", "Output CSV Filename:", value = "thresholds.csv"),

      shiny::actionButton("save_thr", "Save thresholds"),
      shiny::tags$hr(),
      shiny::tableOutput("threshold_table")
    ),
    shiny::mainPanel(
      shiny::plotOutput(
        "gate_plot",
        brush = shiny::brushOpts(id = "gate_brush", resetOnNew = TRUE)
      ),
      shiny::plotOutput("threshold_plot", click = "threshold_click")
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  # initialize reactive values
  rv <- shiny::reactiveValues(
    data = NULL,
    gated_data = NULL,
    thresholds = list(),
    current_fluor = NULL,
    all_fluors = NULL,
    flowcode_fluors = NULL,
    fluor_tags = NULL,
    spectra = NULL,
    cosine_sim = NULL,
    output_path = NULL,
  )

  # Recover environment options
  initial_wd <- getOption("threshold_app_wd") %||% getwd()

  # --- STARTUP PATH LOGIC ---
  default_subfolder <- file.path(initial_wd, "flowcode_spectra")
  if (!dir.exists(default_subfolder)) dir.create(default_subfolder)
  rv$output_path <- default_subfolder
  volumes <- c(
    FlowCodeSpectra = default_subfolder,
    WorkingDir = initial_wd,
    `C:` = "C:/"
  )

  output$output_dir_display <- shiny::renderText({
    if (is.null(rv$output_path)) "No folder selected" else paste("Saving to:", rv$output_path)
  })

  # Use working directory as root for folder chooser + optionally C: drive
  shinyFiles::shinyFileChoose(
    input,
    id = "file",
    roots = volumes,
    session = session,
    filetypes = c("fcs", "rds") # Limit file types to FCS and RDS
  )

  shinyFiles::shinyDirChoose(
    input,
    id = "output_dir",
    roots = volumes,
    session = session
  )

  # --- FILE LOADING ---
  observeEvent(input$file, {
    shiny::req(input$file)
    file_info <- shinyFiles::parseFilePaths(volumes, input$file)
    req(nrow(file_info) > 0)

    file_path <- file_info$datapath
    file_name <- file_info$name
    file_size_mb <- file.info(file_path)$size / 1024^2

    if (file_size_mb > 10) {
      shiny::showModal(modalDialog(
        title = "Large File Warning",
        paste0("This file is ", round(file_size_mb, 1), " MB. Processing..."),
        easyClose = TRUE, footer = modalButton("Dismiss")
      ))
    }

    shiny::withProgress(message = 'Loading data...', value = 0, {
      rv$gated_data <- NULL
      rv$thresholds <- list()
      incProgress(0.2, detail = "Reading file...")

      # Local variables to hold data before assigning to rv
      loaded_raw_data <- NULL
      loaded_fluors <- NULL
      loaded_tags <- NULL

      if (grepl("\\.fcs$", file_name, ignore.case = TRUE)) {
        ff <- flowCore::read.FCS(file_path, transformation = FALSE, truncate_max_range = FALSE, emptyValue = FALSE)
        loaded_raw_data <- base::as.data.frame(flowCore::exprs(ff), check.names = FALSE)

        # Priority 1: Wrapper Option
        opt_fluors <- getOption("flowcode_fluors_list")
        if (!is.null(opt_fluors)) {
          #loaded_fluors <- if(!all(grepl("-A$", opt_fluors))) paste0(opt_fluors, "-A") else opt_fluors
          loaded_fluors <- opt_fluors
          loaded_tags <- names(opt_fluors) %||% opt_fluors
        } else {
          # Priority 2: File Columns
          loaded_fluors <- colnames(loaded_raw_data)
          loaded_tags <- loaded_fluors
        }
        rv$spectra <- getOption("Spectra") %||% NULL

      } else if (grepl("\\.rds$", file_name, ignore.case = TRUE)) {
        rds_obj <- base::readRDS(file_path)
        loaded_raw_data <- base::as.data.frame(rds_obj$Unmixed, check.names = FALSE)
        loaded_fluors <- as.character(rds_obj$Flowcode.fluors)
        loaded_tags <- names(rds_obj$Flowcode.fluors) %||% loaded_fluors
        rv$spectra <- rds_obj$Spectra %||% NULL
      }

      # Clean List (Remove Scatter)
      fcs_cols <- c("FSC-A", "SSC-A", "FSC-H", "SSC-H", "Time")
      clean_list <- base::setdiff(loaded_fluors, fcs_cols)

      # Assign to Reactive Values
      rv$data <- loaded_raw_data
      rv$flowcode_fluors <- loaded_fluors
      rv$all_fluors <- clean_list # This drives the search/radio menu
      rv$fluor_tags <- setNames(loaded_tags, loaded_fluors)
      rv$current_fluor <- clean_list[1]

      # Similarity logic
      if (!is.null(rv$spectra) && all(clean_list %in% rownames(rv$spectra))) {
        incProgress(0.5, detail = "Calculating similarity...")
        rv$cosine_sim <- AutoSpectral::cosine.similarity(rv$spectra[clean_list, , drop = FALSE])
      }

      shiny::updateSelectInput(session, "y_axis_select", choices = clean_list)
      setProgress(1, detail = "Done")
    })
  })

  # --- Fluorophore selection ---
  output$fluor_selector_ui <- shiny::renderUI({
    shiny::req(rv$all_fluors) # Only render once data is actually loaded

    # Safeguard: if for some reason the list is empty, don't crash
    if (length(rv$all_fluors) == 0) return("No channels found.")

    shiny::radioButtons(
      "current_fluor_select",
      label = NULL,
      choices = rv$all_fluors,
      selected = rv$current_fluor %||% rv$all_fluors[1]
    )
  })

  # --- X-AXIS -> Y-AXIS AUTO UPDATE ---
  observeEvent(input$current_fluor_select, {
    shiny::req(input$current_fluor_select)
    rv$current_fluor <- input$current_fluor_select

    if (!is.null(rv$cosine_sim)) {
      sims <- rv$cosine_sim[rv$current_fluor, ]
      best_y <- names(sort(sims, decreasing = TRUE)[2])
      shiny::updateSelectInput(session, "y_axis_select", selected = best_y)
    }
  })

  # ---- Gating plot ----
  output$gate_plot <- shiny::renderPlot({
    req(rv$data)
    df <- rv$data

    if (nrow(df) > 2e5)
      scatter_df <- df[1:2e5,]
    else
      scatter_df <- df

    ggplot(scatter_df, aes(x = `FSC-A`, y = `SSC-A`)) +
      geom_scattermore(size = 0.5, na.rm = TRUE) +
      stat_density_2d(
        aes(fill = after_stat(level)),
        geom = "polygon",
        contour = TRUE,
        na.rm = TRUE) +
      scale_fill_viridis_c(option = "turbo") +
      scale_x_continuous(
        #name = x.lab,
        breaks = seq(
          asp$scatter.data.min.x,
          asp$scatter.data.max.x,
          asp$data.step
        ),
        labels = paste0(
          round(
            seq(
              asp$scatter.data.min.x,
              asp$scatter.data.max.x,
              asp$data.step
            ) / 1e6,
            1
          ),
          "e6"
        ),
        limits = c(asp$scatter.data.min.x, asp$scatter.data.max.x),
        expand = expansion(asp$figure.gate.scale.expand)
      ) +
      scale_y_continuous(
        #name = y.lab,
        breaks = seq(
          asp$scatter.data.min.y,
          asp$scatter.data.max.y,
          asp$data.step
        ),
        labels = paste0(
          round(
            seq(
              asp$scatter.data.min.y,
              asp$scatter.data.max.y,
              asp$data.step
            ) / 1e6, 1
          ),
          "e6"
        ),
        limits = c(asp$scatter.data.min.y, asp$scatter.data.max.y),
        expand = expansion(asp$figure.gate.scale.expand)
      ) +
      labs(title = "Draw a rectangular gate") +
      theme_classic() +
      theme(legend.position = "none")
  })

  # Gate extraction
  observeEvent(input$gate_brush, {
    req(rv$data)
    brush <- input$gate_brush
    df <- rv$data

    # Ensure brush exists (should)
    if (is.null(brush)) return()

    rv$gated_data <- df |>
      dplyr::filter(`FSC-A` >= brush$xmin, `FSC-A` <= brush$xmax,
                    `SSC-A` >= brush$ymin, `SSC-A` <= brush$ymax)

    # When gating, reset thresholds (optional) or keep existing - here we keep them
  })

  # ---- Threshold plotting ----
  output$threshold_plot <- shiny::renderPlot({
    req(rv$gated_data, rv$current_fluor, input$y_axis_select)

    x.fluor <- rv$current_fluor
    # pull out currently selected y fluorophore
    y.fluor <- input$y_axis_select

    # match tags to fluorophores
    x.tag <- rv$fluor_tags[x.fluor]
    y.tag <- rv$fluor_tags[y.fluor]

    p <- AutoSpectral::create.biplot(
      rv$gated_data,
      x.fluor,
      y.fluor,
      asp,
      x.lab = paste(x.tag, x.fluor),
      y.lab = paste(y.tag, y.fluor),
      x.min = asp$default.transformation.param$width,
      y.min = asp$default.transformation.param$width,
      x.width.basis = asp$default.transformation.param$width,
      y.width.basis = asp$default.transformation.param$width,
      title = paste("Thresholding:", x.tag, x.fluor),
      save = FALSE
    )

    # Add vertical threshold line if exists
    if (!is.null(rv$thresholds[[x.fluor]])) {
      p <- p + geom_vline(
        xintercept = rv$thresholds[[x.fluor]],
        color = "red",
        linewidth = 1
      )
    }

    p
  })


  # Record threshold on click
  observeEvent(input$threshold_click, {
    req(rv$current_fluor)
    if (is.null(input$threshold_click$x)) return()
    rv$thresholds[[rv$current_fluor]] <- input$threshold_click$x
  })

  # Show thresholds
  output$threshold_table <- renderTable({
    if (length(rv$thresholds) == 0) return(NULL)
    data.frame(
      Fluor = names(rv$thresholds),
      Transformed_Threshold = unlist(rv$thresholds),
      stringsAsFactors = FALSE
    )
  })

  # Save thresholds
  observeEvent(input$save_thr, {
    req(length(rv$thresholds) > 0)
    # require the user to have chosen an output folder
    if (is.null(rv$output_path) || length(rv$output_path) == 0) {
      showNotification("Please choose an output folder before saving.", type = "error")
      return()
    }

    # Ensure the filename ends in .csv
    fname <- input$out_filename
    if (!grepl("\\.csv$", fname, ignore.case = TRUE)) {
      fname <- paste0(fname, ".csv")
    }

    # Construct the full path using the custom name
    outfile <- file.path(rv$output_path, fname)

    # Ensure asp exists (needed for transforms)
    if (!exists("asp")) {
      showNotification("asp is not available in the app environment. Cannot perform inverse transform.", type = "error")
      return()
    }

    # Build inverse transform functions (bi-exp helpers)
    biexp.transform <- flowWorkspace::flowjo_biexp(
      channelRange = asp$default.transformation.param$length,
      maxValue     = asp$default.transformation.param$max.range,
      pos          = asp$default.transformation.param$pos,
      neg          = asp$default.transformation.param$neg,
      widthBasis   = asp$default.transformation.param$width,
      inverse      = FALSE
    )

    biexp.inverse <- flowWorkspace::flowjo_biexp(
      channelRange = asp$default.transformation.param$length,
      maxValue     = asp$default.transformation.param$max.range,
      pos          = asp$default.transformation.param$pos,
      neg          = asp$default.transformation.param$neg,
      widthBasis   = asp$default.transformation.param$width,
      inverse      = TRUE
    )

    # Invert thresholds before saving
    thresholds_raw <- sapply(rv$thresholds, function(x) {
      # ensure numeric
      as.numeric(biexp.inverse(x))
    })

    df_out <- data.frame(
      Fluor = names(thresholds_raw),
      #Tag = rv$fluor_tags[match(names(thresholds_raw), rv$flowcode_fluors)],
      Tag = rv$fluor_tags[names(thresholds_raw)],
      Threshold_Raw = thresholds_raw,
      Threshold_Transformed = unlist(rv$thresholds),
      stringsAsFactors = FALSE
    )

    outfile <- file.path(rv$output_path, fname)
    utils::write.csv(df_out, outfile, row.names = FALSE)

    # Save per-fluor plots with line
    for (fl in names(rv$thresholds)) {

      x.tag <- rv$fluor_tags[match(fl, rv$flowcode_fluors)]

      y.fluor <- if (!is.null(rv$cosine_sim)) {
        names(sort(rv$cosine_sim[fl, ], decreasing = TRUE)[2])
      } else {
        setdiff(rv$flowcode_fluors, fl)[1]
      }

      y.tag <- rv$fluor_tags[match(y.fluor, rv$flowcode_fluors)]

      p <- AutoSpectral::create.biplot(
        rv$gated_data, fl, y.fluor, asp,
        x.lab = paste(x.tag, fl),
        y.lab = paste(y.tag, y.fluor),
        x.min = asp$default.transformation.param$width,
        y.min = asp$default.transformation.param$width,
        x.width.basis = asp$default.transformation.param$width,
        y.width.basis = asp$default.transformation.param$width,
        title = paste(x.tag, fl, "vs", y.tag, y.fluor),
        save = FALSE
      ) +
        geom_vline(
          xintercept = rv$thresholds[[fl]],
          color = "red", linewidth = 1
        )

      # Save using antigen tag as filename prefix
      outname <- paste0(x.tag, "_", fl, "_threshold_plot.png")

      ggsave(file.path(rv$output_path, outname), p, height = 6, width = 6, dpi = 300)
    }

    shiny::showNotification(paste("Thresholds and plots saved to", rv$output_path))
  })
}

# ---- Run App ----
shiny::shinyApp(ui, server)

