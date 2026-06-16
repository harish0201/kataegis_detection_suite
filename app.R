library(shiny)
library(data.table)
library(DT)
library(shinycssloaders)
options(shiny.maxRequestSize = 500 * 1024^200)

source("kataegis_detect2.R")
source("plots.R")

ui <- fluidPage(
  titlePanel("Kataegis Detection Suite"),
  tags$div(
    id = "status_bar",
    style = "margin: 6px 0 10px 0;",
    uiOutput("status_ui")
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("maf_file", "Upload MAF File",
                accept = c(".maf", ".tsv", ".txt"),
                multiple = FALSE),
      fileInput("chrom_lengths_file", "Chromosome Lengths (optional; .txt/.tsv/.fai)",
                accept = c(".tsv", ".txt", ".fai"),
                multiple = FALSE),

      h5("Analysis Parameters"),
      checkboxInput("dynamic_cutoff", "Use Dynamic Cutoff (chromosome-wise)", FALSE),
      numericInput("min_mutations", "Minimum Mutations per Cluster",
                   value = 6, min = 1, max = 1000000),
      numericInput("max_avg_distance", "Maximum Avg Intermutation Distance (bp)",
                   value = 1000, min = 1, max = 10000000),

      uiOutput("sample_selector_ui"),

      # Run mode selector
      radioButtons(
        "run_mode",
        label = "Detection Mode",
        choices = c(
          "Per-sample (individual results)" = "per_sample",
          "Aggregated (cohort-level pooling)"  = "aggregated"
        ),
        selected = "per_sample"
      ),
      tags$small(
        class = "text-muted",
        tags$em("Aggregated: pools all selected samples into one pseudo-sample (katdetectr-style cohort detection)")
      ),
      hr(),

      actionButton("run_btn", "▶  Run Analysis",
                   class = "btn btn-primary btn-block",
                   style = "width:100%; font-weight:bold;"),
      hr(),

      h5("Plot Settings"),
      checkboxInput("log_scale", "Log Scale Y Axis", TRUE),
      checkboxInput("facet_chrom", "Facet by Chromosome", FALSE),
      hr(),

      downloadButton("download_plot", "Download Plot (PNG)"),
      downloadButton("download_tsv",  "Download Kataegis TSV")
    ),

    mainPanel(
      h4("Summary"),
      verbatimTextOutput("summary_text"),
      h4("Kataegis Plot"),
      uiOutput("stale_plot_banner"),
      withSpinner(plotOutput("kataegis_plot", height = "900px"), type = 6),
      h4("Detected Regions"),
      uiOutput("stale_table_banner"),
      div(
        style = "overflow-x: auto; width: 100%;",
        DT::dataTableOutput("regions_table")
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    stage = "idle",
    message = "Upload a MAF file to begin.",
    cached_df = NULL,
    cached_plot_file = NULL,
    last_run_mode = NULL,      # mode actually used in the last Run click
    last_run_n_samples = NULL,
    last_run_samples = NULL    # character vector of barcodes used in last Run
  )

  # ---- Status bar ----------------------------------------------------------
  output$status_ui <- renderUI({
    col <- switch(rv$stage,
                  idle = "info", maf_loaded = "info",
                  running = "warning", done = "success",
                  error = "danger", "info")
    icon_lbl <- switch(rv$stage,
                       idle = "○  Idle", maf_loaded = "✔  MAF Loaded",
                       running = "⧗  Running...", done = "✔  Complete",
                       error = "✘  Error", "○")
    div(
      class = paste0("alert alert-", col),
      style = "padding: 8px 14px; margin-bottom: 4px; font-size: 13px;",
      tags$strong(icon_lbl), " — ", rv$message
    )
  })

  # ---- 1. Read MAF ---------------------------------------------------------
  maf_data <- reactive({
    req(input$maf_file)
    rv$stage <- "running"; rv$message <- "Reading MAF file..."

    ext <- tolower(tools::file_ext(input$maf_file$name))
    if (!ext %in% c("maf", "tsv", "txt")) {
      rv$stage <- "error"; rv$message <- "Unsupported MAF extension."
      showNotification(rv$message, type = "error"); return(NULL)
    }
    dt <- tryCatch(fread(input$maf_file$datapath),
                   error = function(e) {
                     rv$stage <- "error"
                     rv$message <- paste("Error reading MAF:", e$message)
                     showNotification(rv$message, type = "error"); NULL })
    if (is.null(dt)) return(NULL)

    req_cols <- c("Chromosome", "Start_Position", "End_Position",
                  "Reference_Allele", "Tumor_Seq_Allele2",
                  "Tumor_Sample_Barcode", "Variant_Type", "Variant_Classification")
    missing <- setdiff(req_cols, colnames(dt))
    if (length(missing) > 0) {
      msg <- paste("MAF missing columns:", paste(missing, collapse = ", "))
      rv$stage <- "error"; rv$message <- msg
      showNotification(msg, type = "error"); return(NULL)
    }
    n_samp <- length(unique(dt$Tumor_Sample_Barcode))
    n_muts <- nrow(dt)
    rv$stage <- "maf_loaded"
    rv$message <- sprintf("MAF loaded: %d sample(s), %d mutation(s). Set parameters and click Run.",
                          n_samp, n_muts)
    dt
  })

  # ---- 2. Chromosome sizes ------------------------------------------------
  chrom_sizes <- reactive({
    if (is.null(input$chrom_lengths_file)) return(NULL)
    ext <- tolower(tools::file_ext(input$chrom_lengths_file$name))
    cs <- tryCatch(fread(input$chrom_lengths_file$datapath),
                   error = function(e) {
                     showNotification(paste("Error reading chrom sizes:", e$message), type = "error"); NULL })
    if (is.null(cs)) return(NULL)
    if (ext == "fai") cs <- cs[, 1:2, with = FALSE]
    setnames(cs, c(colnames(cs)[1], colnames(cs)[2]), c("Chromosome", "Size"))
    cs
  })

  # ---- 3. Sample multi-select UI ------------------------------------------
  output$sample_selector_ui <- renderUI({
    df <- maf_data(); if (is.null(df)) return(NULL)
    samples <- sort(unique(df$Tumor_Sample_Barcode))
    tagList(
      selectizeInput(
        "sample_choice", label = "Select Sample(s)",
        choices = samples, selected = NULL, multiple = TRUE,
        options = list(
          placeholder = "Leave blank to include All samples",
          plugins = list("remove_button"),
          closeAfterSelect = FALSE
        )
      ),
      tags$small(class = "text-muted", sprintf("(%d sample(s) in MAF)", length(samples)))
    )
  })

  # ---- 4. Active sample resolution ----------------------------------------
  active_samples <- reactive({
    choice <- input$sample_choice
    if (is.null(choice) || length(choice) == 0) return(NULL)
    choice
  })

  selected_or_all_samples <- reactive({
    df <- maf_data(); req(df)
    sids <- active_samples()
    if (is.null(sids)) sort(unique(df$Tumor_Sample_Barcode)) else sids
  })

  filtered_maf <- reactive({
    df <- maf_data(); req(df)
    sids <- active_samples()
    if (!is.null(sids)) df <- df[Tumor_Sample_Barcode %in% sids]
    if (nrow(df) == 0) return(NULL)
    df
  })

  # ---- 5. Build temp MAF path — apply aggregation label if needed ----------
  filtered_maf_path <- reactive({
    df <- filtered_maf()
    if (is.null(df)) {
      rv$stage <- "error"; rv$message <- "No mutations found for selected sample(s)."
      showNotification(rv$message, type = "error"); return(NULL)
    }
    # In aggregated mode replace all barcodes with a single pseudo-label
    if (input$run_mode == "aggregated") {
      n_samp <- length(unique(df$Tumor_Sample_Barcode))
      df <- copy(df)
      df[, Tumor_Sample_Barcode := paste0("Aggregated (", n_samp, " samples)")]
    }
    tmp <- tempfile(fileext = ".tsv")
    fwrite(df, tmp, sep = "\t")
    tmp
  })

  # ---- 6. Detection --------------------------------------------------------
  detection_results <- eventReactive(input$run_btn, {
    req(maf_data())
    fpath <- filtered_maf_path(); req(fpath)

    if (is.na(input$min_mutations) || input$min_mutations < 1) {
      rv$stage <- "error"; rv$message <- "min_mutations must be >= 1."
      showNotification(rv$message, type = "error"); return(NULL)
    }
    if (is.na(input$max_avg_distance) || input$max_avg_distance < 1) {
      rv$stage <- "error"; rv$message <- "max_avg_distance must be >= 1."
      showNotification(rv$message, type = "error"); return(NULL)
    }

    n_active <- length(selected_or_all_samples())
    mode_lbl <- if (input$run_mode == "aggregated")
      sprintf("Aggregated (%d samples pooled)", n_active)
    else
      sprintf("Per-sample across %d sample(s)", n_active)

    rv$stage <- "running"
    rv$message <- sprintf("Step 1/2: Running kataegis detection — %s...", mode_lbl)
    rv$cached_df <- NULL
    rv$cached_plot_file <- NULL
    rv$last_run_mode <- input$run_mode
    rv$last_run_n_samples <- n_active
    rv$last_run_samples <- selected_or_all_samples()

    results <- tryCatch(
      detect_kataegis(
        maf_file = fpath,
        sample_id = NULL,       # always NULL — per-sample loop OR single pseudo-sample
        min_mutations = input$min_mutations,
        max_avg_distance = input$max_avg_distance,
        use_dynamic_cutoff = input$dynamic_cutoff,
        chrom_sizes = chrom_sizes()
      ),
      error = function(e) {
        rv$stage <- "error"; rv$message <- paste("Detection error:", e$message)
        showNotification(rv$message, type = "error"); NULL
      }
    )
    if (is.null(results)) return(NULL)

    rv$cached_df <- if (length(results) == 0) NULL else kataegis_to_dataframe(results)
    rv$stage <- "running"
    rv$message <- sprintf("Step 2/2: Building plot... (%d region(s) detected)", length(results))
    results
  })

  # ---- 7. Detection dataframe (cached) ------------------------------------
  detection_df <- reactive({
    if (!is.null(rv$cached_df)) return(rv$cached_df)
    res <- detection_results()
    if (is.null(res) || length(res) == 0) return(NULL)
    kataegis_to_dataframe(res)
  })

  # ---- 8. Summary ----------------------------------------------------------
  # Shared reactive: TRUE when settings changed since last Run
  pending_change <- reactive({
    if (is.null(rv$last_run_mode)) return(FALSE)
    mode_changed    <- rv$last_run_mode != input$run_mode
    samples_changed <- !identical(sort(rv$last_run_samples),
                                  sort(selected_or_all_samples()))
    mode_changed || samples_changed
  })

  stale_banner <- function() {
    div(
      class = "alert alert-warning",
      style = "padding: 6px 12px; margin-bottom: 6px; font-size: 12px;",
      tags$strong("⚠️ Settings changed since last run."),
      " Click ▶ Run to recompute."
    )
  }

  output$stale_plot_banner  <- renderUI({ if (pending_change()) stale_banner() })
  output$stale_table_banner <- renderUI({ if (pending_change()) stale_banner() })

  output$summary_text <- renderPrint({
    df_sub <- filtered_maf()
    if (is.null(df_sub)) { cat("No MAF loaded."); return() }
    n_samp <- length(unique(df_sub$Tumor_Sample_Barcode))
    n_muts <- nrow(df_sub)
    df_det <- detection_df()

    # Warn user if mode has changed since last run
    pending_change <- !is.null(rv$last_run_mode) &&
                      (rv$last_run_mode != input$run_mode)

    if (is.null(rv$last_run_mode)) {
      mode_str <- "(not yet run)"
      detected <- "Run analysis first"
    } else {
      mode_str <- if (rv$last_run_mode == "aggregated")
        sprintf("Aggregated (cohort-level, %d samples pooled)", rv$last_run_n_samples)
      else
        sprintf("Per-sample (%d sample(s))", rv$last_run_n_samples)
      detected <- if (is.null(df_det)) 0L else nrow(df_det)
    }

    cat(sprintf(
      "Detection mode  : %s\nSamples in scope: %d\nMutations in scope: %d\nKataegis regions detected: %s",
      mode_str, n_samp, n_muts, as.character(detected)
    ))

    if (pending_change)
      cat("\n\n** Mode changed — click Run to recompute **")
  })

  # ---- 9. Plot helper — draws one sample -----------------------------------
  draw_one_sample <- function(sid, maf_path, plot_mode, chroms) {
    res <- detection_results()
    sample_regions <- if (!is.null(res) && length(res) > 0)
      Filter(function(x) !is.null(x$sample) && identical(x$sample, sid), res)
    else NULL
    plot_rainfall(
      maf_file = maf_path,
      sample_id = sid,
      highlight_kataegis = sample_regions,
      plot_mode = plot_mode,
      log_scale = input$log_scale,
      chromosomes = chroms,
      chrom_sizes = chrom_sizes(),
      include_other_mutations = FALSE
    )
  }

  # ---- 10. Shared plot-building logic (screen + PNG export) ----------------
  build_plot <- function(png_path, width = 1600, res = 150) {
    req(detection_results())
    fpath <- filtered_maf_path(); req(fpath)

    plot_mode <- if (isTRUE(input$facet_chrom)) "facetted" else "genome"
    chroms    <- NULL
    df_det    <- detection_df()
    if (!is.null(df_det) && isTRUE(input$facet_chrom))
      chroms <- unique(df_det$chromosome)

    # Read sample IDs from the actual temp MAF — respects aggregated pseudo-label
    sample_ids <- sort(unique(data.table::fread(
      fpath, select = "Tumor_Sample_Barcode")[[1]]))
    n_samples  <- length(sample_ids)

    # --- Render each sample to its own temp PNG independently ---------------
    # This avoids par(mfrow) being reset by plot_rainfall's internal on.exit.
    sample_pngs <- vapply(sample_ids, function(sid) {
      tmp <- tempfile(fileext = ".png")
      tryCatch({
        png(tmp, width = width, height = 700, res = res)
        draw_one_sample(sid, fpath, plot_mode, chroms)
        dev.off()
      }, error = function(e) {
        try(dev.off(), silent = TRUE)
        rv$stage   <- "error"
        rv$message <- paste("Plot error for", sid, ":", e$message)
        showNotification(rv$message, type = "error")
        tmp <- NA_character_
      })
      tmp
    }, FUN.VALUE = character(1))

    # Drop any failed panels
    sample_pngs <- sample_pngs[!is.na(sample_pngs) & file.exists(sample_pngs)]
    if (length(sample_pngs) == 0) return(invisible(NULL))

    # --- Stitch panels into a single output PNG using grid ------------------
    ncol_layout <- if (n_samples <= 2) 1L else 2L
    nrow_layout <- ceiling(length(sample_pngs) / ncol_layout)
    total_height <- 700 * nrow_layout

    png(png_path, width = width, height = total_height, res = res)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      layout = grid::grid.layout(nrow_layout, ncol_layout)
    ))

    for (i in seq_along(sample_pngs)) {
      img  <- png::readPNG(sample_pngs[i])
      row_i <- ceiling(i / ncol_layout)
      col_i <- ((i - 1) %% ncol_layout) + 1
      grid::pushViewport(grid::viewport(
        layout.pos.row = row_i,
        layout.pos.col = col_i
      ))
      grid::grid.raster(img)
      grid::popViewport()
    }
    dev.off()
    png_path
  }

  # ---- 11. renderPlot — builds once, caches PNG path ----------------------
  output$kataegis_plot <- renderPlot({
    req(detection_results())
    tmp_png <- tempfile(fileext = ".png")
    build_plot(tmp_png)
    rv$cached_plot_file <- tmp_png

    img <- png::readPNG(tmp_png)
    grid::grid.raster(img)

    n_samples <- if (input$run_mode == "aggregated") 1 else length(selected_or_all_samples())
    rv$stage <- "done"
    rv$message <- sprintf(
      "Analysis complete. %d kataegis region(s) detected across %d sample(s).",
      if (is.null(detection_df())) 0L else nrow(detection_df()), n_samples
    )
  })

  # ---- 12. Regions table --------------------------------------------------
  output$regions_table <- DT::renderDataTable({
    df <- detection_df()
    if (is.null(df) || nrow(df) == 0) {
      datatable(data.frame(Message = "No kataegis regions detected."),
                options = list(dom = "t"))
    } else {
      datatable(df,
                options = list(pageLength = 10, scrollX = TRUE,
                               scrollY = "300px", autoWidth = FALSE),
                style = "bootstrap")
    }
  })

  # ---- 13. Download TSV (cached df) ----------------------------------------
  output$download_tsv <- downloadHandler(
    filename = function() paste0("kataegis_regions_", Sys.Date(), ".tsv"),
    content = function(file) {
      df <- rv$cached_df
      if (is.null(df)) writeLines("No kataegis regions detected.", con = file)
      else fwrite(df, file, sep = "\t")
    }
  )

  # ---- 14. Download PNG (cached file copy — no recomputation) -------------
  output$download_plot <- downloadHandler(
    filename = function() paste0("kataegis_plot_", Sys.Date(), ".png"),
    content = function(file) {
      if (is.null(rv$cached_plot_file) || !file.exists(rv$cached_plot_file)) {
        showNotification("No plot available yet. Run analysis first.", type = "warning")
        return()
      }
      file.copy(rv$cached_plot_file, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui, server)
