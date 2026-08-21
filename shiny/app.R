# =============================================================================
# Sensor Placement Visualisation — Shiny App
# =============================================================================
# Run from project root: shiny::runApp("shiny/")
# =============================================================================

library(shiny)
library(bslib)
library(leaflet)
library(DT)

# ── Path setup ────────────────────────────────────────────────────────────────
# Shiny (and withr::with_dir inside reactive evaluations) resets the working
# directory to the app folder. All project paths are therefore resolved as
# absolute paths from project_root to be wd-independent.
#
# Local dev:    getwd() = shiny/       → results/ is one level up
# shinyapps.io: getwd() = app dir      → results/ is copied alongside app.R
#                                         by deploy.R before bundling
project_root <- if (file.exists(file.path(getwd(), "results"))) {
  normalizePath(getwd())                   # deployed: data is right here
} else {
  normalizePath(file.path(getwd(), "..")) # local dev: data is in parent dir
}

# ── Source self-contained mapping helpers ─────────────────────────────────────
# sensor_map_helpers.R contains plot_sensor_selection_map, add_geometries,
# map_traffic_links, and nvdb_objects — no INLA or xdtkit required.
source(file.path(getwd(), "sensor_map_helpers.R"), local = FALSE)

# ── Tell add_geometries() where the GeoJSON lives (absolute path via option) ──
# options() are global and survive Shiny's internal withr::with_dir() calls,
# so this is more robust than patching formals() or using setwd().
options(sensor_app.geojson_path = file.path(
  project_root, "data", "directed_traffic_links_2024_simplified.geojson"
))

# ── Load measured traffic links for map overlay ──────────────────────────────
all_data_measured <- readRDS(file.path(
  project_root, "data", "measured_traffic_links_minimal.rds"
))

# =============================================================================
# Configuration registry
# Add a new configuration here — no other changes needed.
# Each entry: display_label = relative path from project root
# =============================================================================
configs <- list(
  "Funksjonell vegklasse-vekting (kvadratrot, 0,9→0) — 100 per fylke, 100 totalt" =
    "results/production/frc_sqrt(0.9_to_0)_100percounty_100overall.rds",

  "Funksjonsklasse-vekting (kvadratrot, 0,6→0) — 100 per fylke, 100 totalt" =
    "results/production/functionClass_sqrt(0.6_to_0)_100percounty_100overall.rds",

  "Uvektet — 100 per fylke, 100 totalt" =
    "results/production/unweighted_100percounty100overall.rds",

  "ÅDT-vekting (proporsjonal) — 100 per fylke, 100 totalt" =
    "results/production/AADT_prop_100percounty_100overall.rds",

  "ÅDT-vekting (log-proporsjonal) — 100 per fylke, 100 totalt" =
    "results/production/AADT_logprop_100percounty_100overall.rds",

  "ÅDT-vekting (uten endring) — 100 per fylke, 100 totalt" =
    "results/production/AADT_identity_100percounty_100overall.rds",

  "Funksjonsklasse-vekting (kvadratrot, 0,6→0) + log-prop. ÅDT + tyngdeandel (α: 0,8/0,15/0,05) — 100 per fylke, 100 totalt" =
    "results/production/functionClass_sqrt(0.6_to_0)logpropaadt_heavy_08015005_100percounty_100overall.rds",

  "Funksjonell vegklasse-vekting (kvadratrot, 0,9→0) + log-prop. ÅDT + tyngdeandel (α: 0,8/0,15/0,05) — 100 per fylke, 100 totalt" =
    "results/production/frc_sqrt(0.9_to_0)logpropaadt_heavy_08015005_100percounty_100overall.rds",

  "Uvektet bundlet MI (rå, uten ekstra informasjon) — 100 per fylke, 100 totalt" =
    "results/production/unweighted_bundle_100percounty_100overall.rds",

  "Funksjonell vegklasse + log-prop. ÅDT + tyngdeandel (α: 0,6/0,3/0,1) — 100 per fylke, 100 totalt" =
    "results/production/frc_0.6_logaadt_0.3_heavyratio_0.1_100percounty_100overall.rds",

  # ── Fixed-algorithm comparison (dynamic-sum neighbourhood weighting) ──────────────
  "Uvektet (distribuert, sum-vekting) — 100 per fylke, 100 totalt" =
    "results/production/unweighted_distsum_100percounty_100overall.rds",

  "Funksjonell vegklasse (kvadratrot, 0,9→0, distribuert sum) — 100 per fylke, 100 totalt" =
    "results/production/frc_sqrt(0.9_to_0)_distsum_100percounty_100overall.rds",

  "Uvektet bundlet MI (distribuert sum) — 100 per fylke, 100 totalt" =
    "results/production/unweighted_bundle_distsum_100percounty_100overall.rds"
)

# ── Load all results at startup (outside server — shared across sessions) ─────
results <- lapply(configs, function(rel_path) {
  readRDS(file.path(project_root, rel_path))
})
names(results) <- names(configs)

# ── Covariate display labels (for summary dropdown) ──────────────────────────
covariate_labels <- c(
  "functionalRoadClass" = "Funksjonell vegklasse",
  "functionClass"       = "Funksjonsklasse",
  "roadCategory"        = "Vegkategori",
  "county"              = "Fylke",
  "aadt_category"       = "ÅDT-kategori"
)

# ── AADT category display labels (Norwegian, with band thresholds) ────────────
# Keys = English internal labels stored in .rds files (from categorise_aadt()).
# Values = Norwegian display labels shown in the UI.
aadt_level_labels <- c(
  "Very Low"  = "Veldig lav \u00c5DT (<1\u00a0000)",
  "Low"       = "Lav \u00c5DT (1\u00a0000\u20135\u00a0000)",
  "Medium"    = "Middels \u00c5DT (5\u00a0000\u201315\u00a0000)",
  "High"      = "H\u00f8y \u00c5DT (15\u00a0000\u201335\u00a0000)",
  "Very High" = "Veldig h\u00f8y \u00c5DT (>35\u00a0000)"
)

# ── Attach county to measured data from prepared_traffic_links_norway ─────────
# all_data_measured only carries id/aadt/traffic_volume_source.
# prepared_traffic_links_norway.rds has county for all 13 000+ link IDs, so we
# use it as the authoritative lookup instead of pooling from selected_data_entries
# (which only covers the ~800 optimisation-candidate links, leaving 94% as NA).
{
  lkp_prepared <- readRDS(file.path(project_root, "data", "prepared_traffic_links_norway.rds"))
  if ("county" %in% names(lkp_prepared)) {
    all_data_measured$county <- lkp_prepared$county[
      match(all_data_measured$id, lkp_prepared$id)
    ]
  } else {
    message("sensor-app: 'county' column not found in prepared_traffic_links_norway.rds — measured-links county filter disabled.")
  }
  rm(lkp_prepared)
}

# =============================================================================
# UI
# =============================================================================
ui <- page_sidebar(
  title = "Sensorplassering — Norge",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  # ── Sidebar ────────────────────────────────────────────────────────────────
  sidebar = sidebar(
    width = 310,

    h6("Konfigurasjon", class = "text-uppercase text-muted fw-semibold mt-1 mb-1"),
    selectInput(
      inputId  = "config",
      label    = NULL,
      choices  = names(configs),
      selected = names(configs)[[1]]
    ),

    hr(class = "my-2"),

    uiOutput("scalar_summary")
  ),

  # ── Main content ───────────────────────────────────────────────────────────
  navset_card_tab(
    nav_panel(
      title = "Kart",

      # ── Collapsible filter panel ──────────────────────────────────────────
      accordion(
        id   = "map_filter_accordion",
        open = FALSE,
        accordion_panel(
          title = "Filtrer kart",
          layout_column_wrap(
            width = 1/3,

            # ── Col 1: County ──────────────────────────────────────────────
            div(uiOutput("county_filter_ui")),

            # ── Col 2: Vegkategori + Funksjonsklasse ───────────────────────
            div(
              checkboxGroupInput(
                inputId  = "filter_roadCategory",
                label    = "Vegkategori",
                choices  = c(
                  "Kommunal veg" = "KOMMUNAL_VEG",
                  "Fylkesveg"    = "FYLKESVEG",
                  "Riksveg"      = "RIKSVEG",
                  "Europaveg"    = "EUROPAVEG"
                ),
                selected = c("KOMMUNAL_VEG", "FYLKESVEG", "RIKSVEG", "EUROPAVEG")
              ),
              checkboxGroupInput(
                inputId  = "filter_functionClass",
                label    = "Funksjonsklasse",
                choices  = c("A", "B", "C", "D", "E", "Ukjent" = "unknown"),
                selected = c("A", "B", "C", "D", "E", "unknown"),
                inline   = TRUE
              )
            ),

            # ── Col 3: ÅDT-kategori + Funksjonell vegklasse ────────────────
            div(
              checkboxGroupInput(
                inputId  = "filter_aadt",
                label    = "\u00c5DT-kategori",
                choices  = aadt_level_labels,
                selected = names(aadt_level_labels)
              ),
              checkboxGroupInput(
                inputId  = "filter_frc",
                label    = "Funksjonell vegklasse",
                choices  = as.character(0:9),
                selected = as.character(0:9),
                inline   = TRUE
              )
            )
          ),

          div(
            class = "mt-2",
            actionButton(
              inputId = "reset_filters",
              label   = "Tilbakestill filter",
              class   = "btn-sm btn-outline-secondary"
            )
          )
        )
      ),

      leafletOutput("map", height = "640px")
    ),
    nav_panel(
      title = "Sammendragstabell",
      tags$div(
        class = "p-3 pb-0",
        uiOutput("covariate_filter")
      ),
      DT::DTOutput("summary_table")
    )
  )
)

# =============================================================================
# Server
# =============================================================================
server <- function(input, output, session) {

  # ── Central reactive: current result object ─────────────────────────────────
  selected_result <- reactive({
    results[[input$config]]
  })

  # ── County filter UI (rebuilt when config changes) ───────────────────────────
  output$county_filter_ui <- renderUI({
    df       <- selected_result()$selected_data_entries
    counties <- if ("county" %in% names(df) && is.factor(df$county))
      sort(levels(df$county))
    else if ("county" %in% names(df))
      sort(unique(as.character(df$county[!is.na(df$county)])))
    else
      character(0)
    selectizeInput(
      inputId  = "filter_county",
      label    = "Fylke",
      choices  = c("Alle fylker", counties),
      selected = "Alle fylker",
      multiple = TRUE,
      options  = list(plugins = list("remove_button"))
    )
  })

  # ── Filter reactive: subsets selected_data_entries by all active filters ─────
  filtered_entries <- reactive({
    df <- selected_result()$selected_data_entries

    # County
    county_sel <- input$filter_county
    if (!is.null(county_sel) && length(county_sel) > 0 &&
        !"Alle fylker" %in% county_sel && "county" %in% names(df)) {
      df <- df[!is.na(df$county) & df$county %in% county_sel, , drop = FALSE]
    }

    # Vegkategori
    rc_all <- c("KOMMUNAL_VEG", "FYLKESVEG", "RIKSVEG", "EUROPAVEG")
    rc_sel <- input$filter_roadCategory
    if (!is.null(rc_sel) && length(rc_sel) > 0 &&
        !setequal(rc_sel, rc_all) && "roadCategory" %in% names(df)) {
      df <- df[is.na(df$roadCategory) | df$roadCategory %in% rc_sel, , drop = FALSE]
    }

    # ÅDT-kategori — derive category from raw AADT, filter on English internal code
    aadt_all <- names(aadt_level_labels)
    aadt_sel <- input$filter_aadt  # checkboxGroupInput returns English keys
    if (!is.null(aadt_sel) && length(aadt_sel) > 0 &&
        !setequal(aadt_sel, aadt_all) && "lastYearAadt_aadt" %in% names(df)) {
      cats <- as.character(cut(
        df$lastYearAadt_aadt,
        breaks         = c(0, 1000, 5000, 15000, 35000, Inf),
        labels         = names(aadt_level_labels),
        include.lowest = TRUE,
        right          = FALSE
      ))
      df <- df[is.na(cats) | cats %in% aadt_sel, , drop = FALSE]
    }

    # Funksjonell vegklasse
    frc_all <- as.character(0:9)
    frc_sel <- input$filter_frc
    if (!is.null(frc_sel) && length(frc_sel) > 0 &&
        !setequal(frc_sel, frc_all) && "functionalRoadClass" %in% names(df)) {
      df <- df[is.na(df$functionalRoadClass) |
                 as.character(df$functionalRoadClass) %in% frc_sel, , drop = FALSE]
    }

    # Funksjonsklasse
    fc_all <- c("A", "B", "C", "D", "E", "unknown")
    fc_sel <- input$filter_functionClass
    if (!is.null(fc_sel) && length(fc_sel) > 0 &&
        !setequal(fc_sel, fc_all) && "functionClass" %in% names(df)) {
      df <- df[is.na(df$functionClass) | df$functionClass %in% fc_sel, , drop = FALSE]
    }

    df
  })

  # ── Filtered measured data (county filter only) ───────────────────────────────
  filtered_measured <- reactive({
    m          <- all_data_measured
    county_sel <- input$filter_county
    if (!is.null(county_sel) && length(county_sel) > 0 &&
        !"Alle fylker" %in% county_sel && "county" %in% names(m)) {
      m <- m[!is.na(m$county) & m$county %in% county_sel, , drop = FALSE]
    }
    m
  })

  # ── Reset all map filters ─────────────────────────────────────────────────────
  observeEvent(input$reset_filters, {
    updateSelectizeInput(session, "filter_county",
      selected = "Alle fylker")
    updateCheckboxGroupInput(session, "filter_roadCategory",
      selected = c("KOMMUNAL_VEG", "FYLKESVEG", "RIKSVEG", "EUROPAVEG"))
    updateCheckboxGroupInput(session, "filter_functionClass",
      selected = c("A", "B", "C", "D", "E", "unknown"))
    updateCheckboxGroupInput(session, "filter_aadt",
      selected = names(aadt_level_labels))
    updateCheckboxGroupInput(session, "filter_frc",
      selected = as.character(0:9))
  })

  # ── Scalar summary panel ────────────────────────────────────────────────────
  output$scalar_summary <- renderUI({
    res <- selected_result()
    sm  <- res$summary

    items <- list(
      tags$div(class = "d-flex justify-content-between mb-1",
        tags$span(class = "text-muted", "Valgte sensorer"),
        tags$strong(sm$n_selected)
      ),
      tags$div(class = "d-flex justify-content-between mb-1",
        tags$span(class = "text-muted", "Sensorer per fylke (rangeres og filtreres):"),
        tags$strong(res$k_per_county)
      ),
      tags$div(class = "d-flex justify-content-between mb-1",
        tags$span(class = "text-muted", "Antall sensorer i kart (hele Norge, etter filtrering):"),
        tags$strong(if (is.null(res$r)) "—" else res$r)
      )
    )

    # Optional derivable-node columns (only present in some configurations)
    if (!is.null(sm$enables_derivable_count) && sm$enables_derivable_count > 0) {
      items <- c(items, list(
        tags$div(class = "d-flex justify-content-between mb-1",
          tags$span(class = "text-muted", "Sensorer som muliggjør utledbare lenker"),
          tags$strong(sm$enables_derivable_count)
        )
      ))
    }
    if (!is.null(sm$total_derivable_links_enabled) && sm$total_derivable_links_enabled > 0) {
      items <- c(items, list(
        tags$div(class = "d-flex justify-content-between mb-1",
          tags$span(class = "text-muted", "Totalt antall utledbare lenker aktivert"),
          tags$strong(sm$total_derivable_links_enabled)
        )
      ))
    }

    tagList(
      h6("Kjøreparametere", class = "text-uppercase text-muted fw-semibold mb-2"),
      tags$div(class = "small", items)
    )
  })

  # ── Covariate filter (updates when config changes) ────────────────────────
  output$covariate_filter <- renderUI({
    ct              <- selected_result()$summary$counts_table
    available_covs  <- unique(ct$covariate)
    display_labels  <- ifelse(
      available_covs %in% names(covariate_labels),
      covariate_labels[available_covs],
      available_covs
    )
    cov_choices <- c("Alle kovariater" = "all",
                     setNames(available_covs, display_labels))
    selectInput(
      inputId  = "covariate_filter",
      label    = "Vis kovariat:",
      choices  = cov_choices,
      selected = "all"
    )
  })

  # ── Map ─────────────────────────────────────────────────────────────────────
  output$map <- renderLeaflet({
    df <- filtered_entries()
    if (nrow(df) == 0) {
      leaflet::leaflet() |>
        leaflet::addTiles() |>
        leaflet::addControl(
          html     = "<div style='padding:8px 12px;background:white;
                          border-radius:4px;border:1px solid #ccc'>
                        <strong>Ingen data \u00e5 vise for valgte filter.</strong>
                      </div>",
          position = "topright"
        )
    } else {
      plot_sensor_selection_map(df, all_data = filtered_measured())
    }
  })

  # ── Summary table ───────────────────────────────────────────────────────────
  output$summary_table <- DT::renderDT({
    ct <- selected_result()$summary$counts_table

    # Filter by selected covariate (NULL on first render before uiOutput resolves)
    sel_cov <- input$covariate_filter
    if (!is.null(sel_cov) && sel_cov != "all") {
      ct <- ct[ct$covariate == sel_cov, , drop = FALSE]
    }

    # Translate AADT category levels to Norwegian display labels
    aadt_rows <- ct$covariate == "aadt_category"
    if (any(aadt_rows)) {
      ct$level[aadt_rows] <- ifelse(
        ct$level[aadt_rows] %in% names(aadt_level_labels),
        aadt_level_labels[ct$level[aadt_rows]],
        ct$level[aadt_rows]
      )
    }

    # Apply Norwegian labels to covariate column
    ct$covariate_label <- ifelse(
      ct$covariate %in% names(covariate_labels),
      covariate_labels[ct$covariate],
      ct$covariate
    )

    df <- ct[, c("covariate_label", "level", "n", "pct")]
    df <- df[order(df$covariate_label, -df$n), ]
    names(df) <- c("Kovariat", "Nivå", "Antall", "Prosent (%)")

    DT::datatable(
      df,
      rownames  = FALSE,
      options   = list(
        pageLength = 50,
        dom        = "t",      # table only — no search/pagination controls
        ordering   = TRUE
      )
    ) |>
      DT::formatRound("Prosent (%)", digits = 1)
  })
}

# =============================================================================
shinyApp(ui, server)
