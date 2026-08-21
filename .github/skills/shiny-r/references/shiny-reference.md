# Shiny (R) — Full Reference

Source: Mastering Shiny (Hadley Wickham) — https://mastering-shiny.org

---

## Core Structure

Every app lives in `app.R` (or split into `ui.R` + `server.R`):

```r
library(shiny)

ui <- fluidPage(
  # front-end: inputs and output placeholders
)

server <- function(input, output, session) {
  # back-end: reactive logic
}

shinyApp(ui, server)
```

Run with `runApp("myapp/")` or RStudio's Run App button (Cmd/Ctrl+Shift+Enter).

---

## UI — Inputs

All inputs share `inputId` (unique string, letters/numbers/underscores only) and `label`.
Access in server via `input$<inputId>`.

| Function | Use |
|---|---|
| `textInput()` | Single-line text |
| `textAreaInput()` | Multi-line text |
| `numericInput()` | Number box |
| `sliderInput()` | Slider (range slider if `value = c(a,b)`) |
| `selectInput()` | Dropdown (set `multiple = TRUE` for multi-select) |
| `radioButtons()` | Radio group |
| `checkboxGroupInput()` | Checkbox group |
| `checkboxInput()` | Single checkbox (yes/no) |
| `dateInput()` / `dateRangeInput()` | Date picker |
| `fileInput()` | File upload |
| `actionButton()` / `actionLink()` | Trigger events |

`actionButton` classes: `"btn-primary"`, `"btn-success"`, `"btn-danger"`, `"btn-lg"`, `"btn-sm"`, `"btn-block"`.

---

## UI — Outputs (placeholders only; paired with render functions in server)

| UI function | Server render | Use |
|---|---|---|
| `textOutput()` | `renderText()` | Inline text (uses `cat`) |
| `verbatimTextOutput()` | `renderPrint()` | Console-style text (uses `print`) |
| `plotOutput()` | `renderPlot()` | Plots (base/ggplot2) |
| `tableOutput()` | `renderTable()` | Static table |
| `DT::DTOutput()` | `DT::renderDT()` | Interactive table (preferred over deprecated `dataTableOutput`) |
| `downloadButton()` / `downloadLink()` | `downloadHandler()` | File download |
| `leaflet::leafletOutput()` | `leaflet::renderLeaflet()` | Interactive map |
| `plotly::plotlyOutput()` | `plotly::renderPlotly()` | Interactive plot |

`plotOutput()` also accepts `click`, `dblclick`, `hover` → creates `input$<id>` for interactive plots.

---

## Server — Reactivity

**Key principle:** Shiny is *declarative* — you describe recipes, not commands. Code runs when needed, not top-to-bottom.

### Inputs
- `input$id` is read-only; must be accessed inside a reactive context (render function or `reactive()`).
- Each session gets its own independent `input`/`output` environment.

### Outputs
- `output$id <- render*({ ... })` — the render function creates a reactive context that tracks input dependencies.
- Outputs are atomic: if any tracked input changes, the entire output re-runs.

### Reactive expressions — `reactive({})`
- Cache their result; only re-run when their inputs change.
- Called like functions: `my_expr()`.
- Use to share computed values between multiple outputs (avoids duplication and redundant computation).
- **Rule of one**: if you copy reactive code once, extract it into a `reactive()`.

```r
dataset <- reactive({ get(input$dataset, "package:datasets") })
output$summary <- renderPrint({ summary(dataset()) })
output$table   <- renderTable({ dataset() })
```

### Observers — `observeEvent(event, handler)`
- For side effects (writing files, DB updates, console messages).
- NOT assigned to a variable; result cannot be used reactively.
- Similar sibling: `observe({})` (runs whenever any reactive dependency changes).

### Event-driven computation — `eventReactive(event, expr)`
- Like `reactive()` but only updates when `event` fires (e.g. `actionButton`).
- Decouples dependency from value used in computation.

```r
x1 <- eventReactive(input$simulate, { rpois(input$n, input$lambda1) })
```

### Modern alternative — `bindEvent()`
Preferred in Shiny ≥ 1.6 over `eventReactive`/`observeEvent`:

```r
x1 <- reactive({ rpois(input$n, input$lambda1) }) |> bindEvent(input$simulate)
observeEvent(input$simulate, { message("clicked") })
# equivalent to:
observe({ message("clicked") }) |> bindEvent(input$simulate)
```

### Timed invalidation — `reactiveTimer(ms)`
- Creates a reactive input that invalidates every `ms` milliseconds.

### Isolate — `isolate(expr)`
- Read a reactive value without creating a dependency on it.

### Reactive values — `reactiveValues()`
- Mutable R list where assignments trigger downstream updates.

```r
rv <- reactiveValues(count = 0)
observeEvent(input$btn, { rv$count <- rv$count + 1 })
output$n <- renderText({ rv$count })
```

---

## Reactive graph concepts
- **Producers**: inputs + reactive expressions
- **Consumers**: reactive expressions + outputs (+ observers)
- Execution order is determined by the reactive graph, NOT code order.
- Use `reactlog` package to visualise the live reactive graph: `reactlog::reactlog_enable(); shinyApp(ui, server)` then Ctrl+F3.

---

## Layout

### Classic (Bootstrap grid)
```r
fluidPage(
  fluidRow(
    column(4, ...),   # width out of 12
    column(8, ...)
  )
)
```

### Modern — bslib (recommended)
```r
library(bslib)

page_sidebar(
  title = "App title",
  sidebar = sidebar(
    # inputs
  ),
  # main content: cards, plots, tables
  card(
    card_header("Plot"),
    plotOutput("plot")
  )
)
```

Other `bslib` layouts: `page_navbar()`, `page_fillable()`, `layout_columns()`, `layout_sidebar()`.

`value_box()` — KPI-style summary boxes.

### Tabbed layouts
```r
navbarPage("Title",
  tabPanel("Tab 1", ...),
  tabPanel("Tab 2", ...)
)
# or inside fluidPage:
tabsetPanel(
  tabPanel("A", ...),
  tabPanel("B", ...)
)
```

---

## Modules (for larger apps)

Encapsulate UI + server logic into reusable, namespace-isolated units.

```r
# Module UI
myModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    selectInput(ns("var"), "Variable", choices = names(mtcars)),
    plotOutput(ns("plot"))
  )
}

# Module server
myModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$plot <- renderPlot({ hist(mtcars[[input$var]]) })
  })
}

# In app
ui <- fluidPage(myModuleUI("mod1"))
server <- function(input, output, session) { myModuleServer("mod1") }
```

Pass reactive values into modules as function arguments (not `input$x` directly).

---

## Validation & error handling
- `req(input$x)` — silently stop if input is missing/NULL/FALSE.
- `validate(need(condition, "Error message"))` — show user-friendly error in output panel.
- `tryCatch()` / `withCallingHandlers()` for robust error handling in render functions.

---

## Dynamic UI
- `uiOutput()` / `renderUI()` — render UI elements dynamically from the server.
- `updateSelectInput()`, `updateSliderInput()`, etc. — update existing inputs from server.
- `showModal()` / `modalDialog()` — modal dialogs.
- `showNotification()` — toast-style notifications.

---

## File upload / download

```r
# Upload
fileInput("file", "Upload CSV")
# In server:
df <- reactive({ req(input$file); read.csv(input$file$datapath) })

# Download
downloadButton("dl", "Download")
# In server:
output$dl <- downloadHandler(
  filename = "results.csv",
  content = function(file) write.csv(my_data(), file)
)
```

---

## Deployment options
- **shinyapps.io** (cloud, Posit-hosted): `rsconnect::deployApp()`
- **Shiny Server** (self-hosted, open source)
- **Posit Connect** (enterprise)
- **Local**: share `app.R` + `renv.lock`; user runs `runApp()`

---

## Performance tips
- Use `reactive()` to cache expensive computations.
- `bindCache()` (Shiny ≥ 1.6) — persistent caching across sessions, keyed by inputs.
- `renderCachedPlot()` for plot caching.
- Profile with `profvis::profvis({ runApp() })`.
- Move data loading and model fitting **outside** `server()` so it runs once at startup, not per session.

---

## Useful packages

| Package | Purpose |
|---|---|
| `bslib` | Modern Bootstrap themes & layouts |
| `DT` | Interactive tables (`DTOutput()` / `renderDT()`) |
| `plotly` | Interactive plots |
| `leaflet` | Interactive maps |
| `shinyWidgets` | Extended input widgets |
| `waiter` / `shinycssloaders` | Loading spinners |
| `shinyjs` | JS helpers from R |
| `reactlog` | Visualise reactive graph |
| `thematic` | Auto-theme ggplot2/base plots to match app theme |

---

## Common pitfalls
- Accessing `input$x` outside a reactive context → error
- Mismatched IDs between UI and server (typos cause **silent** failures — Shiny is lazy)
- Using plain variables or plain functions instead of `reactive()` for shared reactive values
- Naming reactive expressions the same as base R functions (e.g. `range`, `var`)
- Putting expensive data loading inside `server()` — runs on every new session
- Forgetting `()` when calling a reactive expression: `dataset` vs `dataset()`
