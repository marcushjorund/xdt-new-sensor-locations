---
name: shiny-r
description: 'Shiny R expert. Use when: building Shiny apps in R; writing ui/server code; reactive programming with reactive() observeEvent() eventReactive(); layouts with bslib fluidPage navbarPage; deploying Shiny apps to shinyapps.io; debugging reactivity; creating interactive dashboards for R analysis results; displaying plots tables maps in a Shiny app; Shiny modules; file upload download in Shiny.'
argument-hint: 'Describe the Shiny component or task to build'
---

# Shiny R Expert

## When to Use
- Building or modifying a Shiny app (`app.R`, `ui.R`/`server.R`)
- Writing reactive logic (`reactive()`, `observeEvent()`, `eventReactive()`, `bindEvent()`)
- Designing UI layouts (sidebar, tabs, cards, grid)
- Displaying maps, plots, or tables interactively
- Debugging why outputs don't update (reactivity issues)
- Deploying a Shiny app

## Full Reference

See [shiny-reference.md](./references/shiny-reference.md) for the complete API reference covering:
- All input/output functions and their render-function pairings
- Reactivity patterns and the reactive graph
- Layout with `fluidPage` / `bslib`
- Modules, validation, performance, deployment
- Useful companion packages (`DT`, `leaflet`, `plotly`, `bslib`)
- Common pitfalls

## Quick Start

```r
library(shiny)
library(bslib)

ui <- page_sidebar(
  title = "My App",
  sidebar = sidebar(
    selectInput("choice", "Choose:", choices = c("A", "B"))
  ),
  plotOutput("plot")
)

server <- function(input, output, session) {
  output$plot <- renderPlot({
    plot(1:10, main = input$choice)
  })
}

shinyApp(ui, server)
```

## Key Principles
1. `input$x` is read-only and must be accessed inside a reactive context
2. Use `reactive()` to share computed values between outputs — never copy reactive code
3. Execution order follows the **reactive graph**, not line order
4. `req(input$x)` to silently stop when inputs are missing

## Deployment with a staging deploy.R

When a project uses a `deploy.R` script that manually copies files into `shiny/` before calling `rsconnect::deployApp("shiny/")`, **every data or result file referenced in `app.R` must also be listed in the staging section of `deploy.R`**. Files not listed there are excluded from the bundle and will cause the app to crash on startup (`readRDS` / `read*` on a missing file → `exit status 1`). This includes:
- RDS result files added to a `configs` list
- Any new data files loaded at the top of `app.R`
