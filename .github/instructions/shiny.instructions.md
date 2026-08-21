---
applyTo: "shiny/**,deploy.R"
---
Before making any changes to Shiny app files or `deploy.R`, read `.github/skills/shiny-r/SKILL.md` in full.

## Deploy rule

**Never add an entry to `configs` in `app.R` without simultaneously adding the same filename to the `rds_files` list in `deploy.R`.** These two edits must always happen together.

Files missing from `deploy.R` are excluded from the rsconnect bundle and crash the app on startup (exit status 1 — `readRDS` on a missing file).
