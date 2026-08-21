---
applyTo: "production/**,analysis/**"
---

## Skills to load

- When writing or modifying INLA model code, cross-validation, or covariate selection → load the `r-data-science` skill.
- When working with covariance matrices, derivable nodes, flow conservation, greedy MI selection, or WITH/AGAINST direction pairing → load the `sensor-placement` skill.

## Output naming convention

New result files must follow: `<weighting>_<k>percounty_<r>overall.rds`, saved under `results/production/`.

Examples:
- `results/production/frc_sqrt(0.9_to_0)_100percounty_100overall.rds`
- `results/production/frc_sqrt(0.9_to_0)_and_logAADT_0.9_to_0.1_weight_100percounty_100overall.rds`

## After producing a new result file

Both of these edits are required — neither is optional:
1. Add an entry to the `configs` named list in `shiny/app.R` (path: `results/production/<file>.rds`)
2. Add the same filename (basename only) to `rds_files` in `deploy.R`

Omitting step 2 causes the deployed app to crash on startup (exit status 1 — file not found in bundle).

## Function reference

All function signatures are documented in `/memories/repo/xdt-sensor-locations-structure.md`.
Algorithm internals are documented in `/memories/repo/greedy-mi-sensor-selection-algorithm.md`.
