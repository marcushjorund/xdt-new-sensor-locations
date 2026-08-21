# Code Refactor Plan — scripts_inla_sensor_placement.R

## Overview

Split `serious-testing/scripts_inla_sensor_placement.R` (19 functions, ~2700 lines) into focused
sub-files under `R/`, delete dead/superseded functions, document with roxygen2, and clean up obsolete
testing scripts. Always run the refactor agent in plan mode first — no files are deleted without
explicit per-item confirmation.

---

## Current Function Inventory

All 19 functions in `scripts_inla_sensor_placement.R` with their line numbers and status:

| # | Function | Lines | Status | Notes |
|---|----------|-------|--------|-------|
| 1 | `fit_weighted_inla_model` | 1–463 | ✅ keep | Core INLA fitting |
| 2 | `plot_inla_model` | 464–585 | ✅ keep | Diagnostic plots |
| 3 | `kfold_cv_inla` | 586–701 | ✅ keep | Cross-validation |
| 4 | `plot_kfold_cv` | 702–867 | ✅ keep | CV diagnostic plots |
| 5 | `select_similarity_covariates` | 868–1029 | ✅ keep | Forward covariate selection |
| 6 | `identify_derivable_nodes` | 1030–1159 | ✅ keep | Flow conservation logic |
| 7 | `create_covariance_and_precision_matrix` | 1160–1257 | ✅ keep | Core covariance builder |
| 8 | `create_covariance_and_precision_matrix_both_directions` | 1258–1337 | ⚠️ audit | May be superseded by `with_and_against = TRUE` param |
| 9 | `compute_weights` | 1338–1414 | ✅ keep | Internal weight helper |
| 10 | `greedy_mi_sensor_selection` | 1415–1572 | ✅ keep | Core greedy MI algorithm |
| 11 | `expand_border` | 1573–1592 | ⚠️ duplicate | Same function defined in `partition_norway.R` line 55 |
| 12 | `partition_by_county` | 1593–1631 | ⚠️ duplicate | Same function defined in `partition_norway.R` line 76 |
| 13 | `categorise_aadt` | 1632–1665 | ✅ keep | AADT category labels |
| 14 | `summarise_sensor_selection` | 1666–1785 | ✅ keep | Result summarisation |
| 15 | `greedy_mi_sensor_selection_norway` | 1786–2088 | ✅ keep | Norway-specific wrapper |
| 16 | `plot_sensor_selection_map` | 2089–2442 | ✅ keep | Selection map, used in 12 scripts |
| 17 | `add_geometries` | 2443–2495 | ✅ keep | Attaches sf geometries |
| 18 | `map_traffic_links` | 2496–2665 | ✅ keep | General link map |
| 19 | `greedy_mi_sensor_selection_derivable` | 2666+ | ⚠️ audit | May be superseded by standard fn with derivable weighting |

---

## Step 0: Function Audit (Before Any Deletion)

The refactor agent must audit these 4 functions before the split — do not move them to new files
until the audit determines their fate.

### Audit 1: `create_covariance_and_precision_matrix_both_directions` (line 1258)
**Question**: Is this a standalone utility or entirely superseded by calling
`create_covariance_and_precision_matrix(..., with_and_against = TRUE)`?
- Check: grep all call sites in `serious-testing/**` for `both_directions`
- Check: compare return value structure against the `with_and_against = TRUE` branch
- Decision: if no external call sites AND equivalent output → **delete**, replace any remnant calls
  with `with_and_against = TRUE` argument

### Audit 2: `expand_border` (line 1573) and `partition_by_county` (line 1593)
**Problem**: Both functions are defined twice — once here and once in `partition_norway.R` (lines 55, 76).
- Check: which version was written last (git blame or comment style)
- Check: are the two definitions identical? If different, which is correct?
- Decision: one definition becomes canonical (move to `R/utilities.R`). The duplicate in
  `partition_norway.R` is removed. `partition_norway.R` sources from the canonical location.

### Audit 3: `greedy_mi_sensor_selection_derivable` (line 2666)
**Question**: Is this an active variant or a superseded experiment?
- Check: grep all call sites — only in `testing-sensor-selection-derivable-norway.R`?
- Check: does `greedy_mi_sensor_selection_norway` with derivable-aware weighting already cover this?
- Decision: if no production use and functionality subsumed → **delete**, document removal reason

---

## Step 1: Proposed File Split

After audit, move functions into `R/` (standard R project convention):

### `R/inla-modeling.R`
- `fit_weighted_inla_model`
- `kfold_cv_inla`
- `select_similarity_covariates`

### `R/flow-conservation.R`
- `identify_derivable_nodes`

### `R/covariance-matrices.R`
- `create_covariance_and_precision_matrix`
- `create_covariance_and_precision_matrix_both_directions` (if kept after audit)

### `R/sensor-selection.R`
- `compute_weights`
- `greedy_mi_sensor_selection`
- `greedy_mi_sensor_selection_norway`
- `greedy_mi_sensor_selection_derivable` (if kept after audit)

### `R/diagnostics.R`
- `plot_inla_model`
- `plot_kfold_cv`
- `plot_sensor_selection_map`
- `map_traffic_links`

### `R/utilities.R`
- `expand_border` (canonical version, after audit)
- `partition_by_county` (canonical version, after audit)
- `categorise_aadt`
- `summarise_sensor_selection`
- `add_geometries`

### `serious-testing/scripts_inla_sensor_placement.R` (becomes thin orchestrator)
```r
# Source all sub-modules — replaces the monolithic script
source("R/inla-modeling.R")
source("R/flow-conservation.R")
source("R/covariance-matrices.R")
source("R/sensor-selection.R")
source("R/diagnostics.R")
source("R/utilities.R")
```
All 13 dependent scripts continue to `source("serious-testing/scripts_inla_sensor_placement.R")`
unchanged — zero breaking changes for them.

---

## Step 2: Testing Script Audit (Dead Script Removal)

These scripts in `serious-testing/` are candidates for archival or deletion. Each requires
**explicit confirmation** before removal.

| Script | Status | Reason |
|--------|--------|--------|
| `testing-norway.R` | ⚠️ likely obsolete | Early prototype; superseded by `testing-sensor-selection-norway.R`; confirm no unique logic |
| `testing-sensor-selection-buskerud.R` | ⚠️ likely obsolete | Regional experiment for Buskerud; check if findings were incorporated into the Norway pipeline |
| `partition_norway.R` | ⚠️ likely obsolete | Contains duplicate `expand_border`/`partition_by_county`; if canonical versions move to `R/utilities.R`, this script's purpose is gone |
| `testing-extraconstr-norway.R` | ⚠️ check | Extra constraints experiment; verify if the approach was rejected or merged into current pipeline |
| `covariance_analysis.R` | ⚠️ check | Covariance analysis; check if findings informed the current model or are still needed |
| `pt_lanes_adjacency_comparison.R` | ⚠️ check | PT lanes adjacency study; check if conclusion was incorporated and script is now a dead end |
| `testing-sensor-selection-filtered.R` | ⚠️ check | Filtered variant; unclear if active or superseded |
| `testing-sensor-selection-derivable-norway.R` | ⚠️ depends | Keep if `greedy_mi_sensor_selection_derivable` is kept; review together with Audit 3 above |
| `testing_inla_rbf_model.R` | ⚠️ check | RBF model testing; may contain useful diagnostics or be fully superseded by `run_cv_comparison.R` |
| `demonstration-traffic-links.R` | 🔵 keep | Demonstration/presentation script; has presentation value even if not in the pipeline |
| `testing-derivable-nodes.R` | 🔵 keep | Active exploration of derivable node logic; `identify_derivable_nodes` is actively developed |
| `run_cv_comparison.R` | ✅ keep | Active CV model comparison; produces model selection decisions |
| `testing-sensor-selection-norway.R` | ✅ keep | Primary pipeline script |

**Process**: The `/cleanup` prompt handles this step. It reads this plan, lists candidates, and
asks for per-item confirmation before any deletion.

---

## Step 3: Documentation (Post-Split)

After the split, invoke the documentation agent to:
1. Add roxygen2 `@title`, `@description`, `@param`, `@return` blocks to each function in `R/`
2. Update `/memories/repo/xdt-sensor-locations-structure.md` if any function signatures changed
3. Update `docs/functions/*.Rmd` narrative files to reflect the new file structure
4. Remove roxygen2 for any deleted functions

---

## Execution Order

```
1. Step 0: Function audit (refactor agent, read-only phase)
   → confirm audit decisions before proceeding

2. Step 1: File split (refactor agent, execute phase)
   → create R/ sub-files
   → convert scripts_inla_sensor_placement.R to source() orchestrator
   → verify all 13 dependent scripts still work (source the orchestrator, check no errors)

3. Step 2: Script cleanup (/cleanup prompt, one item at a time)
   → per-item confirmation for each candidate script

4. Step 3: Documentation (documentation agent, plan mode)
   → review proposed roxygen2 additions
   → approve per-file before writing
```

---

## What Does NOT Change

- All `source("serious-testing/scripts_inla_sensor_placement.R")` call sites — zero changes needed
- The `shiny/` app — it does not source this file
- `deploy.R` — does not source this file
- Result `.rds` file structure — purely internal refactor, no output changes
