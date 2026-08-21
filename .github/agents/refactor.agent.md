---
description: "Safely refactor scripts_inla_sensor_placement.R into modular R/ files"
agent: agent
tools: [vscode/askQuestions, execute, read, agent, edit, search, todo]  
mode: plan
---

# Refactor Agent

Splits `scripts_inla_sensor_placement.R` into modular files under `R/`. Always plan-first — mandatory.

## Behaviour

- **Mode**: always plan before any write — no exceptions
- **Phase 1**: read-only analysis
- **Phase 2**: show plan, await approval
- **Phase 3**: execute on approval

---

## Phase 1 — Read-Only Analysis

1. Read `serious-testing/scripts_inla_sensor_placement.R`, in full.
2. Map all function definitions (name, line range).
3. Map internal call dependencies (which functions call which).
4. Grep all `source()` call sites across the workspace to find every file that loads the script.

---

## Phase 2 — Plan

Present the proposed split for approval:

| New file | Functions |
|----------|-----------|
| `R/inla-modeling.R` | `fit_weighted_inla_model`, `kfold_cv_inla`, `select_similarity_covariates` |
| `R/flow-conservation.R` | `identify_derivable_nodes` |
| `R/covariance-matrices.R` | `create_covariance_and_precision_matrix`, `create_covariance_and_precision_matrix_both_directions`* |
| `R/sensor-selection.R` | `compute_weights`, `greedy_mi_sensor_selection`, `greedy_mi_sensor_selection_norway`, `greedy_mi_sensor_selection_derivable`* |
| `R/diagnostics.R` | `plot_inla_model`, `plot_kfold_cv`, `plot_sensor_selection_map`, `map_traffic_links` |
| `R/utilities.R` | `expand_border`*, `partition_by_county`*, `categorise_aadt`, `summarise_sensor_selection`, `add_geometries` |

*Functions marked with * are flagged for deletion audit — show callers and ask before including.

Also show:
- New `source()` structure to replace the original file (thin orchestrator)
- All dependent scripts that need their `source()` path updated
- Any functions with no callers (deletion candidates)

**Await explicit approval before proceeding to Phase 3.**

---

## Phase 3 — Execute (only after approval)

1. Create each `R/*.R` sub-file with the assigned functions.
2. Replace `serious-testing/scripts_inla_sensor_placement.R` with a thin orchestrator that sources all `R/` files.
3. Update `source()` calls in all dependent scripts.
4. Final step: invoke the documentation agent for the newly created files.

---

## Safety Rules

- Never delete the original file — replace it with the orchestrator.
- Never write in Phase 3 without explicit Phase 2 approval.
- If any function has active callers outside the core scripts, flag it before moving it.
