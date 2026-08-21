# AI Infrastructure Plan — xdt-new-sensor-locations

## Overview

Build a complete AI infrastructure layer on top of the existing r-data-science and shiny-r skills.
13 new/modified files across 6 phases.

---

## Current State

| Status | File | Purpose |
|--------|------|---------|
| ✅ exists | `.github/skills/r-data-science/SKILL.md` | R/stats expert skill |
| ✅ exists | `.github/skills/shiny-r/SKILL.md` | Shiny expert skill |
| ✅ exists | `.github/instructions/shiny.instructions.md` | Shiny file instruction (minimal) |
| ✅ exists | `/memories/repo/xdt-sensor-locations-structure.md` | Function signatures reference |
| ✅ exists | `/memories/repo/greedy-mi-sensor-selection-algorithm.md` | Algorithm deep-dive |
| ❌ missing | `.github/copilot-instructions.md` | Always-on project context |
| ❌ missing | `.github/instructions/sensor-selection.instructions.md` | Core scripts instruction |
| ❌ missing | `.github/skills/sensor-placement/SKILL.md` | Algorithm domain skill |
| ❌ missing | `.github/prompts/*.prompt.md` | Repeatable task prompts |
| ❌ missing | `.github/agents/*.agent.md` | Custom agents |
| ❌ missing | `.github/hooks/doc-queue.json` | PostToolUse hook |
| ❌ missing | `docs/` | Narrative documentation folder |

---

## Phase 1 — Global Instructions

### Step 1: Create `.github/copilot-instructions.md`

Always loaded for every chat interaction in this workspace.

Content:
- **Project**: Norway AADT sensor placement optimisation. Goal: select k new traffic sensor locations
  that maximally reduce prediction uncertainty on unmeasured links. Two-stage pipeline:
  (1) INLA Besag-proper RBF spatial model → posterior predictions, (2) greedy MI sensor selection.
- **Key files**:
  - `serious-testing/scripts_inla_sensor_placement.R` — core function library (19 functions, 2700+ lines)
  - `serious-testing/testing-sensor-selection-norway.R` — main pipeline script (run this to produce results)
  - `results/` — output .rds files; each contains `selected_data_entries` + `summary` list
  - `shiny/app.R` — interactive visualisation of results
  - `deploy.R` — staging script that copies files into `shiny/` before rsconnect upload
- **Critical rule (deploy)**: When adding a new result file to `configs` in `app.R`, you MUST also
  add the filename to `rds_files` in `deploy.R`. Omitting this causes the deployed app to crash on
  startup with exit status 1.
- **R conventions**: Use `|>` (not `%>%`). snake_case. tidyverse-first. `req()` / `validate(need())`
  in Shiny reactive contexts.
- **Skills**: For statistical modelling → r-data-science skill. For Shiny/deploy → shiny-r skill.
  For greedy MI / covariance / flow conservation → sensor-placement skill.

---

## Phase 2 — File-Level Instructions

### Step 2: Create `.github/instructions/sensor-selection.instructions.md`

```yaml
---
applyTo: "serious-testing/**"
---
```

Content:
- Load the `r-data-science` skill when writing or modifying INLA model code, CV, or covariate selection.
- Load the `sensor-placement` skill when working with covariance matrices, derivable nodes, flow
  conservation, greedy MI selection, or WITH/AGAINST direction pairing.
- Output naming convention for new results: `<weighting>_<k>percounty_<r>overall.rds`
- After producing a new .rds in `results/`, always do both: add to `configs` in `shiny/app.R`
  AND add to `rds_files` in `deploy.R`.
- Function signatures are documented in `/memories/repo/xdt-sensor-locations-structure.md`.

### Step 3: Strengthen `.github/instructions/shiny.instructions.md`

Keep `applyTo: "shiny/**,deploy.R"`. Add inline rule (not just a pointer to the skill):

> **Deploy rule**: Never add an entry to `configs` in `app.R` without simultaneously adding the
> same filename to the `rds_files` list in `deploy.R`. These two edits must always happen together.
> Files missing from `deploy.R` are excluded from the rsconnect bundle and crash the app on startup.

---

## Phase 3 — New Domain Skill

### Step 4: Create `.github/skills/sensor-placement/SKILL.md`

Frontmatter:
```yaml
---
name: sensor-placement
description: "Use when: greedy MI sensor selection; mutual information criterion; covariance matrix;
  precision matrix; Schur complement update; flow conservation; derivable nodes; incidence matrix;
  WITH/AGAINST direction pairing; INLA RBF hyperparameters (tau d sigma); Krause et al 2008;
  greedy_mi_sensor_selection(); greedy_mi_sensor_selection_norway(); identify_derivable_nodes();
  create_covariance_and_precision_matrix(); weighting_bias; neighbourhood_decay; frc_aadt_alpha"
argument-hint: "Describe the sensor selection problem or algorithm step"
---
```

Content sections:
1. **Domain overview**: Norway road network, budget of k sensors per county + r nationwide, physical
   constraint that one sensor measures both WITH and AGAINST directions simultaneously.
2. **Two-stage pipeline**: `fit_weighted_inla_model` → `identify_derivable_nodes` →
   `create_covariance_and_precision_matrix` → `greedy_mi_sensor_selection_norway`
3. **Key function signatures**: pulled from `/memories/repo/xdt-sensor-locations-structure.md`
4. **WITH/AGAINST pattern**: C matrix (m×n sparse), Σ_road = C·Σ_link·Cᵀ
5. **Derivable nodes logic**: incidence matrix rows, 1-unknown → derivable, 2-unknown → enables_derivable
6. **MI scoring formula**: score ∝ cond_var1[s] / cond_var2[s] × w[s]
7. **Weighting system**: FRC weights, AADT weights, blend factor alpha, neighbourhood_decay tuple (f0,f1,f2)
8. **Result .rds structure**: `list(selected_data_entries, summary = list(n_selected, k_per_county, r, counts_table))`
9. **Output naming convention**: `<weighting>_<k>percounty_<r>overall.rds`

---

## Phase 4 — Prompt Files

### Step 5: Create `.github/prompts/add-config.prompt.md`

```yaml
---
description: "Add a new sensor selection result file to the Shiny app and deploy script"
agent: agent
argument-hint: "Filename in results/ and display label for the UI dropdown"
---
```

Task (atomic — all 3 must succeed or none):
1. Verify the named file exists in `results/`
2. Add an entry to the `configs` named list in [shiny/app.R](../../shiny/app.R)
3. Add the filename to `rds_files` in [deploy.R](../../deploy.R)

### Step 6: Create `.github/prompts/run-experiment.prompt.md`

```yaml
---
description: "Set up and run a new greedy MI sensor selection experiment"
agent: agent
argument-hint: "Describe the weighting scheme: FRC weights, AADT blend, neighbourhood decay"
---
```

Task:
1. Specify `weighting_bias`, `frc_aadt_alpha`, `neighbourhood_decay` values
2. Add a new call block to [testing-sensor-selection-norway.R](../../serious-testing/testing-sensor-selection-norway.R)
   mirroring the pattern of existing `norway_sensors_*` blocks
3. Define output filename following the convention `<weighting>_<k>percounty_<r>overall.rds`
4. Add `saveRDS(...)` call to write to `results/`

### Step 7: Create `.github/prompts/cleanup.prompt.md`

```yaml
---
description: "Safely remove dead code: unused functions or obsolete testing scripts"
agent: agent
---
```

Task:
1. Read the refactor plan at [.github/plans/refactor.md](./refactor.md) for the pre-identified candidates
2. For each candidate file or function, show what would be deleted and ask for explicit confirmation
3. Never delete anything without per-item approval
4. After deletion, check for any remaining references (grep) and flag broken call sites

---

## Phase 5 — Agents

### Step 8: Create `.github/agents/shiny-deploy.agent.md`

- Mode: plan (always shows proposed changes before writing)
- Tool restrictions: block `rm`, block `git push`
- Workflow: verify results/ file → update app.R configs → update deploy.R rds_files → confirm plan → execute

### Step 9: Create `.github/agents/experiment-designer.agent.md`

- Tools: read-only only (no write_file)
- Reads all `.rds` files in `results/`, inspects `summary$counts_table` and `mi_scores` distributions
- Compares weighting schemes, identifies gaps in coverage (counties, FRC levels, road categories)
- Output: structured recommendation for the next experiment → user hands off to `/run-experiment`

### Step 10: Create `.github/agents/refactor.agent.md`

- Mode: always plan first (mandatory)
- Phase 1 (read-only): read `scripts_inla_sensor_placement.R` fully; map all 19 function definitions
  and their internal call dependencies; grep all 13 source() call sites
- Phase 2 (plan): propose exact file split; list functions flagged for deletion audit;
  show new source() structure; await approval
- Phase 3 (execute): create sub-files under `R/`; replace original with thin source() orchestrator;
  update all dependent scripts
- Final step: invoke documentation agent for newly created files

Proposed split:
| New file | Functions |
|----------|-----------|
| `R/inla-modeling.R` | `fit_weighted_inla_model`, `kfold_cv_inla`, `select_similarity_covariates` |
| `R/flow-conservation.R` | `identify_derivable_nodes` |
| `R/covariance-matrices.R` | `create_covariance_and_precision_matrix`, `create_covariance_and_precision_matrix_both_directions`* |
| `R/sensor-selection.R` | `compute_weights`, `greedy_mi_sensor_selection`, `greedy_mi_sensor_selection_norway`, `greedy_mi_sensor_selection_derivable`* |
| `R/diagnostics.R` | `plot_inla_model`, `plot_kfold_cv`, `plot_sensor_selection_map`, `map_traffic_links` |
| `R/utilities.R` | `expand_border`*, `partition_by_county`*, `categorise_aadt`, `summarise_sensor_selection`, `add_geometries` |

*Flagged for deletion audit — see refactor plan.

### Step 11: Create `.github/agents/documentation.agent.md`

- Mode: always plan (never writes without explicit approval)
- Three-layer scope:
  1. Roxygen2 `@param` / `@return` / `@description` blocks in all `.R` files
  2. `/memories/repo/` files — flag if function signatures changed since last update
  3. `docs/*.Rmd` narrative files — flag sections that describe changed functions
- On startup: read `.github/doc-review-queue.txt` to know which files were modified
- Plan output: per-file diff of proposed documentation changes → await approval → clear queue
- Can be invoked manually or as final step by refactor agent

---

## Phase 6 — Hook + docs/ Scaffold

### Step 12: Create `.github/hooks/doc-queue.json`

PostToolUse hook that fires after every `write_file` call on `*.R` files.
Shell command: `echo "$FILE" >> .github/doc-review-queue.txt`
(deterministic — just records what changed; AI judgment happens inside documentation agent)

### Step 13: Create `docs/` seed files

- `docs/pipeline-overview.Rmd` — narrative walkthrough of the full two-stage pipeline
- `docs/functions/inla-modeling.Rmd` — reference for INLA fitting functions
- `docs/functions/sensor-selection.Rmd` — reference for greedy MI functions
- `docs/functions/flow-conservation.Rmd` — reference for derivable nodes
- `docs/functions/utilities.Rmd` — reference for mapping and post-processing utilities

---

## Verification Checklist

After implementing, verify:
1. Open `shiny/app.R` → deploy rule appears inline in shiny.instructions.md (not just deferred to skill)
2. Open `serious-testing/testing-sensor-selection-norway.R` → sensor-selection instruction active
3. Type `/add-config` in chat → 3-step atomic prompt is shown
4. Type `/run-experiment` in chat → parameterized template shown
5. Fresh chat, ask "what is this project?" → describes Norway AADT pipeline without searching
6. Write to any `.R` file → `.github/doc-review-queue.txt` is appended to
7. Invoke documentation agent → reads queue and proposes only what changed since last review
