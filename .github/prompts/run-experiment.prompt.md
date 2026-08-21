---
description: "Set up and run a new greedy MI sensor selection experiment"
agent: agent
argument-hint: "Describe the weighting scheme: FRC weights, AADT blend, neighbourhood decay"
---

Set up a new sensor selection experiment in the Norway pipeline. Load the `sensor-placement` skill before proceeding.

## Steps

1. **Specify parameters** — confirm the following with the user or infer from context:
   - `weighting_bias`: numeric vector indexed by FRC level (or `NULL` for unweighted)
   - `frc_aadt_alpha`: blend factor between 0 and 1 (1 = FRC only, 0 = AADT only)
   - `neighbourhood_decay`: tuple `c(f0, f1, f2)` for 0/1/2-hop smoothing (or `NULL`)
   - `k_per_county`: sensors per county budget
   - `r`: total sensors nationwide

2. **Define output filename** following the naming convention:
   ```
   <weighting>_<k>percounty_<r>overall.rds
   ```

3. **Add a new experiment block** to [production/run_norway_sensor_selection.R](../../production/run_norway_sensor_selection.R). Mirror the exact structure of existing `norway_sensors_*` blocks:
   - Call `greedy_mi_sensor_selection_norway()` with the new parameters
   - Assign result to a descriptively named variable
   - Add a `saveRDS(result, file.path("results", "<filename>.rds"))` call

4. **Add the result to the app** — after the experiment produces output, run `/add-config` with the new filename and a UI label.

## Constraints

- Do not remove or modify any existing experiment blocks.
- Follow R conventions: `|>` pipes, snake_case, no `%>%`.
