---
name: sensor-placement
description: "Use when: greedy MI sensor selection; mutual information criterion; covariance matrix;
  precision matrix; Schur complement update; flow conservation; derivable nodes; incidence matrix;
  WITH/AGAINST direction pairing; INLA RBF hyperparameters (tau d sigma); Krause et al 2008;
  greedy_mi_sensor_selection(); greedy_mi_sensor_selection_norway(); identify_derivable_nodes();
  create_covariance_and_precision_matrix(); weighting_bias; neighbourhood_decay; frc_aadt_alpha"
argument-hint: "Describe the sensor selection problem or algorithm step"
---

# Sensor Placement Skill

## Domain Overview

Norway road network. The goal is to allocate a budget of k sensors per county plus r sensors
nationwide to unmeasured road links. Physical constraint: one sensor always measures both the
WITH and AGAINST direction simultaneously (bidirectional sensor). The dataset splits each road
into two directional links (`road-id-WITH` and `road-id-AGAINST`), so the covariance matrix
must be transformed from link-level to road-level before selection.

---

## Two-Stage Pipeline

```
fit_weighted_inla_model()
  └─► identify_derivable_nodes()
        └─► create_covariance_and_precision_matrix()
              └─► greedy_mi_sensor_selection_norway()
```

1. **`fit_weighted_inla_model`** — Besag-proper weighted INLA spatial model. Returns posterior predictions
   on all links plus hyperparameters: `tau`, `d`, `sigma`, `distances`.
2. **`identify_derivable_nodes`** — marks which unmeasured links can be derived from flow
   conservation (incidence matrix A); adds `derivable`, `n_derivable`, `enables_derivable` columns.
3. **`create_covariance_and_precision_matrix`** — builds Σ from INLA hyperparameters; when
   `with_and_against=TRUE` applies C·Σ·Cᵀ transformation to get road-level covariances.
4. **`greedy_mi_sensor_selection_norway`** — county-constrained greedy selection; calls
   `greedy_mi_sensor_selection` per county then globally.

---

## Key Function Signatures

```r
identify_derivable_nodes(
  data,             # traffic links data frame with aadt column
  nodes,            # traffic nodes sf object (NULL = backward-compatible mode)
  adjacency_matrix  # sparse n×n undirected adjacency
)

create_covariance_and_precision_matrix(
  adjacency_matrix,           # sparse n×n
  tau, d, sigma,              # INLA RBF hyperparameters
  distances,                  # edge distance vector
  data          = NULL,       # traffic links (needed for with_and_against)
  with_and_against = FALSE    # apply C·Σ·Cᵀ road-level transformation
)

greedy_mi_sensor_selection(
  data,                       # traffic links with aadt (NA = unmeasured)
  covariance_matrix,          # m×m (m = # base road IDs)
  ids,                        # names for rows/cols of covariance_matrix
  k,                          # sensors to select
  adjacency_matrix = NULL,
  weighting_bias = NULL,      # FRC-based weight vector (length = # FRC levels)
  weighting_bias_aadt = FALSE,# use inla_pred column for AADT weighting
  frc_aadt_alpha = 1,         # blend α: w = w_frc^α · w_aadt^(1-α)
  neighbourhood_decay = NULL, # tuple (f0, f1, f2) for 0/1/2-hop smoothing
  append_neighbours = FALSE   # add 1-hop neighbours to output
)

greedy_mi_sensor_selection_norway(
  data,
  covariance_matrix,
  ids,
  k_per_county,               # sensors per county (budget)
  r,                          # total sensors nationwide
  county_column = "county",
  ...                         # passed to greedy_mi_sensor_selection
)
```

---

## WITH/AGAINST Direction Pairing (C Matrix)

Traffic sensors measure both directions simultaneously. The dataset has directional links
(`road-id-WITH`, `road-id-AGAINST`). The C matrix (m×n sparse) maps n links → m roads:

```
C is m × n where:
  m = length(paired_bases) + length(unpaired_bases)
  n = total number of directional links

For paired roads (both WITH and AGAINST present):
  C[i, WITH_idx]    = 1
  C[i, AGAINST_idx] = 1

For unpaired roads (WITH only):
  C[i, WITH_idx] = 1
```

Road-level covariance: **Σ_road = C · Σ_link · Cᵀ**

Variance of a bidirectional sensor at road i:
`Var(WITH + AGAINST) = Var(WITH) + Var(AGAINST) + 2·Cov(WITH, AGAINST)`

---

## Derivable Nodes Logic

Incidence matrix A (flow conservation equations):
- Each row = a flow node (weighted sum of flows through that node = 0)
- Each column = a traffic link

**Three row categories:**
| Row type | Unknowns | Action |
|----------|----------|--------|
| Derivable | exactly 1 unmeasured link | that link is already derivable; no sensor needed |
| Sensor-value (exact) | exactly 2 unknowns, exact row | measuring one makes the other derivable |
| Approximate | passthrough, node with =2 links | excluded — unobserved diversions possible |

**Exactness rules:**
- Passthrough at node with >2 links → EXACT (other roads in separate rows)
- Passthrough at node with =2 links → APPROXIMATE
- All mixing rows → EXACT

**Output columns added to data:**
- `derivable`: TRUE if link is unmeasured and in a derivable row
- `n_derivable`: count of additional links that become derivable if this link is measured
- `enables_derivable`: TRUE if measuring this link enables derivations
- `enables_derivable_links`: list of link IDs that become derivable

---

## MI Scoring Formula

At each greedy step, the score for candidate location s is:

$$\text{score}[s] = \frac{\text{cond\_var1}[s]}{\text{cond\_var2}[s]} \times w[s]$$

- **`cond_var1[s]`** — marginal variance of s given already-measured S₀ (how uncertain s is marginally)
- **`cond_var2[s]`** — variance of s given all other unmeasured locations (how much uniquely known)
- **`w[s]`** — importance weight

Higher ratio = larger MI gain from placing a sensor at s.

---

## Weighting System

**FRC weights** (`weighting_bias` vector, indexed by `functionalRoadClass + 1`):
```r
w_frc[s] = weighting_bias[functionalRoadClass[s] + 1]
```

**AADT weights** (`weighting_bias_aadt = TRUE`, requires `inla_pred` column):
```r
w_aadt[s] = inla_pred[s] / median(inla_pred)
```

**Blend** (`frc_aadt_alpha`):
```r
w[s] = w_frc[s]^alpha * w_aadt[s]^(1 - alpha)
```

**Neighbourhood smoothing** (`neighbourhood_decay = c(f0, f1, f2)`):
```r
w_smooth[s] = f0*w[s] + f1*mean(w[neighbours_1hop(s)]) + f2*mean(w[neighbours_2hop(s)])
```
Only unmeasured neighbours are included.

---

## Result `.rds` Structure

```r
list(
  selected_data_entries = <sf data frame>,  # WITH/AGAINST rows + optional neighbours
  summary = list(
    n_selected   = <integer>,
    k_per_county = <integer>,
    r            = <integer>,
    counts_table = <data frame>             # per-county sensor counts
  )
)
```

`selected_data_entries` columns: all original data columns + `selected` (logical) + `mi_score` (numeric).

---

## Output Naming Convention

```
<weighting>_<k>percounty_<r>overall.rds
```

Examples:
- `frc_sqrt(0.9_to_0)_100percounty_100overall.rds`
- `frc_sqrt(0.9_to_0)_and_logAADT_0.9_to_0.1_weight_100percounty_100overall.rds`
- `unweighted_100percounty100overall.rds`

---

## Precision Update (Schur Complement)

After selecting sensor at index j, the precision matrix for remaining unmeasured shrinks:

$$\text{prec}_\text{new} = \text{prec}_{-j,-j} - \frac{\text{prec}_{-j,j} \cdot \text{prec}_{j,-j}}{\text{prec}_{j,j}}$$

Avoids re-computing Cholesky each iteration. `cond_var2` updated as `1/diag(prec_new)`.

## Incremental Cholesky (cond_var1 update)

Maintains a factor L₂ (m×k) across greedy iterations:

$$\text{cond\_var1}^{(i)} = \text{cond\_var1}^{(i-1)} - L_2[\cdot,i]^2$$

where `L2[•,i] = col_vec / sqrt(pivot)` and `col_vec = Σ[•, best] - AᵀA[•, best]`.
