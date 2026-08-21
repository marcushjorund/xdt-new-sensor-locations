# Source all sub-modules — replaces the monolithic script.
# Split into R/ on 2026-08-20 (refactor agent). See .github/plans/refactor.md.
source("R/inla-modeling.R")
source("R/flow-conservation.R")
source("R/covariance-matrices.R")
source("R/sensor-selection.R")
source("R/diagnostics.R")
source("R/utilities.R")
