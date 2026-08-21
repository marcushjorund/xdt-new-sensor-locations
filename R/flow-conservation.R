# ════════════════════════════════════════════════════════════════════════════════
# identify_derivable_nodes() ----
# ════════════════════════════════════════════════════════════════════════════════

#' Identify derivable traffic links via flow conservation
#'
#' For each flow node (row) in the incidence matrix, counts how many of its
#' nonzero-column links are unmeasured (\code{is.na(aadt)}).  A flow node is
#' *derivable* when exactly 2 of its links are unmeasured -- measuring one
#' allows the other to be derived from flow conservation.
#'
#' A link may participate in multiple derivable flow nodes, meaning measuring
#' it simultaneously unlocks derivation in all of them.
#'
#' @param incidence_matrix Sparse or dense incidence matrix from
#'   \code{build_incidence_matrix()}.  Rows are flow nodes (named), columns
#'   are traffic link IDs (named via \code{colnames}).
#' @param traffic_links Data frame with columns \code{id}, \code{aadt}, and
#'   optionally \code{traffic_volume_source}.
#'   Links with \code{is.na(aadt)} or whose source is not in
#'   \code{reliable_sources} are treated as unmeasured.
#' @param reliable_sources Character vector of \code{traffic_volume_source}
#'   values considered trustworthy enough to anchor flow conservation.
#'   Defaults to \code{c("Trafikkdata_continuous", "AutoPASS", "Derived")}.
#'   Set to \code{NULL} to fall back to \code{is.na(aadt)} only.
#'
#' @return \code{traffic_links} with three additional columns:
#'   \describe{
#'     \item{\code{derivable}}{Logical. \code{TRUE} if the link is one of the
#'       two unmeasured members of at least one derivable flow node.}
#'     \item{\code{n_derivable}}{Integer. Number of additional links that become
#'       derivable if this link is measured (0 for non-derivable links).}
#'     \item{\code{derivable_flow_nodes}}{List-column of character vectors.
#'       Names of all derivable flow nodes this link participates in.
#'       \code{NULL} for non-derivable links.  Group on this column (after
#'       \code{tidyr::unnest}) to find paired links.}
#'   }
identify_derivable_nodes <- function(
    incidence_matrix,
    traffic_links,
    nodes            = NULL,
    reliable_sources = c("Trafikkdata_continuous", "AutoPASS", "Derived")) {
  if (!"id" %in% names(traffic_links))
    stop("traffic_links must have an 'id' column.")
  if (!"aadt" %in% names(traffic_links))
    stop("traffic_links must have an 'aadt' column.")
  if (!is.null(nodes) && !"number_of_traffic_links" %in% names(nodes))
    stop("nodes must contain a 'number_of_traffic_links' column (produced by identify_unbalanceable_nodes()).")

  A <- as(incidence_matrix, "dgCMatrix")

  # A link is "unmeasured" if aadt is NA OR its source is not in reliable_sources.
  col_link_idx <- match(colnames(A), traffic_links$id)
  has_src      <- "traffic_volume_source" %in% names(traffic_links)
  is_unmeasured <- is.na(traffic_links$aadt[col_link_idx]) |
    (has_src & !traffic_links$traffic_volume_source[col_link_idx] %in% reliable_sources)

  # Sparse triplet: (row, col) of every nonzero entry
  trip       <- Matrix::summary(A)
  rows_split <- split(trip$j, trip$i)

  # Row integer indices → row names (Matrix::summary returns integer row indices)
  row_names <- rownames(A)

  # ── Build exactness flag per row ──────────────────────────────────────────────
  # A _passthrough row is approximate when the underlying node has exactly 2
  # known traffic links: no other incidence rows cover additional roads there,
  # so an unobserved diversion cannot be ruled out.
  # A _passthrough row at a node with >2 traffic links is exact: the other roads
  # are fully captured in separate _mixing rows at the same node.
  # All _mixing rows are exact by construction.
  # When nodes = NULL every row is treated as exact (backward compatible).
  if (!is.null(nodes)) {
    node_links_lookup <- stats::setNames(
      nodes$number_of_traffic_links,
      as.character(nodes$id)
    )
    row_node_ids    <- sub("_component_.*", "", row_names)
    is_passthrough  <- grepl("_passthrough", row_names)
    n_links_at_node <- node_links_lookup[row_node_ids]
    row_is_exact    <- !is_passthrough |
      (!is.na(n_links_at_node) & n_links_at_node > 2L)
  } else {
    row_is_exact <- rep(TRUE, length(row_names))
  }
  # Named by row name for safe lookup
  names(row_is_exact) <- row_names

  # ── 1. Currently derivable (== 1 unmeasured link in row) ─────────────────────
  # Flow conservation has exactly 1 unknown → uniquely determined from existing
  # measurements. No new sensor needed.
  derivable_rows <- Filter(function(cols) sum(is_unmeasured[cols]) == 1L, rows_split)

  link_to_nodes <- vector("list", ncol(A))
  for (row_idx in names(derivable_rows)) {
    cols     <- derivable_rows[[row_idx]]
    col      <- cols[is_unmeasured[cols]]   # exactly 1
    row_name <- row_names[as.integer(row_idx)]
    link_to_nodes[[col]] <- c(link_to_nodes[[col]], row_name)
  }

  traffic_links$derivable_flow_nodes <- vector("list", nrow(traffic_links))
  traffic_links$n_derivable          <- 0L
  traffic_links$derivable            <- FALSE

  for (col_idx in seq_along(link_to_nodes)) {
    fn <- link_to_nodes[[col_idx]]
    if (is.null(fn)) next
    link_row <- col_link_idx[col_idx]
    if (is.na(link_row)) next
    traffic_links$derivable_flow_nodes[[link_row]] <- fn
    traffic_links$n_derivable[link_row]            <- length(fn)
    traffic_links$derivable[link_row]              <- TRUE
  }

  # ── 2. Sensor-placement value (== 2 unmeasured links in row, exact rows only) ─
  # Measuring either of the two unmeasured links reduces the row to 1 unknown,
  # making the other derivable via flow conservation.
  # Excluded:
  #   - Rows where the two unmeasured links are WITH/AGAINST of the same road
  #     (terminus nodes — both directions measured simultaneously, nothing unlocked).
  #   - Approximate passthrough rows (node has only 2 known traffic links).
  same_road <- function(id1, id2) {
    sub("-WITH$|-AGAINST$", "", id1) == sub("-WITH$|-AGAINST$", "", id2)
  }

  enables_rows <- list()
  for (row_idx in names(rows_split)) {
    cols <- rows_split[[row_idx]]
    if (sum(is_unmeasured[cols]) != 2L) next
    ums      <- cols[is_unmeasured[cols]]
    row_name <- row_names[as.integer(row_idx)]
    if (same_road(colnames(A)[ums[1L]], colnames(A)[ums[2L]])) next
    if (!isTRUE(row_is_exact[row_name])) next
    enables_rows[[row_idx]] <- cols
  }

  link_enables <- vector("list", ncol(A))
  for (row_idx in names(enables_rows)) {
    cols <- enables_rows[[row_idx]]
    ums  <- cols[is_unmeasured[cols]]   # exactly 2: u1, u2
    u1 <- ums[1L]; u2 <- ums[2L]
    link_enables[[u1]] <- c(link_enables[[u1]], colnames(A)[u2])
    link_enables[[u2]] <- c(link_enables[[u2]], colnames(A)[u1])
  }

  traffic_links$enables_derivable_links <- vector("list", nrow(traffic_links))
  traffic_links$n_enables_derivable     <- 0L
  traffic_links$enables_derivable       <- FALSE

  for (col_idx in seq_along(link_enables)) {
    els <- unique(link_enables[[col_idx]])
    if (is.null(els)) next
    link_row <- col_link_idx[col_idx]
    if (is.na(link_row)) next
    traffic_links$enables_derivable_links[[link_row]] <- els
    traffic_links$n_enables_derivable[link_row]       <- length(els)
    traffic_links$enables_derivable[link_row]         <- TRUE
  }

  traffic_links
}


