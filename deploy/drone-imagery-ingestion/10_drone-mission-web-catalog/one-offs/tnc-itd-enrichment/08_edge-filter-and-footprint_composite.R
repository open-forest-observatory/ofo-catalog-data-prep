# Purpose: Remove edge-effect trees/crowns and compute analysis-ready forested footprints. For each
# composite, computes a convex hull around crowns, buffers inward by 1 m, removes crowns that
# intersect the hull boundary (and their corresponding tree points), then derives a forested mission
# footprint by buffering out from tree points by 25 m, unioning, and buffering back in by 25 m.

source("deploy/drone-imagery-ingestion/00_set-constants.R")
library(sf)

DELIVERABLES_DIR = "/ofo-share/project-data/tnc-yuba-deliverables"
ALL_DETECTED_CROWNS_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/all-detected-tree-crowns.gpkg")
ALL_DETECTED_TREES_WITH_DBH_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/all-detected-trees-with-dbh.gpkg")

ANALYSIS_READY_TREES_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/all-detected-trees_analysis-ready.gpkg")
ANALYSIS_READY_CROWNS_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/all-detected-crowns_analysis-ready.gpkg")
ANALYSIS_READY_FOOTPRINTS_FILEPATH = file.path(DELIVERABLES_DIR, "composite-missions/overall/forested-mission-footprints_analysis-ready.gpkg")

EDGE_BUFFER_M = 1
FOOTPRINT_SIMPLIFY_TOLERANCE_M = 5
FOOTPRINT_SMOOTH_BUFFER_M = 3
FOOTPRINT_BUFFER_M = 25


# --- Load data ---

all_crowns = st_read(ALL_DETECTED_CROWNS_FILEPATH, quiet = TRUE)
all_trees = st_read(ALL_DETECTED_TREES_WITH_DBH_FILEPATH, quiet = TRUE)


# --- Remove trees without a crown ---

n_before = nrow(all_trees)
all_trees = all_trees |> dplyr::filter(!is.na(crown_area_sqm))
cat("Removed", n_before - nrow(all_trees), "trees with no crown (NA crown area);", nrow(all_trees), "remain\n")


# --- Process each composite: edge filtering and footprint ---

composite_ids = unique(all_trees$composite_id)

filtered_trees_list = list()
filtered_crowns_list = list()
footprints_list = list()

for (i in seq_along(composite_ids)) {
  cid = composite_ids[i]
  cat("  Processing", cid, "(", i, "/", length(composite_ids), ")\n")

  trees_foc = all_trees |> dplyr::filter(composite_id == cid)
  crowns_foc = all_crowns |> dplyr::filter(composite_id == cid)

  # Keep only crowns that have a matching tree point (after NA crown area removal)
  crowns_foc = crowns_foc |>
    dplyr::filter(treetop_unique_ID %in% trees_foc$unique_ID)

  if (nrow(crowns_foc) < 3) {
    cat("    Fewer than 3 crowns, skipping edge filtering\n")
    filtered_trees_list[[cid]] = trees_foc
    filtered_crowns_list[[cid]] = crowns_foc
    next
  }

  ## Need to compute forested area in a way that does not bias density estimate upward by clipping at the edge of the crowns (because that would bias for higher density because the plot would be designed to include edge trees)

  # Compute convex hull around crowns, buffer inward
  crowns_union = st_union(crowns_foc)
  hull = st_convex_hull(crowns_union)
  hull_buffered = st_buffer(hull, -EDGE_BUFFER_M)

  # Remove crowns that intersect the boundary of the buffered hull (i.e. not fully contained)
  crowns_contained = st_within(crowns_foc, hull_buffered, sparse = FALSE)[, 1]
  n_edge_crowns = sum(!crowns_contained)

  crowns_kept = crowns_foc[crowns_contained, ]
  trees_kept = trees_foc |>
    dplyr::filter(unique_ID %in% crowns_kept$treetop_unique_ID)

  cat("    Convex hull edge filter: removed", n_edge_crowns,
      "crowns (and trees) not fully within hull;",
      nrow(crowns_kept), "crowns and", nrow(trees_kept), "trees remain\n")

  # Compute forested mission footprint: buffer out from tree points, union, buffer back in,
  # then simplify and smooth to get a clean boundary
  if (nrow(trees_kept) > 0) {
    trees_buffered = st_buffer(trees_kept, FOOTPRINT_BUFFER_M)
    footprint = st_union(trees_buffered)
    footprint = st_buffer(footprint, -FOOTPRINT_BUFFER_M)
    footprint = st_simplify(footprint, dTolerance = FOOTPRINT_SIMPLIFY_TOLERANCE_M)
    footprint = st_buffer(footprint, -FOOTPRINT_SMOOTH_BUFFER_M)
    footprint = st_simplify(footprint, dTolerance = FOOTPRINT_SIMPLIFY_TOLERANCE_M)
    footprint_sf = st_sf(composite_id = cid, geometry = footprint)
    footprints_list[[cid]] = footprint_sf

    # Remove trees and crowns that fall outside the final footprint
    inside = st_within(trees_kept, footprint_sf, sparse = FALSE)[, 1]
    n_outside = sum(!inside)
    trees_kept = trees_kept[inside, ]
    crowns_kept = crowns_kept |>
      dplyr::filter(treetop_unique_ID %in% trees_kept$unique_ID)
    if (n_outside > 0) {
      cat("    Footprint clip: removed", n_outside,
          "trees (and crowns) outside smoothed footprint;",
          nrow(trees_kept), "trees and", nrow(crowns_kept),
          "crowns remain\n")
    }
  }

  filtered_trees_list[[cid]] = trees_kept
  filtered_crowns_list[[cid]] = crowns_kept
}


# --- Combine and save ---

analysis_ready_trees = dplyr::bind_rows(filtered_trees_list)
analysis_ready_crowns = dplyr::bind_rows(filtered_crowns_list)
analysis_ready_footprints = dplyr::bind_rows(footprints_list)

# Temporary fixup until re-run upstream: make numeric cols numeric again
analysis_ready_trees = analysis_ready_trees |>
  mutate(across(c(height, ends_with("_frac_matching_mode"), ends_with("_n_preds")), as.numeric))


cat("\nFinal counts:\n")
cat("  Trees:", nrow(analysis_ready_trees), "\n")
cat("  Crowns:", nrow(analysis_ready_crowns), "\n")
cat("  Footprints:", nrow(analysis_ready_footprints), "\n")

st_write(analysis_ready_trees, ANALYSIS_READY_TREES_FILEPATH, delete_dsn = TRUE)
st_write(analysis_ready_crowns, ANALYSIS_READY_CROWNS_FILEPATH, delete_dsn = TRUE)
st_write(analysis_ready_footprints, ANALYSIS_READY_FOOTPRINTS_FILEPATH, delete_dsn = TRUE)

cat("Saved analysis-ready files\n")
