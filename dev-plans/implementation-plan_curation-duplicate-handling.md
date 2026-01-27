# Implementation Plan: Duplicate Image Handling in Curation Workflow

## Background and Problem Statement

### Current Workflow
1. **05_parse-exif-metadata-per-image.R** - Parses EXIF metadata
2. **05b_find-duplicate-images.R** - Finds duplicate images (same lat/lon + file_size) per project
3. **05c_remove-duplicate-images.R** - Removes duplicates, keeping the one with alphabetically first `image_path_contrib`
4. **01_apply-curation-filters.R** - Applies curation notes to filter extraneous images

### Two Phases of Curation (Why Crosswalks Are Needed)

Curation of drone imagery has occurred in two phases:

1. **Formerly Curated Missions**: Some missions were curated using an earlier version of the image metadata (before image metadata processing was re-run). The curation notes for these missions reference **former image IDs** that no longer exist in the current image metadata. These missions are listed in `FORMERLY_CURATED_MISSION_LIST` and their old metadata is preserved in `FORMERLY_CURATED_IMAGE_METADATA`.

2. **Currently Curated Missions**: Other missions were curated after image metadata was regenerated. The curation notes for these missions reference **current image IDs** that match the existing metadata directly.

**The crosswalk** is needed only for formerly curated missions. It maps former image IDs to current image IDs by matching images based on geographic coordinates and attributes (file size, camera settings, etc.). This allows curation notes written with old IDs to be applied to the current metadata.

### Problem 1: Duplicate Images and Curation

Curators viewed imagery **before** duplicate removal. When they marked an image for removal, they may have marked only one of a duplicate pair (the "top" one visible in the curation interface). After duplicate removal, the image the curator marked may have been:
- **The kept image** - In this case, removing it by ID works correctly
- **The removed image** - In this case, the curation filter can't find it, and its duplicate remains

This affects both formerly curated and currently curated missions, since duplicate removal happens after curation viewing in both cases.

### Problem 2: Crosswalk Limitations

The crosswalk for formerly-curated missions currently records `"different sub-mission"` as a placeholder instead of the actual image ID when a match is found in a different sub-mission. This prevents:
- Looking up the image in the duplicate log
- Properly resolving which duplicate was kept
- Correct handling in the extraneous images column vs anomaly columns

### Key Distinction: Within-Mission vs Cross-Mission Duplicates

The duplicate detection script (05b) finds two types of duplicates:
- **Within-mission duplicates**: Same image appearing multiple times within a single mission (same `mission_id`, same coordinates, same file size)
- **Cross-mission duplicates**: Same image appearing across different missions (different `mission_id`, same coordinates, same file size)

**For this implementation, we only care about within-mission duplicates.** Only within-mission duplicates are removed by 05c_remove-duplicate-images.R. Cross-mission duplicates are logged for reference but are NOT removed and are not relevant to curation or crosswalking.

### Goals
1. Create a catalog-wide merged duplicate log from per-project logs
2. Update crosswalk creation to record actual image IDs (not `"different sub-mission"`)
3. When applying curation filters, account for duplicate pairs so that marking either duplicate for removal removes the remaining one
4. Handle crosswalk targets that were duplicates (use the remaining duplicate)
5. Expand image ranges before crosswalk application (since consecutiveness isn't preserved)
6. Mark deleted images appropriately in non-extraneous columns

---

## Implementation Steps

### Step 1: Add New Constants to `00_set-constants.R`

**File:** [00_set-constants.R](deploy/drone-imagery-ingestion/00_set-constants.R)

Add the following constants for the merged duplicate log:

```r
# Catalog-wide merged duplicate images log (combined from all per-project logs)
CATALOG_DUPLICATE_IMAGES_LOG_PATH = file.path(RAW_IMAGERY_INGESTION_PATH, "metadata-outputs/3_pre-curation-final/catalog-duplicate-images-log.csv")
```

**Location:** Add near line 68 (after `FULL_METADATA_PER_IMAGE_COMBINED_FILEPATH`)

---

### Step 2: Create New Script to Merge Per-Project Duplicate Logs

**File:** `deploy/drone-imagery-ingestion/01a_finalize-raw-imagery-metadata-prep/merge-duplicate-logs.R` (new file)

**Purpose:** Combine all per-project `*_duplicate-images_within-mission.csv` files into a single catalog-wide log.

**Implementation Details:**

1. **List all within-mission duplicate files** in `DUPLICATE_IMAGES_PATH` matching pattern `*_duplicate-images_within-mission.csv`

2. **Read and combine** all files using `map_dfr()` with `read_csv()`

3. **Reassign `duplicate_group_id`** since each project file has its own group IDs starting from 1:
   - Group by `mission_id`, `lat`, `lon`, `file_size_gb` across the entire combined dataset
   - Assign new sequential `duplicate_group_id` values
   - Keep a column `project_name` indicating the source project

4. **Write output** to `CATALOG_DUPLICATE_IMAGES_LOG_PATH`

**Key columns to preserve:**
- `duplicate_group_id` (reassigned)
- `project_name` (new - extracted from source filename)
- `image_id`
- `mission_id`
- `sub_mission_id`
- `lat`, `lon`, `file_size_gb`
- `image_path_contrib` (needed to determine which was kept)

---

### Step 3: Update Control Script to Include New Merge Script

The merge script is located in `01a_finalize-raw-imagery-metadata-prep/` which runs after all per-project processing is complete. It should be called from the appropriate control script in that directory, or run manually after all projects have been processed through the raw imagery metadata pipeline.

---

### Step 4: Update `create_image_id_crosswalk()` in `curation-utils.R`

**File:** [curation-utils.R](src/curation-utils.R)

**Current behavior (lines 276-280):** When exactly one match is found at a location but it's in a different sub-mission, the function returns `"different sub-mission"` instead of the actual image ID.

**Required change:** Always return the actual `current_image_id`, regardless of sub-mission.

**Modified code section (replace lines 266-281):**

```r
  # Handle single matches - now always return actual image ID
  single_match = match_counts == 1
  if (any(single_match)) {
    single_idx = which(single_match)
    single_current_ids = map_chr(matches[single_idx], ~ .x[1])
    crosswalk$current_image_id[single_idx] = single_current_ids
  }
```

**Reasoning:** Downstream logic can compute whether an image is in a different sub-mission by comparing sub-mission prefixes. Recording the actual ID enables:
- Looking up the image in the duplicate log
- Finding the correct remaining duplicate if needed
- Proper filtering in the extraneous images column

---

### Step 5: Create Helper Function for Duplicate Lookup

**File:** [curation-utils.R](src/curation-utils.R)

**Add new function** `get_remaining_duplicate()`:

```r
#' Find the remaining duplicate for an image that was removed during deduplication
#'
#' Given an image ID that may have been removed as a duplicate, look up its
#' duplicate group and return the image ID that was kept (alphabetically first
#' by image_path_contrib).
#'
#' @param image_id The image ID to look up
#' @param duplicate_log Tibble containing the catalog-wide duplicate log with columns:
#'   duplicate_group_id, image_id, image_path_contrib
#' @return Character: the kept image ID if found in a duplicate group,
#'   NA if not a duplicate, or the original image_id if it was the one kept
get_remaining_duplicate = function(image_id, duplicate_log) {
  # Find the duplicate group containing this image
  image_record = duplicate_log |>
    filter(image_id == !!image_id)

  if (nrow(image_record) == 0) {
    return(NA_character_)
  }

  group = image_record$duplicate_group_id[1]

  # Get all images in this group and find the one that was kept
  # (alphabetically first by image_path_contrib)
  group_images = duplicate_log |>
    filter(duplicate_group_id == group) |>
    arrange(image_path_contrib)

  # The first one was kept
  kept_id = group_images$image_id[1]
  return(kept_id)
}
```

**Note:** The existing `parse_extraneous_images()` function already handles range expansion, so no additional helper is needed. Use `parse_extraneous_images()` directly wherever range expansion is required.

---

### Step 6: Create `update_anomaly_image_references()` in `curation-utils.R`

**File:** [curation-utils.R](src/curation-utils.R)

**Purpose:** Update image ID references in anomaly columns (not extraneous_images) using crosswalk and duplicate resolution. Anomaly columns contain individual image IDs (not ranges), so no range expansion is needed.

**Note:** Extraneous images are handled separately in Step 9, where `parse_extraneous_images()` handles range expansion and then crosswalk + duplicate resolution are applied.

**Add new function** `update_anomaly_image_references()`:

```r
#' Update image ID references in anomaly columns
#'
#' Updates image IDs in anomaly columns using the crosswalk (if provided) and
#' duplicate resolution. This function is for anomaly columns only (not extraneous_images).
#' Anomaly columns contain individual image IDs, not ranges.
#'
#' For formerly curated missions: applies crosswalk AND duplicate resolution.
#' For currently curated missions: applies only duplicate resolution (crosswalk = NULL).
#'
#' @param text String containing image ID references (individual IDs, not ranges)
#' @param crosswalk Tibble with former_image_id and current_image_id columns,
#'   or NULL for currently curated missions (no crosswalking needed)
#' @param duplicate_log Tibble with catalog-wide duplicate information
#' @return String with updated image IDs. IDs that no longer exist are marked as "(image deleted)".
update_anomaly_image_references = function(text, crosswalk, duplicate_log) {
  if (is.na(text) || text == "") {
    return(text)
  }

  # Find all image ID patterns in the text (format: NNNNNN-NN_NNNNNN)
  pattern = "\\d{6}-\\d{2}_\\d{6}"
  matches = str_extract_all(text, pattern)[[1]]

  if (length(matches) == 0) {
    return(text)
  }

  has_crosswalk = !is.null(crosswalk) && nrow(crosswalk) > 0
  result = text

  for (old_id in matches) {
    new_id = old_id  # Default: keep original

    # Step 1: Apply crosswalk if available (formerly curated missions)
    if (has_crosswalk) {
      crosswalk_result = crosswalk |>
        filter(former_image_id == old_id) |>
        pull(current_image_id)

      if (length(crosswalk_result) > 0) {
        new_id = crosswalk_result[1]
      }
    }

    # Step 2: Handle crosswalk results and check for duplicates
    if (new_id == "none") {
      # Crosswalk says image doesn't exist in current data.
      # No duplicate lookup needed: if duplicates existed at this location,
      # one would remain after deduplication and the crosswalk would have found it.
      result = str_replace(result, fixed(old_id), "(image deleted)")
    } else if (new_id == "multiple") {
      # Multiple matches - ambiguous, keep original with warning
      warning(paste("Multiple matches found for image", old_id, "- keeping original"))
    } else {
      # new_id is a valid ID (either crosswalked or original)
      # Check if this ID was a removed duplicate and resolve to remaining
      remaining = get_remaining_duplicate(new_id, duplicate_log)

      if (!is.na(remaining) && remaining != new_id) {
        # The ID we have was a removed duplicate; use the remaining one
        result = str_replace(result, fixed(old_id), remaining)
      } else if (new_id != old_id) {
        # Crosswalked to a different ID
        result = str_replace(result, fixed(old_id), new_id)
      }
      # If new_id == old_id and not a removed duplicate, no change needed
    }
  }

  return(result)
}
```

---

### Step 7: Update `01_apply-curation-filters.R` - Load Duplicate Log

**File:** [01_apply-curation-filters.R](deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/01_apply-curation-filters.R)

**Add after line 27** (after loading curation notes):

```r
# ============================================================================
# Load catalog-wide duplicate log
# ============================================================================

cat("Loading catalog-wide duplicate log...\n")

if (!file.exists(CATALOG_DUPLICATE_IMAGES_LOG_PATH)) {
  stop("Catalog duplicate log not found. Run 05d_merge-duplicate-logs.R first: ",
       CATALOG_DUPLICATE_IMAGES_LOG_PATH)
}

duplicate_log = read_csv(CATALOG_DUPLICATE_IMAGES_LOG_PATH, col_types = cols(.default = "c"))
cat(sprintf("  Loaded %d duplicate image records\n", nrow(duplicate_log)))
```

---

### Step 8: Update `01_apply-curation-filters.R` - Modify Curation Note Update Logic

**File:** [01_apply-curation-filters.R](deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/01_apply-curation-filters.R)

**Important distinction:** Two types of missions require different handling:
- **Formerly curated missions**: Need crosswalking (old ID → new ID) AND duplicate resolution
- **Currently curated missions**: Only need duplicate resolution (IDs are already current, but may reference removed duplicates)

**Replace `update_curation_row()` function (lines 115-140)** with:

```r
update_curation_row = function(row, crosswalks, formerly_curated_missions, duplicate_log) {
  mission_id = row$mission_id
  is_formerly_curated = mission_id %in% formerly_curated_missions

  # Get crosswalk if this is a formerly curated mission
  crosswalk = NULL
  if (is_formerly_curated) {
    crosswalk = crosswalks[[mission_id]]
    if (is.null(crosswalk) || nrow(crosswalk) == 0) {
      warning(paste("No crosswalk available for formerly curated mission", mission_id))
      # Continue anyway - duplicate resolution can still be applied
    }
  }

  # Anomaly columns to update (extraneous_images handled separately in Step 9)
  anomaly_cols = c("collection_time_anomaly", "altitude_anomaly",
                   "spatial_anomaly", "camera_pitch_anomaly", "excess_images_anomaly",
                   "missing_images_anomaly", "other_anomalies", "anomaly_notes")

  for (col in anomaly_cols) {
    if (col %in% names(row) && !is.na(row[[col]])) {
      row[[col]] = update_anomaly_image_references(row[[col]], crosswalk, duplicate_log)
    }
  }

  return(row)
}
```

**Update the function call (around line 146)** to pass `duplicate_log`:

```r
updated = update_curation_row(row_list, crosswalks, formerly_curated_missions, duplicate_log)
```

---

### Step 9: Update `01_apply-curation-filters.R` - Handle Extraneous Image Filtering with Crosswalk and Duplicates

**File:** [01_apply-curation-filters.R](deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/01_apply-curation-filters.R)

**Purpose:** Parse extraneous images, apply crosswalk (for formerly curated missions), and resolve duplicates. This step handles all extraneous_images processing since Step 8 only handles anomaly columns.

**Replace the exclusion list parsing section (lines 156-182)** with enhanced logic:

```r
# ============================================================================
# Parse extraneous images and create exclusion lists per mission
# (with crosswalk application and duplicate resolution)
# ============================================================================

cat("Parsing extraneous images, applying crosswalks, and resolving duplicates...\n")

exclusion_lists = list()

# Helper to crosswalk and resolve an image ID
crosswalk_and_resolve = function(image_id, crosswalk, duplicate_log, current_image_ids) {
  working_id = image_id

  # Step 1: Apply crosswalk if available (formerly curated missions)
  if (!is.null(crosswalk) && nrow(crosswalk) > 0) {
    crosswalk_result = crosswalk |>
      filter(former_image_id == image_id) |>
      pull(current_image_id)

    if (length(crosswalk_result) > 0 && crosswalk_result[1] != "none" && crosswalk_result[1] != "multiple") {
      working_id = crosswalk_result[1]
    } else if (length(crosswalk_result) > 0 && crosswalk_result[1] == "none") {
      # Crosswalk says image doesn't exist in current data.
      # No duplicate lookup needed: if duplicates existed at this location,
      # one would remain after deduplication and the crosswalk would have found it.
      return(NA_character_)
    }
    # If "multiple" or not found in crosswalk, continue with original ID
  }

  # Step 2: If the working ID exists in current data, use it directly
  if (working_id %in% current_image_ids) {
    return(working_id)
  }

  # Step 3: Otherwise, try to find its remaining duplicate
  remaining = get_remaining_duplicate(working_id, duplicate_log)

  if (!is.na(remaining) && remaining %in% current_image_ids) {
    return(remaining)
  }

  # Neither exists - return NA (will be filtered out)
  return(NA_character_)
}

for (i in seq_len(nrow(curation_notes))) {
  row = curation_notes[i, ]
  mission_id = row$mission_id

  if (!is.na(row$extraneous_images) && row$extraneous_images != "") {
    tryCatch({
      # Parse extraneous images (handles range expansion)
      excluded_ids = parse_extraneous_images(row$extraneous_images, mission_id)

      if (length(excluded_ids) > 0) {
        # Load current image IDs for this mission to check existence
        current_filepath = file.path(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
                                     paste0(mission_id, "_image-metadata.gpkg"))
        if (file.exists(current_filepath)) {
          current_images = st_read(current_filepath, quiet = TRUE)
          current_image_ids = current_images$image_id

          # Get crosswalk if this is a formerly curated mission
          crosswalk = NULL
          is_formerly_curated = mission_id %in% formerly_curated_missions
          if (is_formerly_curated) {
            crosswalk = crosswalks[[mission_id]]
          }

          # Apply crosswalk and resolve duplicates for each ID
          resolved_ids = map_chr(excluded_ids, ~crosswalk_and_resolve(.x, crosswalk, duplicate_log, current_image_ids))

          # Remove NAs (images that don't exist)
          n_missing = sum(is.na(resolved_ids))
          resolved_ids = resolved_ids[!is.na(resolved_ids)]
          resolved_ids = unique(resolved_ids)

          if (n_missing > 0) {
            cat(sprintf("  Mission %s: %d image(s) not found in current data (skipped)\n",
                        mission_id, n_missing))
          }

          # Count how many required resolution (crosswalk or duplicate)
          n_resolved = sum(excluded_ids != resolved_ids[seq_along(excluded_ids)], na.rm = TRUE)
          if (n_resolved > 0) {
            cat(sprintf("  Mission %s: resolved %d image(s) via crosswalk/duplicates\n",
                        mission_id, n_resolved))
          }

          if (length(resolved_ids) > 0) {
            if (!(mission_id %in% names(exclusion_lists))) {
              exclusion_lists[[mission_id]] = character(0)
            }
            exclusion_lists[[mission_id]] = unique(c(exclusion_lists[[mission_id]], resolved_ids))

            cat(sprintf("  Mission %s: %d images to exclude\n", mission_id, length(resolved_ids)))
          }
        }
      }
    }, error = function(e) {
      stop(paste("Error parsing extraneous_images for mission", mission_id, ":", e$message))
    })
  }
}

cat(sprintf("  Total: %d missions with images to exclude\n", length(exclusion_lists)))
```

---

### Step 10: Handle Crosswalk Targets That Were Duplicates

**File:** [01_apply-curation-filters.R](deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/01_apply-curation-filters.R)

This is already handled by the new `update_image_references_with_duplicates()` function from Step 6, which:
1. Looks up each former image ID in the crosswalk
2. If the result is "none", calls `get_remaining_duplicate()` to find the kept duplicate
3. Uses the remaining duplicate's ID if found

**No additional changes needed** for this step - it's covered by the new function.

---
