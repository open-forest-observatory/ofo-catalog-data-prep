# Implementation Plan: Post-Curation Metadata Viewer for Curation Web Pages

## Overview

Add post-curation metadata visualization to the drone data curation webpages, displaying pre-curation and post-curation image point maps and mission metadata tables side by side with a JavaScript toggle to switch between views.

## User Requirements Summary

- Display two image point maps and two mission metadata tables side by side (pre-curation on left, post-curation on right)
- JavaScript toggle to switch between: Pre-curation only, Post-curation only, Both (default)
- For missions without post-curation data: show both panels but display "No data available" in post-curation panel
- Read post-curation metadata from staged location following existing pattern
- Reuse existing functions where possible

---

## Files to Modify

| File | Change Type |
|------|-------------|
| `deploy/drone-imagery-ingestion/00_set-constants.R` | Add new constants |
| `deploy/drone-imagery-ingestion/10_drone-mission-web-catalog/11_create-drone-data-curation-webpages.R` | Major modifications |
| `deploy/drone-imagery-ingestion/10_drone-mission-web-catalog/templates/drone-mission-details_curation.md` | Major modifications |
| `src/web-catalog-creation_drone-imagery-catalog.R` | Minor refactoring (if needed) |

---

## Step-by-Step Implementation

### Step 1: Add New Constants for Post-Curation Metadata Paths

**File:** [00_set-constants.R](deploy/drone-imagery-ingestion/00_set-constants.R)

**Location:** Add after line 246 (after `S3_LISTING_FILEPATH`)

**What to add:**
- `POST_CURATION_MISSION_METADATA_FILEPATH` - path to staged post-curation mission metadata
- `POST_CURATION_IMAGE_METADATA_FILEPATH` - path to staged post-curation image metadata

**Pattern to follow:**
```
BASE_DATA_PATH/05_drone-imagery-web-catalog/01_metadata/post-curation/mission-metadata.gpkg
BASE_DATA_PATH/05_drone-imagery-web-catalog/01_metadata/post-curation/image-metadata.gpkg
```

Also add new constants for post-curation output directories to keep widgets separate:
- `CURATION_POST_MISSION_DETAILS_DATATABLE_DIR` - for post-curation datatables
- `CURATION_POST_MISSION_DETAILS_MAP_DIR` - for post-curation maps

---

### Step 2: Update the Curation Webpage Generation Script

**File:** [11_create-drone-data-curation-webpages.R](deploy/drone-imagery-ingestion/10_drone-mission-web-catalog/11_create-drone-data-curation-webpages.R)

#### 2a. Load Post-Curation Metadata

**Location:** After line 34 (after loading primary image points)

**Logic:**
1. Attempt to load post-curation mission metadata from `POST_CURATION_MISSION_METADATA_FILEPATH`
2. Attempt to load post-curation image metadata from `POST_CURATION_IMAGE_METADATA_FILEPATH`
3. If files don't exist, set to `NULL` and log a warning (script should still work)
4. Create a vector of `missions_with_post_curation_data` by checking which mission_ids exist in post-curation metadata

**Guidance:**
- Use `file.exists()` checks before loading
- Store loaded data in variables like `post_curation_mission_metadata` and `post_curation_image_points`

#### 2b. Generate Post-Curation Widgets

**Location:** Inside `make_mission_curation_page` function (around line 180)

**Logic:**
1. Check if current mission has post-curation data (mission_id exists in post-curation metadata)
2. If yes:
   - Call `make_mission_details_map()` for post-curation image points, saving to `CURATION_POST_MISSION_DETAILS_MAP_DIR`
   - Call `make_mission_details_datatable()` for post-curation mission metadata, saving to `CURATION_POST_MISSION_DETAILS_DATATABLE_DIR`
   - Store the returned paths
3. If no:
   - Set paths to `NA` or a special marker value

**Key consideration:** The existing functions `make_mission_details_map()` and `make_mission_details_datatable()` should work as-is since post-curation data has the same format. You're just calling them with different data and saving to different directories.

#### 2c. Update Template Rendering Call

**Location:** Inside `make_mission_curation_page` function, at the `render_mission_details_page()` call (line 233)

**Changes needed:**
- Pass additional template variables for post-curation:
  - `post_curation_map_html_path` - path to post-curation map (or NA)
  - `post_curation_datatable_html_path` - path to post-curation datatable (or NA)
  - `has_post_curation_data` - boolean flag

**Note:** You may need to modify `render_mission_details_page()` function signature to accept these new parameters, OR create a separate rendering function for curation pages.

---

### Step 3: Update the Jinjar Template

**File:** [drone-mission-details_curation.md](deploy/drone-imagery-ingestion/10_drone-mission-web-catalog/templates/drone-mission-details_curation.md)

#### 3a. Add View Toggle UI

**Location:** After line 24 (after the "Internal curation view" paragraph)

**Implementation:**
- Add a `<div>` containing three radio buttons or toggle buttons:
  - "Pre-curation"
  - "Post-curation"
  - "Both" (checked by default)
- Style buttons inline and centered
- Each button should have an `onclick` handler that calls a JavaScript function

#### 3b. Restructure Content Layout for Side-by-Side Display

**Current structure (lines 26-35):**
```
## Image locations and attributes
<iframe src="...map...">

## Mission metadata
<iframe src="...datatable...">
```

**New structure:**
```
## Image locations and attributes
<div class="row">
  <div class="col pre-curation-panel">
    <h4>Pre-curation</h4>
    <iframe src="{* map_html_path *}">
  </div>
  <div class="col post-curation-panel">
    <h4>Post-curation</h4>
    {% if has_post_curation_data %}
    <iframe src="{* post_curation_map_html_path *}">
    {% else %}
    <p class="no-data-message">No post-curation data available</p>
    {% endif %}
  </div>
</div>

## Mission metadata
<div class="row">
  <div class="col pre-curation-panel">
    <h4>Pre-curation</h4>
    <iframe src="{* datatable_html_path *}">
  </div>
  <div class="col post-curation-panel">
    <h4>Post-curation</h4>
    {% if has_post_curation_data %}
    <iframe src="{* post_curation_datatable_html_path *}">
    {% else %}
    <p class="no-data-message">No post-curation data available</p>
    {% endif %}
  </div>
</div>
```

**Note:** Use Bootstrap's grid classes (`row`, `col-sm-6`) for responsive side-by-side layout.

#### 3c. Add JavaScript Toggle Logic

**Location:** At the bottom of the template (after line 45)

**Implementation:**
```javascript
<script>
function setView(view) {
  const prePanels = document.querySelectorAll('.pre-curation-panel');
  const postPanels = document.querySelectorAll('.post-curation-panel');

  if (view === 'pre') {
    prePanels.forEach(p => { p.style.display = 'block'; p.className = 'col-12 pre-curation-panel'; });
    postPanels.forEach(p => p.style.display = 'none');
  } else if (view === 'post') {
    prePanels.forEach(p => p.style.display = 'none');
    postPanels.forEach(p => { p.style.display = 'block'; p.className = 'col-12 post-curation-panel'; });
  } else { // both
    prePanels.forEach(p => { p.style.display = 'block'; p.className = 'col-sm-6 pre-curation-panel'; });
    postPanels.forEach(p => { p.style.display = 'block'; p.className = 'col-sm-6 post-curation-panel'; });
  }
}

// Initialize to 'both' view on page load
document.addEventListener('DOMContentLoaded', function() {
  setView('both');
});
</script>
```

**Styling considerations:**
- Add CSS for the "no data available" message (gray italic text, centered)
- Ensure iframes resize appropriately in both full-width and half-width modes

---

### Step 4: Update render_mission_details_page Function (if needed)

**File:** [web-catalog-creation_drone-imagery-catalog.R](src/web-catalog-creation_drone-imagery-catalog.R)

**Location:** `render_mission_details_page()` function starting at line 596

**Decision point:** The current function is used for both public catalog pages and curation pages. You have two options:

**Option A (Recommended):** Add optional parameters with default values
- Add `post_curation_map_html_path = NA`
- Add `post_curation_datatable_html_path = NA`
- Add `has_post_curation_data = FALSE`
- Pass these to `jinjar::render()` call

**Option B:** Create a separate `render_curation_details_page()` function
- Copy and modify for curation-specific needs
- More code duplication but cleaner separation
