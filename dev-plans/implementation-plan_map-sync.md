# Implementation Plan: Sync Pan/Zoom Between Pre and Post-Curation Maps

## Goal
Make the two side-by-side maps on curation pages sync their pan and zoom, so when you pan/zoom one map, the other follows.

## Overview
Since the maps are in separate iframes, we'll use the PostMessage API to communicate between them. Each map will broadcast its view changes to the parent page, which relays them to the other map.

---

## Step 1: Add relay script to the template

**File to edit:** `deploy/drone-imagery-ingestion/10_drone-mission-web-catalog/templates/drone-mission-details_curation.md`

**What to do:**
Add the following script inside the existing `<script>` tag (before the closing `</script>`), after the `setView` function:

```javascript
// === Map Sync Functionality ===
// Relay pan/zoom messages between the two map iframes
window.addEventListener('message', function(e) {
  // Only handle mapSync messages
  if (!e.data || e.data.type !== 'mapSync') return;

  // Find both iframes
  const preFrame = document.querySelector('.pre-curation-panel iframe');
  const postFrame = document.querySelector('.post-curation-panel iframe');

  if (!preFrame || !postFrame) return;

  // Relay message to the other iframe
  if (e.source === preFrame.contentWindow) {
    postFrame.contentWindow.postMessage(e.data, '*');
  } else if (e.source === postFrame.contentWindow) {
    preFrame.contentWindow.postMessage(e.data, '*');
  }
});
```

---

## Step 2: Create a JavaScript snippet file for map sync

**File to create:** `deploy/drone-imagery-ingestion/10_drone-mission-web-catalog/templates/map-sync-snippet.js`

**Content:**
```javascript
<script>
(function() {
  // Wait for the map to be ready
  // Leaflet htmlwidgets store the map in HTMLWidgets.widgets
  function getLeafletMap() {
    if (typeof HTMLWidgets !== 'undefined' && HTMLWidgets.widgets && HTMLWidgets.widgets.length > 0) {
      for (var i = 0; i < HTMLWidgets.widgets.length; i++) {
        var widget = HTMLWidgets.widgets[i];
        if (widget.getMap) {
          return widget.getMap();
        }
      }
    }
    // Alternative: check for Leaflet map on the widget element
    var mapEl = document.querySelector('.leaflet-container');
    if (mapEl && mapEl._leaflet_map) {
      return mapEl._leaflet_map;
    }
    return null;
  }

  function initSync() {
    var map = getLeafletMap();
    if (!map) {
      // Retry if map not ready yet
      setTimeout(initSync, 100);
      return;
    }

    var syncInProgress = false;

    // Broadcast view changes to parent
    map.on('moveend', function() {
      if (syncInProgress) return;
      parent.postMessage({
        type: 'mapSync',
        center: [map.getCenter().lat, map.getCenter().lng],
        zoom: map.getZoom()
      }, '*');
    });

    // Receive sync from other map via parent
    window.addEventListener('message', function(e) {
      if (!e.data || e.data.type !== 'mapSync') return;
      syncInProgress = true;
      map.setView(e.data.center, e.data.zoom, {animate: false});
      // Reset flag after a short delay to allow the moveend event to fire
      setTimeout(function() { syncInProgress = false; }, 50);
    });
  }

  // Start initialization when DOM is ready
  if (document.readyState === 'complete') {
    initSync();
  } else {
    window.addEventListener('load', initSync);
  }
})();
</script>
```

---

## Step 3: Modify the R function to inject the sync script into map HTML files

**File to edit:** `src/web-catalog-creation_drone-imagery-catalog.R`

**What to do:**
Find the `make_mission_details_map` function and modify it to append the sync JavaScript to the saved HTML file.

1. First, at the top of the file (after the other source/library statements), add code to read the snippet:

```r
# Load map sync JavaScript snippet for curation pages
MAP_SYNC_SNIPPET_PATH <- "deploy/drone-imagery-ingestion/10_drone-mission-web-catalog/templates/map-sync-snippet.js"
MAP_SYNC_SNIPPET <- if (file.exists(MAP_SYNC_SNIPPET_PATH)) {
  readLines(MAP_SYNC_SNIPPET_PATH, warn = FALSE) |> paste(collapse = "\n")
} else {
  ""
}
```

2. Find where `saveWidget` is called in `make_mission_details_map`. After the widget is saved, add code to append the sync script:

```r
# Append map sync script for iframe communication
if (nzchar(MAP_SYNC_SNIPPET)) {
  html_content <- readLines(output_path, warn = FALSE)
  # Insert before closing </body> tag
  html_content <- sub("</body>", paste0(MAP_SYNC_SNIPPET, "\n</body>"),
                       paste(html_content, collapse = "\n"))
  writeLines(html_content, output_path)
}
```

