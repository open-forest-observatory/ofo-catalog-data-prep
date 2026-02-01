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
