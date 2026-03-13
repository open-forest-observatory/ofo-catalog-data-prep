# Turn a data store path into a public HTTP url
cyverse_url = function(data_store_path) {
  # data_store_path is the path to a file in the CyVerse Data Store, e.g. "/iplant/home/shared/ofo/public/missions/000001/processed_000001-0001/full/chm-mesh.tif"
  # Returns a public HTTP URL to the file, e.g. "https://data.cyverse.org/dav-anon/iplant/home/shared/ofo/public/missions/000001/processed_000001-0001/full/chm-mesh.tif"

  base_url = "https://data.cyverse.org/dav-anon/"
  url = paste0(base_url, data_store_path)
  return(url)
}

# Query cyverse data store for a list of files matching a provided directory pattern and file
# pattern (with % as wildcard)
cyverse_list_files = function(dir_pattern, file_pattern) {
  # dir_pattern is the directory pattern to search for, e.g. "/iplant/projects/ofo/public/missions/%/processed_%/full"
  # file_pattern is the file pattern to search for, e.g. "chm-mesh.tif"
  call = paste0("iquest --no-page '%s/%s' \"select COLL_NAME, DATA_NAME where COLL_NAME like '", dir_pattern, "' and DATA_NAME like '", file_pattern, "'\"")
  output = system(call, intern = TRUE)
  return(output)
}

cyverse_list_dirs_recursive = function(dir_pattern) {
  call = paste0("iquest --no-page '%s' \"select COLL_NAME where COLL_NAME like '", dir_pattern, "'\"")
  output = system(call, intern = TRUE)
  return(output)
}

cyverse_list_subdirs = function(dir) {
  call = paste0("ils ", dir, "")
  output = system(call, intern = TRUE)

  # Remove the first line which is the header
  output = output[-1]

  # Remove the leading and trailing whitespace
  output = str_trim(output)

  # Remove the leading "C- "
  output = str_replace(output, "^C- ", "")

  return(output)
}

df_from_url = function(url) {
  temp_file = tempfile(fileext = ".csv")
  download.file(url, temp_file, quiet = TRUE, method = "wget")
  read_csv(temp_file, col_types = cols(.default = "c"))
}

sf_from_url = function(url) {
  temp_file = tempfile(fileext = ".gpkg")
  download.file(url, temp_file, quiet = TRUE, method = "wget")
  sf = st_read(temp_file, quiet = TRUE)
  return(sf)
}

# Pull in all attributes to the plot summary table, compute relevant columns for display, (like
# plot area in ha instead of m2), round numeric columns, create links to project and dataset pages
compile_mission_summary_data = function(mission_level_metadata, base_ofo_url, mission_details_dir, dataset_type = "mission") {
  # dataset_type can be "mission" or "sub-mission"

  # TODO: The dataset_id determination is flexible to either mission or sub-mission, but there is
  # other hardcoded references to mission that would need to be updated if sub-mission is a
  # possibility

  if (dataset_type == "mission") {
    dataset_id_column = "mission_id"
  } else if (dataset_type == "sub-mission") {
    dataset_id_column = "sub_mission_id"
  } else {
    stop("dataset_type must be either 'mission' or 'sub-mission'")
  }

  mission_level_metadata$dataset_id = pull(mission_level_metadata, dataset_id_column)

  # Pre-process the display text for the relevant attributes, based on the metadata in the database
  d = mission_level_metadata |>
    dplyr::mutate(
      # delete if still commented: overlap_combined_nominal = paste(overlap_front_nominal, overlap_side_nominal, sep = "/"),
      dataset_id_link = paste0('<a href="', base_ofo_url, mission_details_dir, dataset_id, '/"', ' target="_PARENT">', dataset_id, "</a>"),
      time_range_local_derived = paste0(earliest_time_local_derived, " to ", latest_time_local_derived),
      overlap_front_side_nominal = paste0(overlap_front_nominal, "/", overlap_side_nominal),
      image_dimensions_derived = paste0(resolution_x_derived, " x ", resolution_y_derived),
      altitude_agl_mean_derived = NA,
      embargoed = FALSE, # Dummy, need to remove if add this attribute to the database
      display_message = NA) # Dummy, need to remove if add this attribute to the database

  return(d)

}



make_mission_catalog_datatable = function(mission_summary,
                                       website_static_path,
                                       datatable_header_files_dir,
                                       mission_catalog_datatable_dir,
                                       mission_catalog_datatable_filename) {

  d = mission_summary |>
    sf::st_drop_geometry() |>
    select("ID" = dataset_id_link,
           "Area (ha)" = area_derived,
           "Date" = earliest_date_derived,
           "Altitude (m) (N)" = altitude_agl_nominal,
           "Overlap (N)" = overlap_front_side_nominal,
           "Camera pitch" = camera_pitch_derived,
           "Terrain follow (N)" = terrain_follow,
           "Flight pattern" = flight_pattern,
           "Image count" = image_count_derived,
           "RTK (N)" = rtk_nominal,
           "RTK images (%)" = percent_images_rtk_derived,
           "Contributor dataset name" = contributor_dataset_name,
           "Project" = project_id,
           "Aircraft" = aircraft_model_name,
           dataset_id) |>
    arrange(dataset_id) |>
    select(-dataset_id)

  # Prep formatting code to pass to datatable creation
  format_js = DT::JS("function(settings, json) {",
                     "$('body').css({'font-family': 'Arial'});",
                     "}")

  dt = DT::datatable(d,
                     rownames = FALSE,
                     escape = FALSE,
                     options = list(paging = FALSE, initComplete = format_js))

  # Save the datatable HTML to the website repo
  save_widget_html(dt,
                   website_static_path = website_static_path,
                   header_files_dir = datatable_header_files_dir,
                   html_dir = mission_catalog_datatable_dir,
                   html_filename = mission_catalog_datatable_filename,
                   delete_folder_first = TRUE)

  return(dt)

}



make_mission_catalog_map = function(mission_summary,
                                   website_static_path,
                                   leaflet_header_files_dir,
                                   mission_catalog_map_dir,
                                   mission_catalog_map_filename) {

  mission_summary = mission_summary |>
    dplyr::mutate(popup = paste0(
      "<b>Mission ID: </b>", dataset_id_link, "<br>",
      "<b>Date: </b>", earliest_date_derived, "<br>",
      "<b>Altitude (m): </b>", altitude_agl_nominal, "<br>",
      "<b>Camera pitch (deg): </b>", camera_pitch_derived, "<br>",
      "<b>Overlap (front/side): </b>", overlap_front_side_nominal, "<br>"
    ))

  mission_centroids = sf::st_centroid(mission_summary)

  m = leaflet() |>
    addTiles(options = providerTileOptions(maxZoom = 16)) |>
    addMarkers(data = mission_centroids, popup = ~popup, clusterOptions = markerClusterOptions(freezeAtZoom = 16)) |>
    addPolygons(data = mission_summary, popup = ~popup, group = "bounds") |>
    addProviderTiles(providers$Esri.WorldTopoMap, group = "Topo", options = providerTileOptions(maxZoom = 16)) |>
    addProviderTiles(providers$Esri.WorldImagery, group = "Imagery") |>
    groupOptions("bounds", zoomLevels = 13:20) |>
    addLayersControl(baseGroups = c("Topo", "Imagery"),
                     options = layersControlOptions(collapsed = FALSE))

  save_widget_html(m,
                   website_static_path = website_static_path,
                   header_files_dir = leaflet_header_files_dir,
                   html_dir = mission_catalog_map_dir,
                   html_filename = mission_catalog_map_filename,
                   delete_folder_first = TRUE)

  return(m)

}



sub_mission_from_image_id = function(image_id) {
  mission_w_sub_mission = str_split(image_id, fixed("_")) |> purrr::map(1) |> unlist()
  sub_mission = str_split(mission_w_sub_mission, fixed("-")) |> purrr::map(2) |> unlist()
  return(sub_mission)
}

#### Make mission-level leaflet map of image points and flight path
make_mission_details_map = function(mission_summary_foc,
                                    mission_points_foc,
                                    mission_polygons_for_mission_details_map,
                                    mission_centroids,
                                    website_static_path,
                                    leaflet_header_files_dir,
                                    mission_details_map_dir) {

  dataset_id = mission_summary_foc$dataset_id

  # Recast some attributes. TODO: This may have become necessary because when compiling all
  # mission metadata, we set all columns to character so they could be combined. Consider fixing
  # that further upstream by using readr's infer types function
  mission_points_foc = mission_points_foc |>
    mutate(
      altitude_asl_drone = as.numeric(altitude_asl_drone),
      camera_pitch = as.numeric(camera_pitch),
      datetime_local = as.POSIXct(datetime_local, tz = "UTC")
    )

  mission_points_foc = mission_points_foc |>
    # Get which sub-mission each image is from
    mutate(sub_mission = sub_mission_from_image_id(image_id)) |>
    # Create a column "hours elapsed since mission start" for legend coloring
    mutate(time_secs = as.numeric(datetime_local))

  initial_time = min(mission_points_foc$time_secs)

  mission_points_foc = mission_points_foc |>
    mutate(hours_elapsed = (time_secs - initial_time)/60/60) |>
    # Create a popup text
    mutate(popup = paste0("<b>Image ID: </b>", image_id, "<br>"))

  # Optoinal addl rows for popup
    #   "<b>Altitude ASL: </b>", altitude_asl_drone, " m<br>",
    # "<b>Camera pitch: </b>", camera_pitch, " deg<br>",
    # "<b>RTK fix: </b>", rtk_fix, "<br>",
    # "<b>Capture datetime: </b>", datetime_local, "<br>",
    # "<b>Hours elapsed: </b>", round(hours_elapsed, 2), "<br>",
    # "<b>Sub-mission: </b>", sub_mission

  # Compute flightpath
  flightpath = mission_points_foc |>
    summarize(do_union = FALSE) |>
    st_cast("LINESTRING")

  # Make leaflet map

  js_for_legend = function(x) {
    htmlwidgets::onRender(x, "
      function(el, x) {
        var updateLegend = function () {
          var selectedGroup = document.querySelectorAll('input:checked')[0].nextSibling.innerText.substr(1).replace(/[^a-zA-Z]+/g, '');
          document.querySelectorAll('.legend').forEach( a => a.hidden=true );
          document.querySelectorAll('.legend').forEach( l => { if (l.classList.contains(selectedGroup)) l.hidden=false; } );
        };
        updateLegend();
        this.on('baselayerchange', el => updateLegend());
        }"
    )
  }

  # Define color palettes
  pal_hourselapsed = colorNumeric("viridis", domain = mission_points_foc$hours_elapsed)
  pal_pitch = colorNumeric("viridis", domain = mission_points_foc$camera_pitch)
  pal_rtk = colorFactor("viridis", domain = mission_points_foc$rtk)
  pal_sub_mission = colorFactor("viridis", domain = mission_points_foc$sub_mission)
  # All the values of the altitude may be the same. In this case, use a categorical color map to
  # ensure the value is still shown on the legend.
  if (min(mission_points_foc$altitude_asl_drone) == max(mission_points_foc$altitude_asl_drone)) {
    pal_asl = colorFactor("viridis", domain = mission_points_foc$altitude_asl_drone)
  } else {
    pal_asl = colorNumeric("viridis", domain = mission_points_foc$altitude_asl_drone)
  }

  m = leaflet() |>
    addPolygons(data = mission_summary_foc, group = "bounds",
                fillOpacity = 0) |>
    addProviderTiles(providers$Esri.WorldTopo, group = "Topo",
                     options = providerTileOptions(minZoom = 1, maxZoom = 20)) |>
    addProviderTiles(providers$Esri.WorldImagery, group = "Imagery",
                     options = providerTileOptions(minZoom = 1, maxZoom = 20)) |>
    addLayersControl(overlayGroups = c("Imagery", "Topo"), #  ", Nearby missions"
                     baseGroups = c("Hours elapsed", "Altitude ASL (m)", "Camera pitch (deg)", "RTK fix", "Sub-mission"),
                     options = layersControlOptions(collapsed = FALSE)) |>
    # # Nearby mission polygons and centroids
    # addMarkers(data = mission_centroids, popup = ~dataset_id_link, clusterOptions = markerClusterOptions(freezeAtZoom = 16), group = "Nearby missions") |>
    # addPolygons(data = mission_polygons_for_mission_details_map, popup = ~dataset_id_link, group = "Nearby missions") |>
    # Flight lines
    addPolylines(data = flightpath, color = "black", weight = 1, group = "Flight path") |>
    # Altitude ASL
    addCircleMarkers(data = mission_points_foc,
                    radius = 3,
                      stroke = FALSE,
                      fillOpacity = 1,
                    #  popup = ~popup,
                      color = pal_asl(mission_points_foc$altitude_asl_drone),
                      group = "Altitude ASL (m)") |>
    addLegend(pal = pal_asl,
              values = mission_points_foc$altitude_asl_drone,
              title = "Altitude ASL (m)", opacity = 1,
              group = "AltitudeASLm",
              className = "info legend AltitudeASLm") |>
    # Camera pitch
    addCircleMarkers(data = mission_points_foc,
                    radius = 3,
                      stroke = FALSE,
                      fillOpacity = 1,
                    #  popup = ~popup,
                      color = pal_pitch(mission_points_foc$camera_pitch),
                      group = "Camera pitch (deg)") |>
    addLegend(pal = pal_pitch,
              values = mission_points_foc$camera_pitch,
              title = "Camera pitch (deg)", opacity = 1,
              group = "Camerapitchdeg",
              className = "info legend Camerapitchdeg") |>
    # RTK
    addCircleMarkers(data = mission_points_foc,
                    radius = 3,
                      stroke = FALSE,
                      fillOpacity = 1,
                    #  popup = ~popup,
                      color = pal_rtk(mission_points_foc$rtk),
                      group = "RTK fix") |>
    addLegend(pal = pal_rtk,
              values = mission_points_foc$rtk,
              title = "RTK fix", opacity = 1,
              group = "RTKfix",
              className = "info legend RTKfix") |>
    # Hours elapsed
    addCircleMarkers(data = mission_points_foc,
                    radius = 3,
                      stroke = FALSE,
                      fillOpacity = 1,
                    #  popup = ~popup,
                      color = pal_hourselapsed(mission_points_foc$hours_elapsed),
                      group = "Hours elapsed") |>
    addLegend(pal = pal_hourselapsed,
              values = mission_points_foc$hours_elapsed,
              title = "Hours elapsed", opacity = 1,
              group = "Hourselapsed",
              className = "info legend Hourselapsed") |>
    # Sub-mission
    addCircleMarkers(data = mission_points_foc,
                    radius = 3,
                      stroke = FALSE,
                      fillOpacity = 1,
                    #  popup = ~popup,
                      color = pal_sub_mission(mission_points_foc$sub_mission),
                      group = "Sub-mission") |>
    addLegend(pal = pal_sub_mission,
              values = mission_points_foc$sub_mission,
              title = "Sub-mission", opacity = 1,
              group = "Submission",
              className = "info legend Submission") |>
    # Invisible markers on top of all for popup
    addCircleMarkers(data = mission_points_foc,
                    radius = 10,
                    stroke = FALSE,
                    fillOpacity = 0,
                    popup = ~popup,
                    group = "dummyforpopup") |>
    hideGroup("Imagery") |>
    hideGroup("Topo") |>
    hideGroup("Nearby missions") |>
    js_for_legend()

  # Customize background color (can also use this to make transparent)
  backg <- htmltools::tags$style(".leaflet-container { background: rgba(200,200,200,1) }")
  m = prependContent(m, backg)

  # -- Save map HTML to website repo
  mission_details_map_filename = paste0(dataset_id, ".html")
  save_widget_html(m,
                    website_static_path = website_static_path,
                    header_files_dir = leaflet_header_files_dir,
                    html_dir = mission_details_map_dir,
                    html_filename = mission_details_map_filename)

  # Record where it was saved to
  map_html_path = paste(mission_details_map_dir, mission_details_map_filename, sep = "/")

  return(map_html_path)

}



#### Make mission-level leaflet map of image points and flight path
make_itd_map = function(mission_summary_foc,
                                    itd_points_foc,
                                    mission_polygons_for_mission_details_map,
                                    mission_centroids,
                                    website_static_path,
                                    leaflet_header_files_dir,
                                    itd_map_dir) {

  dataset_id = mission_summary_foc$dataset_id

  #mission_summary_foc is the polygon. Buffer it in by 10 because that's the area if detected trees
  #that we retained

  mission_summary_foc = mission_summary_foc |>
    transform_to_local_utm() |>
    st_buffer(-10) |>
    st_transform(4326)

  itd_points_foc = itd_points_foc |>
    mutate(height = round(Z, 1)) |>
    # Create a popup text
    mutate(popup = paste0("<b>Height: </b>", height, " m<br> 
                          <b>Predicted species: </b> <i>coming soon</i><br>")) |>
    st_transform(4326)

  # Optoinal addl rows for popup
    #   "<b>Altitude ASL: </b>", altitude_asl_drone, " m<br>",
    # "<b>Camera pitch: </b>", camera_pitch, " deg<br>",
    # "<b>RTK fix: </b>", rtk_fix, "<br>",
    # "<b>Capture datetime: </b>", datetime_local, "<br>",
    # "<b>Hours elapsed: </b>", round(hours_elapsed, 2), "<br>",
    # "<b>Sub-mission: </b>", sub_mission

  # Make leaflet map

  js_for_legend = function(x) {
    htmlwidgets::onRender(x, "
      function(el, x) {
        var updateLegend = function () {
          var selectedGroup = document.querySelectorAll('input:checked')[0].nextSibling.innerText.substr(1).replace(/[^a-zA-Z]+/g, '');
          document.querySelectorAll('.legend').forEach( a => a.hidden=true );
          document.querySelectorAll('.legend').forEach( l => { if (l.classList.contains(selectedGroup)) l.hidden=false; } );
        };
        updateLegend();
        this.on('baselayerchange', el => updateLegend());
        }"
    )
  }

  # Define color palettes
  pal_height = colorNumeric("viridis", domain = itd_points_foc$height)

  m = leaflet() |>
    addPolygons(data = mission_summary_foc, group = "bounds",
                fillOpacity = 0) |>
    addProviderTiles(providers$Esri.WorldTopo, group = "Topo",
                     options = providerTileOptions(minZoom = 1, maxZoom = 20)) |>
    addProviderTiles(providers$Esri.WorldImagery, group = "Imagery",
                     options = providerTileOptions(minZoom = 1, maxZoom = 20)) |>
    addLayersControl(overlayGroups = c("Imagery", "Topo"), #  ", Nearby missions"
                     options = layersControlOptions(collapsed = FALSE)) |>
    # # Nearby mission polygons and centroids
    # addMarkers(data = mission_centroids, popup = ~dataset_id_link, clusterOptions = markerClusterOptions(freezeAtZoom = 16), group = "Nearby missions") |>
    # addPolygons(data = mission_polygons_for_mission_details_map, popup = ~dataset_id_link, group = "Nearby missions") |>
    # ITD points
    addCircleMarkers(data = itd_points_foc,
                    radius = itd_points_foc$height/5,
                      stroke = FALSE,
                      fillOpacity = 1,
                    #  popup = ~popup,
                      color = pal_height(itd_points_foc$height),
                      group = "Height") |>
    addLegend(pal = pal_height,
              values = itd_points_foc$height,
              title = "Height", opacity = 1,
              group = "Height",
              className = "info legend Height") |>
    # Invisible markers on top of all for popup
    addCircleMarkers(data = itd_points_foc,
                    radius = 10,
                    stroke = FALSE,
                    fillOpacity = 0,
                    popup = ~popup,
                    group = "dummyforpopup") |>
    hideGroup("Imagery") |>
    hideGroup("Topo") |>
    hideGroup("Nearby missions") |>
    js_for_legend()

  # Customize background color (can also use this to make transparent)
  backg <- htmltools::tags$style(".leaflet-container { background: rgba(200,200,200,1) }")
  m = prependContent(m, backg)

  # -- Save map HTML to website repo
  itd_map_filename = paste0(dataset_id, ".html")
  save_widget_html(m,
                    website_static_path = website_static_path,
                    header_files_dir = leaflet_header_files_dir,
                    html_dir = itd_map_dir,
                    html_filename = itd_map_filename)

  # Record where it was saved to
  map_html_path = paste(itd_map_dir, itd_map_filename, sep = "/")

  return(map_html_path)

}



make_mission_details_datatable = function(mission_summary_foc,
                                          website_static_path,
                                          datatable_header_files_dir,
                                          mission_details_datatable_dir) {

  # Select the desired columns, name them, and pivot longer
  d = mission_summary_foc |>
    st_drop_geometry() |>
    # Select just what's needed for a datatable
    dplyr::select(
      "Mission ID" = dataset_id,
      "Sub-mission IDs" = sub_mission_ids,
      "Date" = earliest_date_derived,
      "Altitude (nominal) (m)" = altitude_agl_nominal,
      "Altitude AGL mean (actual) (m)" = altitude_agl_mean_derived,
      "Flight pattern" = flight_pattern,
      "Overlap (nominal) (front/side)" = overlap_front_side_nominal,
      "Camera pitch (deg up from nadir)" = camera_pitch_derived,
      "Smart oblique" = smart_oblique_derived,
      "Terrain follow (nominal)" = terrain_follow,
      "Terrain follow fidelity" = flight_terrain_correlation_derived,
      "Percent RTK" = percent_images_rtk_derived,
      "Flight speed (m/s)" = flight_speed_derived,
      "Exposure median (sec)" = exposure_median_derived,
      "Exposure CV" = exposure_cv_derived,
      "Exposure compensation (nominal)" = exposure_compensation,
      "White balance mode" = white_balance_mode_derived,
      "White balance percent" = white_balance_pct_mode_derived,
      "Flight time range (local)" = time_range_local_derived,
      "Aircraft model" = aircraft_model_name,
      "Sensor model" = sensor_name,
      "Flight planner" = flight_planner_name,
      "Base station latitude" = base_lat,
      "Base station longitude" = base_lon,
      "Base station altitude (m)" = base_alt,
      "Permanent base marker" = base_marked_permanently,
      "Flight footprint area (ha)" = area_derived,
      "Image count" = image_count_derived,
      "Image dimensions" = image_dimensions_derived,
      "Dataset size (GB)" = file_size_derived,
      "File format" = file_format_derived,
      "Project ID" = project_id,
      "Contributor dataset name" = contributor_dataset_name,
      "Creator" = contributor_names,
      "License" = license
    ) |>
      # Pivot longer
      dplyr::mutate(across(everything(), as.character)) |>
      tidyr::pivot_longer(cols = everything(), names_to = "Attribute", values_to = "Value")

  # Prep the JS formatting code
  format_js = DT::JS("function(settings, json) {",
                     "$('body').css({'font-family': 'Arial'});",
                     "}")

  dt = datatable(d, rownames = FALSE, escape = TRUE,
                 colnames = NULL,
                 options = list(paging = FALSE, scrollY = "100%",
                                dom = 't',
                                bSort = FALSE,
                                autoWidth = TRUE,
                                columnDefs = list(list(width = '40%', targets = "Attribute")),
                                #headerCallback = headerCallbackJS,
                                initComplete = format_js))

  dt$sizingPolicy$browser$padding = 0
  dt$sizingPolicy$browser$fill = FALSE

  # Save the datatable HTML to the website repo
  mission_details_datatable_filename = paste0(mission_summary_foc$dataset_id, ".html")
  save_widget_html(dt,
                   website_static_path = website_static_path,
                   header_files_dir = datatable_header_files_dir,
                   html_dir = mission_details_datatable_dir,
                   html_filename = mission_details_datatable_filename)

  # Record where it was saved to
  datatable_html_path = paste(mission_details_datatable_dir, mission_details_datatable_filename, sep = "/")

  return(datatable_html_path)

}


# Compute S3 product URLs for a given dataset (mission or composite)
# Returns a named list of product existence flags and URLs
compute_s3_product_urls = function(dataset_id, s3_file_listing_foc, data_server_base_url) {

  processed_folder = paste0("photogrammetry_", PHOTOGRAMMETRY_CONFIG_ID)
  processed_products = s3_file_listing_foc |> filter(str_detect(filepath, processed_folder))

  filepath_parts = str_split(processed_products$filepath, fixed("/"))
  part_3 = purrr::map_chr(filepath_parts, 3)
  itd_folders = part_3[str_which(part_3, "^(itd_|detected_trees_|detected-trees_)")] |> unique() |> sort(decreasing = TRUE)
  itd_folder_mostrecent = itd_folders[1]
  itd_path_mostrecent = file.path(dataset_id, processed_folder, itd_folder_mostrecent)
  ttops_file_path = file.path(itd_path_mostrecent, paste0(dataset_id, "_treetops.gpkg"))

  # Check if products exist and if so, get the URLs needed to add them to the page

  # Orthomosaic
  ortho_path_thumb = file.path(dataset_id, processed_folder, "thumbnails", paste0(dataset_id, "_ortho-dsm-ptcloud.png"))
  ortho_path_full = file.path(dataset_id, processed_folder, "full", paste0(dataset_id, "_ortho-dsm-ptcloud.tif"))
  ortho_exists = ortho_path_thumb %in% processed_products$filepath
  ortho_url_thumb = paste(data_server_base_url, ortho_path_thumb, sep = "/") |> strip_double_slashes()
  ortho_url_full = paste(data_server_base_url, ortho_path_full, sep = "/") |> strip_double_slashes()

  # CHM
  chm_path_thumb = file.path(dataset_id, processed_folder, "thumbnails", paste0(dataset_id, "_chm-mesh.png"))
  chm_path_full = file.path(dataset_id, processed_folder, "full", paste0(dataset_id, "_chm-mesh.tif"))
  chm_exists = chm_path_thumb %in% processed_products$filepath
  chm_url_thumb = paste(data_server_base_url, chm_path_thumb, sep = "/") |> strip_double_slashes()
  chm_url_full = paste(data_server_base_url, chm_path_full, sep = "/") |> strip_double_slashes()


  # DSM
  dsm_path_thumb = file.path(dataset_id, processed_folder, "thumbnails", paste0(dataset_id, "_dsm-mesh.png"))
  dsm_path_full = file.path(dataset_id, processed_folder, "full", paste0(dataset_id, "_dsm-mesh.tif"))
  dsm_exists = dsm_path_thumb %in% processed_products$filepath
  dsm_url_thumb = paste(data_server_base_url, dsm_path_thumb, sep = "/") |> strip_double_slashes()
  dsm_url_full = paste(data_server_base_url, dsm_path_full, sep = "/") |> strip_double_slashes()

  # DTM
  dtm_path_thumb = file.path(dataset_id, processed_folder, "thumbnails", paste0(dataset_id, "_dtm-ptcloud.png"))
  dtm_path_full = file.path(dataset_id, processed_folder, "full", paste0(dataset_id, "_dtm-ptcloud.tif"))
  dtm_exists = dtm_path_thumb %in% processed_products$filepath
  dtm_url_thumb = paste(data_server_base_url, dtm_path_thumb, sep = "/") |> strip_double_slashes()
  dtm_url_full = paste(data_server_base_url, dtm_path_full, sep = "/") |> strip_double_slashes()

  # Point cloud
  pc_path_full = file.path(dataset_id, processed_folder, "full", paste0(dataset_id, "_points.copc.laz"))
  pc_exists = pc_path_full %in% processed_products$filepath
  pc_url_full = paste(data_server_base_url, pc_path_full, sep = "/") |> strip_double_slashes()

  # Mesh model
  mesh_path_full = file.path(dataset_id, processed_folder, "full", paste0(dataset_id, "_mesh.ply"))
  mesh_exists = mesh_path_full %in% processed_products$filepath
  mesh_url_full = paste(data_server_base_url, mesh_path_full, sep = "/") |> strip_double_slashes()

  # Raw images
  image_example_path = file.path(dataset_id, "images", "examples", "thumbnails", "example_4.JPG")
  images_example_exists = image_example_path %in% s3_file_listing_foc$filepath
  images_example_url_thumb = paste(data_server_base_url, dataset_id, "images/examples/thumbnails/", sep = "/") |> strip_double_slashes()
  images_example_url_full = paste(data_server_base_url, dataset_id, "images/examples/fullsize/", sep = "/") |> strip_double_slashes()

  image_zip_path =  file.path(dataset_id, "images", paste0(dataset_id, "_images.zip"))
  images_zip_exists = image_zip_path %in% s3_file_listing_foc$filepath
  images_zip_url = paste(data_server_base_url, image_zip_path, sep = "/") |> strip_double_slashes()

  # Image metadata
  image_metadata_path = file.path(dataset_id, "metadata-images", paste0(dataset_id, "_image-metadata.gpkg"))
  image_metadata_exists = image_metadata_path %in% s3_file_listing_foc$filepath
  image_metadata_url = paste(data_server_base_url, image_metadata_path, sep = "/") |> strip_double_slashes()

  # Mission footprint
  footprint_path = file.path(dataset_id, "metadata-mission", paste0(dataset_id, "_mission-metadata.gpkg"))
  footprint_exists = footprint_path %in% s3_file_listing_foc$filepath
  footprint_url = paste(data_server_base_url, footprint_path, sep = "/") |> strip_double_slashes()

  # Cameras
  cameras_path = file.path(dataset_id, processed_folder, "full", paste0(dataset_id, "_cameras.xml"))
  cameras_exists = cameras_path %in% s3_file_listing_foc$filepath
  cameras_url = paste(data_server_base_url, cameras_path, sep = "/") |> strip_double_slashes()

  # Log
  log_path = file.path(dataset_id, processed_folder, "full", paste0(dataset_id, "_log.txt"))
  log_exists = log_path %in% s3_file_listing_foc$filepath
  log_url = paste(data_server_base_url, log_path, sep = "/") |> strip_double_slashes()

  # ITD
  ttops_exists = ttops_file_path %in% processed_products$filepath
  ttops_url = paste(data_server_base_url, ttops_file_path, sep = "/") |> strip_double_slashes()

  # Return a named list of all the computed values
  return(list(
    ortho_exists = ortho_exists,
    ortho_url_thumb = ortho_url_thumb,
    ortho_url_full = ortho_url_full,
    chm_exists = chm_exists,
    chm_url_thumb = chm_url_thumb,
    chm_url_full = chm_url_full,
    dsm_exists = dsm_exists,
    dsm_url_thumb = dsm_url_thumb,
    dsm_url_full = dsm_url_full,
    dtm_exists = dtm_exists,
    dtm_url_thumb = dtm_url_thumb,
    dtm_url_full = dtm_url_full,
    pc_exists = pc_exists,
    pc_url_full = pc_url_full,
    mesh_exists = mesh_exists,
    mesh_url_full = mesh_url_full,
    images_example_exists = images_example_exists,
    images_example_url_thumb = images_example_url_thumb,
    images_example_url_full = images_example_url_full,
    images_zip_exists = images_zip_exists,
    images_zip_url = images_zip_url,
    image_metadata_exists = image_metadata_exists,
    image_metadata_url = image_metadata_url,
    footprint_exists = footprint_exists,
    footprint_url = footprint_url,
    cameras_exists = cameras_exists,
    cameras_url = cameras_url,
    log_exists = log_exists,
    log_url = log_url,
    ttops_exists = ttops_exists,
    ttops_url = ttops_url,
    itd_path_mostrecent = itd_path_mostrecent
  ))
}


# Generic function to render dataset details page (works for missions, composites, etc.)
# Reduces code duplication between render_mission_details_page and render_composite_details_page
render_dataset_details_page = function(
    template_filepath,
    dataset_summary_foc,
    s3_file_listing_foc,
    data_server_base_url,
    details_map_path,
    itd_map_path,
    details_datatable_path,
    next_dataset_page_path,
    previous_dataset_page_path,
    website_repo_content_path,
    details_page_dir,
    display_data = FALSE,
    additional_jinjar_vars = list(),  # For dataset-type-specific variables
    # Optional post-curation parameters for curation pages (side-by-side display)
    has_post_curation_data = FALSE,
    post_curation_map_html_path = NA,
    post_curation_datatable_html_path = NA
  ) {

  dataset_id = unique(dataset_summary_foc$dataset_id)

  if(length(dataset_id) != 1) {
    stop("dataset_summary_foc should have exactly one unique dataset_id. Found: ", paste(dataset_id, collapse = ", "))
  }

  # Compute S3 product URLs if display_data is enabled, otherwise initialize to defaults
  if (display_data) {
    product_urls = compute_s3_product_urls(dataset_id, s3_file_listing_foc, data_server_base_url)
    ortho_exists = product_urls$ortho_exists
    ortho_url_thumb = product_urls$ortho_url_thumb
    ortho_url_full = product_urls$ortho_url_full
    chm_exists = product_urls$chm_exists
    chm_url_thumb = product_urls$chm_url_thumb
    chm_url_full = product_urls$chm_url_full
    dsm_exists = product_urls$dsm_exists
    dsm_url_thumb = product_urls$dsm_url_thumb
    dsm_url_full = product_urls$dsm_url_full
    dtm_exists = product_urls$dtm_exists
    dtm_url_thumb = product_urls$dtm_url_thumb
    dtm_url_full = product_urls$dtm_url_full
    pc_exists = product_urls$pc_exists
    pc_url_full = product_urls$pc_url_full
    mesh_exists = product_urls$mesh_exists
    mesh_url_full = product_urls$mesh_url_full
    images_example_exists = product_urls$images_example_exists
    images_example_url_thumb = product_urls$images_example_url_thumb
    images_example_url_full = product_urls$images_example_url_full
    images_zip_exists = product_urls$images_zip_exists
    images_zip_url = product_urls$images_zip_url
    image_metadata_exists = product_urls$image_metadata_exists
    image_metadata_url = product_urls$image_metadata_url
    footprint_exists = product_urls$footprint_exists
    footprint_url = product_urls$footprint_url
    cameras_exists = product_urls$cameras_exists
    cameras_url = product_urls$cameras_url
    log_exists = product_urls$log_exists
    log_url = product_urls$log_url
    ttops_exists = product_urls$ttops_exists
    ttops_url = product_urls$ttops_url
    itd_path_mostrecent = product_urls$itd_path_mostrecent
  } else {
    # Initialize all to defaults when not displaying data
    ortho_exists = FALSE
    ortho_url_thumb = NULL
    ortho_url_full = NULL
    chm_exists = FALSE
    chm_url_thumb = NULL
    chm_url_full = NULL
    dsm_exists = FALSE
    dsm_url_thumb = NULL
    dsm_url_full = NULL
    dtm_exists = FALSE
    dtm_url_thumb = NULL
    dtm_url_full = NULL
    pc_exists = FALSE
    pc_url_full = NULL
    mesh_exists = FALSE
    mesh_url_full = NULL
    images_example_exists = FALSE
    images_example_url_thumb = NULL
    images_example_url_full = NULL
    images_zip_exists = FALSE
    images_zip_url = NULL
    image_metadata_exists = FALSE
    image_metadata_url = NULL
    footprint_exists = FALSE
    footprint_url = NULL
    cameras_exists = FALSE
    cameras_url = NULL
    log_exists = FALSE
    log_url = NULL
    ttops_exists = FALSE
    ttops_url = NULL
    itd_path_mostrecent = NULL
  }

  # Build base jinjar variables (common to all dataset types)
  base_vars = list(
    dataset_id = dataset_id,
    map_html_path = details_map_path,
    itd_map_html_path = itd_map_path,
    itd_map_path = itd_map_path,  # Composite template uses this name
    datatable_html_path = details_datatable_path,
    next_dataset_page_path = next_dataset_page_path,
    previous_dataset_page_path = previous_dataset_page_path,
    ortho_exists = ortho_exists,
    ortho_url_thumb = ortho_url_thumb,
    ortho_url_full = ortho_url_full,
    chm_exists = chm_exists,
    chm_url_thumb = chm_url_thumb,
    chm_url_full = chm_url_full,
    dsm_exists = dsm_exists,
    dsm_url_thumb = dsm_url_thumb,
    dsm_url_full = dsm_url_full,
    dtm_exists = dtm_exists,
    dtm_url_thumb = dtm_url_thumb,
    dtm_url_full = dtm_url_full,
    pc_exists = pc_exists,
    pc_url_full = pc_url_full,
    mesh_exists = mesh_exists,
    mesh_url_full = mesh_url_full,
    images_example_exists = images_example_exists,
    images_example_url_thumb = images_example_url_thumb,
    images_example_url_full = images_example_url_full,
    images_zip_exists = images_zip_exists,
    images_zip_url = images_zip_url,
    image_metadata_exists = image_metadata_exists,
    image_metadata_url = image_metadata_url,
    footprint_exists = footprint_exists,
    footprint_url = footprint_url,
    cameras_exists = cameras_exists,
    cameras_url = cameras_url,
    log_exists = log_exists,
    log_url = log_url,
    ttops_exists = ttops_exists,
    ttops_url = ttops_url,
    itd_path_mostrecent = itd_path_mostrecent,
    has_post_curation_data = has_post_curation_data,
    post_curation_map_html_path = post_curation_map_html_path,
    post_curation_datatable_html_path = post_curation_datatable_html_path
  )

  # Merge with additional dataset-type-specific variables
  all_vars = c(base_vars, additional_jinjar_vars)
  all_vars$.config = jinjar_config(variable_open = "{*", variable_close = "*}")

  # Render template
  rendered = do.call(jinjar::render, c(list(template_filepath), all_vars))

  write_path = file.path(website_repo_content_path,
                         details_page_dir,
                         paste0(dataset_id, ".md"))
  writeLines(rendered, write_path)

  return(write_path)
}


# Compose and render the {dataset_id}.md page for a mission based on the Jinjar template
# This is now a wrapper around render_dataset_details_page for backward compatibility
render_mission_details_page = function(
    template_filepath,
    mission_summary_foc,
    s3_file_listing_foc,
    mission_details_map_path,
    itd_map_path,
    mission_details_datatable_path,
    next_dataset_page_path,
    previous_dataset_page_path,
    website_repo_content_path,
    mission_details_page_dir,
    display_data = FALSE,
    # Optional post-curation parameters for curation pages (side-by-side display)
    has_post_curation_data = FALSE,
    post_curation_map_html_path = NA,
    post_curation_datatable_html_path = NA
  ) {

  # Determine if this is an oblique mission for the template
  oblique = abs(as.numeric(mission_summary_foc$camera_pitch_derived)) > 10

  # Call generic function with mission-specific parameters
  render_dataset_details_page(
    template_filepath = template_filepath,
    dataset_summary_foc = mission_summary_foc,
    s3_file_listing_foc = s3_file_listing_foc,
    data_server_base_url = DATA_SERVER_MISSIONS_BASE_URL,
    details_map_path = mission_details_map_path,
    itd_map_path = itd_map_path,
    details_datatable_path = mission_details_datatable_path,
    next_dataset_page_path = next_dataset_page_path,
    previous_dataset_page_path = previous_dataset_page_path,
    website_repo_content_path = website_repo_content_path,
    details_page_dir = mission_details_page_dir,
    display_data = display_data,
    additional_jinjar_vars = list(oblique = oblique),
    has_post_curation_data = has_post_curation_data,
    post_curation_map_html_path = post_curation_map_html_path,
    post_curation_datatable_html_path = post_curation_datatable_html_path
  )

  return(TRUE)

}


## Loop through each mission and make a details page, including its media (map and datatable). TODO: There are some constants used in here that shoud be getting passed in as arguments from where this function is called
make_mission_details_page = function(
    mission_id_foc,
    all_mission_ids, # Needed for making the next and previous links
    mission_summaries,
    mission_points,
    s3_file_listing,
    website_static_path,
    website_content_path,
    leaflet_header_files_dir,
    datatable_header_files_dir,
    mission_details_datatable_dir,
    mission_details_map_dir,
    itd_map_dir,
    mission_details_template_filepath,
    mission_details_page_dir
  ) {

    print(paste0("Making mission details page for ", mission_id_foc))

    # Get the S3 file listing for this mission (so we know what data products exist)
    s3_file_listing_foc = s3_file_listing |> filter(mission_id == mission_id_foc)

    # Extract the mission-level metadata that's associated with that dataset
    mission_summary_foc = mission_summaries |> filter(mission_id == mission_id_foc)

    # Get the mission points for this mission
    mission_points_foc = mission_points |> filter(mission_id == mission_id_foc)

    # TODO: Instead of the above, could pull from S3 or the database directly on the fly.

    # Make details map and datatable
    mission_details_map_path = make_mission_details_map(
      mission_summary_foc = mission_summary_foc,
      mission_points_foc = mission_points_foc,
      mission_polygons_for_mission_details_map = mission_summary,
      mission_centroids = mission_centroids,
      website_static_path = website_static_path,
      leaflet_header_files_dir = leaflet_header_files_dir,
      mission_details_map_dir = mission_details_map_dir
    )

    mission_details_datatable_path = make_mission_details_datatable(
      mission_summary_foc = mission_summary_foc,
      website_static_path = website_static_path,
      datatable_header_files_dir = datatable_header_files_dir,
      mission_details_datatable_dir = mission_details_datatable_dir
    )
    
    # Make detected tree map, if ITD data exists. For now using the most recent ITD folder if there
    # are multiple.
    
    processed_products = s3_file_listing_foc |> filter(str_detect(filepath, paste0("processed_", PHOTOGRAMMETRY_CONFIG_ID)))

    filepath_parts = str_split(processed_products$filepath, fixed("/"))
    part_3 = purrr::map_chr(filepath_parts, 3)
    itd_folders = part_3[str_which(part_3, "^(itd_|detected_trees_|detected-trees_)")] |> unique() |> sort(decreasing = TRUE)
    itd_folder_mostrecent = itd_folders[1]
    itd_path_mostrecent = file.path(mission_id_foc, paste0("processed_", PHOTOGRAMMETRY_CONFIG_ID), itd_folder_mostrecent)
    ttops_file_path = file.path(itd_path_mostrecent, paste0(mission_id_foc, "_treetops.gpkg"))
    ttops_exists = ttops_file_path %in% processed_products$filepath

    # TODO: We could streamline the above. For example, we could just check if any path with the
    # file name above exists, and if so, sort them and keep the alphabetically last (which should be
    # the most recent).

    if (ttops_exists) {
      ttops_url = paste(DATA_SERVER_MISSIONS_BASE_URL, ttops_file_path, sep = "/") |> strip_double_slashes()

      # Download ttops from S3
      remote_file = paste0(RCLONE_REMOTE, ":", REMOTE_MISSIONS_DIR, ttops_file_path)
      temp_ttops_file = tempfile(paste0("ttops", mission_id_foc), fileext = ".gpkg")
      command = paste("rclone copyto", remote_file, temp_ttops_file, "--progress --transfers 32 --checkers 32 --stats 1s --retries 5 --retries-sleep=15s --s3-upload-cutoff 100Mi --s3-chunk-size 100Mi --s3-upload-concurrency 16 --multi-thread-streams 2", sep = " ")
      system(command)

      itd_points = st_read(temp_ttops_file)

      # Make ITD map
      itd_map_path = make_itd_map(
        mission_summary_foc = mission_summary_foc,
        itd_points_foc = itd_points,
        mission_polygons_for_mission_details_map = mission_summary,
        mission_centroids = mission_centroids,
        website_static_path = website_static_path,
        leaflet_header_files_dir = leaflet_header_files_dir,
        itd_map_dir = itd_map_dir
      )
    } else {
      itd_map_path = NA
    }



    # Compute previous and next dataset. We have the current
    # mission_id_foc the list of all_mission_ids (assume it is sorted) and we need to find the previous
    # and next, wrapping around as needed.
    current_index = which(all_mission_ids == mission_id_foc)
    if (current_index == 1) {
      previous_mission_id = all_mission_ids[length(all_mission_ids)]
    } else {
      previous_mission_id = all_mission_ids[current_index - 1]
    }
    if (current_index == length(all_mission_ids)) {
      next_mission_id = all_mission_ids[1]
    } else {
      next_mission_id = all_mission_ids[current_index + 1]
    }

    next_dataset_page_path = paste0("/", mission_details_page_dir, "/", next_mission_id)
    previous_dataset_page_path = paste0("/", mission_details_page_dir, "/", previous_mission_id)

    # Render plot details page from template
    render_mission_details_page(
      template_filepath = mission_details_template_filepath,
      mission_summary_foc = mission_summary_foc,
      s3_file_listing = s3_file_listing_foc,
      mission_details_map_path = mission_details_map_path,
      itd_map_path = itd_map_path,
      mission_details_datatable_path = mission_details_datatable_path,
      next_dataset_page_path = next_dataset_page_path,
      previous_dataset_page_path = previous_dataset_page_path,
      website_repo_content_path = website_content_path,
      mission_details_page_dir = mission_details_page_dir,
      display_data = TRUE
    )


  gc()

}


# ==== COMPOSITE MISSION FUNCTIONS ====


# Compile composite summary data for catalog generation
compile_composite_summary_data = function(composite_mission_metadata,
                                          base_ofo_url,
                                          composite_details_dir,
                                          individual_mission_details_dir) {

  # First, call existing compile_mission_summary_data to get standard derived columns
  d = compile_mission_summary_data(composite_mission_metadata, base_ofo_url, individual_mission_details_dir, dataset_type = "mission")

  # Override dataset_id with composite_id
  d$dataset_id = d$composite_id

  # Override dataset_id_link to point to composite detail pages
  d$dataset_id_link = paste0('<a href="', base_ofo_url, composite_details_dir, d$dataset_id, '/"', ' target="_PARENT">', d$dataset_id, "</a>")

  # Add mission link (mission_type already has "hn" or "lo")
  d = d |>
    group_by(composite_id) |>
    mutate(
      individual_mission_link = paste0('<a href="', base_ofo_url, individual_mission_details_dir, mission_id, '/"', ' target="_PARENT">', mission_id, "</a>")
    ) |>
    arrange(composite_id, desc(altitude_agl_nominal)) |>  # Higher altitude first within each composite
    ungroup()

  return(d)
}


# Helper function to collapse composite fields
# If values are the same for both datasets, use the single value
# If different, paste higher and lower values with separator
collapse_composite_field = function(field, mission_type) {
  # mission_type is "hn" (high nadir) or "lo" (low oblique)
  # if (5.63 %in% field) browser()
  if_else(
    n_distinct(field) == 1,
    first(field),
    paste(field[mission_type == "hn"], field[mission_type == "lo"], sep = "   |   ")
  )
}


# Make composite catalog datatable with separator for differing values
make_composite_catalog_datatable = function(composite_summary,
                                            website_static_path,
                                            datatable_header_files_dir,
                                            composite_catalog_datatable_dir,
                                            composite_catalog_datatable_filename) {

  # Group by composite_id and collapse values
  d = composite_summary |>
    sf::st_drop_geometry() |>
    group_by(composite_id) |>
    summarise(
      dataset_id_link = first(dataset_id_link),  # Link is same for both
      area_derived = collapse_composite_field(area_derived, mission_type),
      earliest_date_derived = collapse_composite_field(earliest_date_derived, mission_type),
      altitude_agl_nominal = collapse_composite_field(altitude_agl_nominal, mission_type),
      overlap_front_side_nominal = collapse_composite_field(overlap_front_side_nominal, mission_type),
      camera_pitch_derived = collapse_composite_field(camera_pitch_derived, mission_type),
      terrain_follow = collapse_composite_field(terrain_follow, mission_type),
      flight_pattern = collapse_composite_field(flight_pattern, mission_type),
      image_count_derived = collapse_composite_field(image_count_derived, mission_type),
      rtk_nominal = collapse_composite_field(rtk_nominal, mission_type),
      percent_images_rtk_derived = collapse_composite_field(percent_images_rtk_derived, mission_type),
      contributor_dataset_name = collapse_composite_field(contributor_dataset_name, mission_type),
      project_id = collapse_composite_field(project_id, mission_type),
      aircraft_model_name = collapse_composite_field(aircraft_model_name, mission_type)
    ) |>
    ungroup() |>
    select("ID" = dataset_id_link,
           "Area (ha)" = area_derived,
           "Date" = earliest_date_derived,
           "Altitude (m) (N)" = altitude_agl_nominal,
           "Overlap (N)" = overlap_front_side_nominal,
           "Camera pitch" = camera_pitch_derived,
           "Terrain follow (N)" = terrain_follow,
           "Flight pattern" = flight_pattern,
           "Image count" = image_count_derived,
           "RTK (N)" = rtk_nominal,
           "RTK images (%)" = percent_images_rtk_derived,
           "Contributor dataset name" = contributor_dataset_name,
           "Project" = project_id,
           "Aircraft" = aircraft_model_name) |>
    arrange(ID)

  # Prep formatting code to pass to datatable creation
  format_js = DT::JS("function(settings, json) {",
                     "$('body').css({'font-family': 'Arial'});",
                     "}")

  dt = DT::datatable(d,
                     rownames = FALSE,
                     escape = FALSE,
                     options = list(paging = FALSE, initComplete = format_js))

  # Save the datatable HTML to the website repo
  save_widget_html(dt,
                   website_static_path = website_static_path,
                   header_files_dir = datatable_header_files_dir,
                   html_dir = composite_catalog_datatable_dir,
                   html_filename = composite_catalog_datatable_filename,
                   delete_folder_first = TRUE)

  return(dt)
}


# Make composite catalog map showing one polygon per composite (the higher-altitude one)
make_composite_catalog_map = function(composite_summary,
                                      website_static_path,
                                      leaflet_header_files_dir,
                                      composite_catalog_map_dir,
                                      composite_catalog_map_filename) {

  # Keep only the high nadir (hn) mission from each composite for the map
  composite_summary_map = composite_summary |>
    filter(mission_type == "hn")

  # Get both rows for popup content
  composite_summary_popup = composite_summary |>
    group_by(composite_id) |>
    summarise(
      dataset_id_link = first(dataset_id_link),
      popup = paste0(
        "<b>Composite: </b>", first(dataset_id_link), "<br>",
        "<b>Mission ", mission_id[mission_type == "hn"], " (high nadir):</b><br>",
        "&nbsp;&nbsp;Date: ", earliest_date_derived[mission_type == "hn"], "<br>",
        "&nbsp;&nbsp;Alt: ", altitude_agl_nominal[mission_type == "hn"], " m<br>",
        "&nbsp;&nbsp;Overlap: ", overlap_front_side_nominal[mission_type == "hn"], "<br>",
        "&nbsp;&nbsp;Pitch: ", camera_pitch_derived[mission_type == "hn"], " deg<br>",
        "<b>Mission ", mission_id[mission_type == "lo"], " (low oblique):</b><br>",
        "&nbsp;&nbsp;Date: ", earliest_date_derived[mission_type == "lo"], "<br>",
        "&nbsp;&nbsp;Alt: ", altitude_agl_nominal[mission_type == "lo"], " m<br>",
        "&nbsp;&nbsp;Overlap: ", overlap_front_side_nominal[mission_type == "lo"], "<br>",
        "&nbsp;&nbsp;Pitch: ", camera_pitch_derived[mission_type == "lo"], " deg<br>"
      ),
      # Dynamically use the actual geometry column name (could be "geometry" or "geom")
      !!attr(composite_summary, "sf_column") := first(!!sym(attr(composite_summary, "sf_column")))
    ) |>
    ungroup() |>
    st_as_sf()

  composite_centroids = sf::st_centroid(composite_summary_popup)

  m = leaflet() |>
    addTiles(options = providerTileOptions(maxZoom = 16)) |>
    addMarkers(data = composite_centroids, popup = ~popup, clusterOptions = markerClusterOptions(freezeAtZoom = 16)) |>
    addPolygons(data = composite_summary_popup, popup = ~popup, group = "bounds") |>
    addProviderTiles(providers$Esri.WorldTopoMap, group = "Topo", options = providerTileOptions(maxZoom = 16)) |>
    addProviderTiles(providers$Esri.WorldImagery, group = "Imagery") |>
    groupOptions("bounds", zoomLevels = 13:20) |>
    addLayersControl(baseGroups = c("Topo", "Imagery"),
                     options = layersControlOptions(collapsed = FALSE))

  save_widget_html(m,
                   website_static_path = website_static_path,
                   header_files_dir = leaflet_header_files_dir,
                   html_dir = composite_catalog_map_dir,
                   html_filename = composite_catalog_map_filename,
                   delete_folder_first = TRUE)

  return(m)
}


# Make composite details map with both polygons and all image points, color-coded by mission_id
make_composite_details_map = function(composite_summary_foc,
                                      composite_points_foc,
                                      website_static_path,
                                      leaflet_header_files_dir,
                                      composite_details_map_dir) {

  composite_id = unique(composite_summary_foc$composite_id)
  mission_ids = composite_summary_foc$mission_id

  # Create color palette for two missions
  mission_colors = colorFactor(c("#E41A1C", "#377EB8"), domain = mission_ids)

  # Create flight paths per mission
  flight_paths = composite_points_foc |>
    arrange(mission_id, datetime_local) |>
    group_by(mission_id) |>
    summarise(do_union = FALSE) |>
    st_cast("LINESTRING")

  # Create map
  m = leaflet() |>
    addProviderTiles(providers$Esri.WorldTopoMap, group = "Topo", options = providerTileOptions(maxZoom = 22)) |>
    addProviderTiles(providers$Esri.WorldImagery, group = "Imagery", options = providerTileOptions(maxZoom = 22))

  # Add polygons, flight paths, and image points per mission (grouped together for per-mission toggling)
  for (i in 1:nrow(composite_summary_foc)) {
    mission = composite_summary_foc[i, ]
    mission_id_i = mission$mission_id
    group_name = as.character(mission_id_i)

    # Footprint
    m = m |>
      addPolygons(data = mission,
                  fillOpacity = 0,
                  color = mission_colors(mission_id_i),
                  weight = 2,
                  group = group_name)

    # Flight path
    path = flight_paths[flight_paths$mission_id == mission_id_i, ]
    if (nrow(path) > 0) {
      m = m |>
        addPolylines(data = path,
                     color = mission_colors(mission_id_i),
                     weight = 1,
                     opacity = 0.7,
                     group = group_name)
    }

    # Image points (visible) + larger invisible markers on top for easier clicking
    points_i = composite_points_foc |> filter(mission_id == mission_id_i)
    m = m |>
      addCircleMarkers(data = points_i,
                       radius = 3,
                       fillColor = mission_colors(mission_id_i),
                       fillOpacity = 0.6,
                       color = mission_colors(mission_id_i),
                       weight = 1,
                       group = group_name) |>
      addCircleMarkers(data = points_i,
                       radius = 10,
                       fillOpacity = 0,
                       stroke = FALSE,
                       group = group_name,
                       popup = ~image_id)
  }

  # Add legend
  m = m |>
    addLegend("bottomright",
              colors = c("#E41A1C", "#377EB8"),
              labels = mission_ids,
              title = "Mission ID",
              opacity = 1) |>
    addLayersControl(baseGroups = c("Topo", "Imagery"),
                     overlayGroups = as.character(mission_ids),
                     options = layersControlOptions(collapsed = FALSE))

  # Save widget
  composite_details_map_filename = paste0(composite_id, ".html")
  save_widget_html(m,
                   website_static_path = website_static_path,
                   header_files_dir = leaflet_header_files_dir,
                   html_dir = composite_details_map_dir,
                   html_filename = composite_details_map_filename,
                   delete_folder_first = FALSE)

  composite_details_map_path = paste(composite_details_map_dir, composite_details_map_filename, sep = "/")

  return(composite_details_map_path)
}


# Make composite details datatable with 3 columns: Attribute, Mission A (higher), Mission B (lower)
make_composite_details_datatable = function(composite_summary_foc,
                                            website_static_path,
                                            datatable_header_files_dir,
                                            composite_details_datatable_dir) {

  composite_id = unique(composite_summary_foc$composite_id)

  # Ensure higher altitude (HN) is first
  composite_summary_foc = composite_summary_foc |>
    arrange(mission_type)  # "hn" is higher and comes first alphabetically

  mission_id_a = composite_summary_foc$mission_id[1]
  mission_id_b = composite_summary_foc$mission_id[2]

  # Select same attributes as individual mission details datatable, plus composite-specific columns
  d = composite_summary_foc |>
    sf::st_drop_geometry() |>
    dplyr::select(
      mission_id,
      individual_mission_link,
      # sub_mission_ids,
      earliest_date_derived,
      time_range_local_derived,
      altitude_agl_nominal,
      terrain_follow,
      camera_pitch_derived,
      smart_oblique_derived,
      overlap_front_side_nominal,
      flight_pattern,
      any_of(c(# "altitude_agl_mean_derived", <- I blieve this is just an obsolete placeholder that usess the old name we have deprecated
              #  "photogrammetry_altitude_agl_median_derived",
               "photogrammetry_altitude_agl_mean_derived",
               "photogrammetry_terrain_fidelity_derived")),
      gps_terrain_fidelity_derived,
      percent_images_rtk_derived,
      flight_speed_derived,
      exposure_median_derived,
      exposure_cv_derived,
      exposure_compensation,
      white_balance_mode_derived,
      white_balance_pct_mode_derived,
      aircraft_model_name,
      sensor_name,
      flight_planner_name,
      base_lat,
      base_lon,
      base_alt,
      base_marked_permanently,
      area_derived,
      image_count_derived,
      image_dimensions_derived,
      file_size_derived,
      file_format_derived,
      project_id,
      contributor_dataset_name,
      contributor_names,
      license
    )

  # Pivot to long format for each mission
  d_long_a = d |>
    filter(mission_id == mission_id_a) |>
    select(-mission_id) |>
    pivot_longer(everything(), names_to = "Attribute", values_to = "High nadir mission")

  d_long_b = d |>
    filter(mission_id == mission_id_b) |>
    select(-mission_id) |>
    pivot_longer(everything(), names_to = "Attribute", values_to = "Low oblique mission")

  # Join on Attribute
  d_combined = d_long_a |>
    left_join(d_long_b, by = "Attribute")

  # Format attribute names for display
  d_combined = d_combined |>
    mutate(Attribute = case_when(
      Attribute == "individual_mission_link" ~ "Mission ID",
      # Attribute == "sub_mission_ids" ~ "Sub-mission IDs",
      Attribute == "earliest_date_derived" ~ "Date",
      Attribute == "time_range_local_derived" ~ "Flight time range (local)",
      Attribute == "altitude_agl_nominal" ~ "Altitude (nominal) (m)",
      Attribute == "altitude_agl_mean_derived" ~ "Photogrammetry altitude mean 2 (m)",
      # Attribute == "photogrammetry_altitude_agl_median_derived" ~ "Photogrammetry altitude median (m)",
      Attribute == "photogrammetry_altitude_agl_mean_derived" ~ "Photogrammetry altitude mean (m)",
      Attribute == "photogrammetry_terrain_fidelity_derived" ~ "Photogrammetry terrain fidelity",
      Attribute == "flight_pattern" ~ "Flight pattern",
      Attribute == "overlap_front_side_nominal" ~ "Overlap (nominal) (front/side)",
      Attribute == "camera_pitch_derived" ~ "Camera pitch (deg up from nadir)",
      Attribute == "smart_oblique_derived" ~ "Smart oblique",
      Attribute == "terrain_follow" ~ "Terrain follow (nominal)",
      Attribute == "gps_terrain_fidelity_derived" ~ "GPS terrain fidelity",
      Attribute == "percent_images_rtk_derived" ~ "Percent RTK",
      Attribute == "flight_speed_derived" ~ "Flight speed (m/s)",
      Attribute == "exposure_median_derived" ~ "Exposure median (sec)",
      Attribute == "exposure_cv_derived" ~ "Exposure CV",
      Attribute == "exposure_compensation" ~ "Exposure compensation (nominal)",
      Attribute == "white_balance_mode_derived" ~ "White balance mode",
      Attribute == "white_balance_pct_mode_derived" ~ "White balance percent",
      Attribute == "aircraft_model_name" ~ "Aircraft model",
      Attribute == "sensor_name" ~ "Sensor model",
      Attribute == "flight_planner_name" ~ "Flight planner",
      Attribute == "base_lat" ~ "Base station latitude",
      Attribute == "base_lon" ~ "Base station longitude",
      Attribute == "base_alt" ~ "Base station altitude (m)",
      Attribute == "base_marked_permanently" ~ "Permanent base marker",
      Attribute == "area_derived" ~ "Flight footprint area (ha)",
      Attribute == "image_count_derived" ~ "Image count",
      Attribute == "image_dimensions_derived" ~ "Image dimensions",
      Attribute == "file_size_derived" ~ "Dataset size (GB)",
      Attribute == "file_format_derived" ~ "File format",
      Attribute == "project_id" ~ "Project ID",
      Attribute == "contributor_dataset_name" ~ "Contributor dataset name",
      Attribute == "contributor_names" ~ "Creator",
      Attribute == "license" ~ "License",
      TRUE ~ Attribute
    ))

  # Remove "Attribute" header from first column
  d_combined = d_combined |>
    rename(" " = Attribute)

  # Create datatable
  format_js = DT::JS("function(settings, json) {",
                     "$('body').css({'font-family': 'Arial'});",
                     "}")

  dt = DT::datatable(d_combined,
                     rownames = FALSE,
                     escape = FALSE,
                     options = list(
                       paging = FALSE,
                       dom = 't',
                       bSort = FALSE,
                       initComplete = format_js,
                       columnDefs = list(
                         list(width = '35%', targets = 0),
                         list(width = '32.5%', targets = 1),
                         list(width = '32.5%', targets = 2)
                       )
                     ))

  dt$sizingPolicy$browser$padding = 0
  dt$sizingPolicy$browser$fill = FALSE

  # Save widget
  composite_details_datatable_filename = paste0(composite_id, ".html")
  save_widget_html(dt,
                   website_static_path = website_static_path,
                   header_files_dir = datatable_header_files_dir,
                   html_dir = composite_details_datatable_dir,
                   html_filename = composite_details_datatable_filename,
                   delete_folder_first = FALSE)

  composite_details_datatable_path = paste(composite_details_datatable_dir, composite_details_datatable_filename, sep = "/")

  return(composite_details_datatable_path)
}


# Render composite details page template
# This is now a wrapper around render_dataset_details_page
render_composite_details_page = function(
    template_filepath,
    composite_summary_foc,
    s3_file_listing_foc,
    composite_details_map_path,
    itd_map_path,
    composite_details_datatable_path,
    next_dataset_page_path,
    previous_dataset_page_path,
    website_repo_content_path,
    composite_details_page_dir,
    display_data = TRUE) {

  composite_id = unique(composite_summary_foc$composite_id)

  # Get mission IDs (higher altitude first)
  composite_summary_foc = composite_summary_foc |>
    arrange(mission_type)  # "hn" is higher and comes first alphabetically
  mission_id_a = composite_summary_foc$mission_id[1]
  mission_id_b = composite_summary_foc$mission_id[2]

  # Build individual mission page paths
  individual_mission_page_path_a = paste0("/", MISSION_DETAILS_PAGE_DIR, mission_id_a, "/")
  individual_mission_page_path_b = paste0("/", MISSION_DETAILS_PAGE_DIR, mission_id_b, "/")

  # Call generic function with composite-specific parameters
  render_dataset_details_page(
    template_filepath = template_filepath,
    dataset_summary_foc = composite_summary_foc,
    s3_file_listing_foc = s3_file_listing_foc,
    data_server_base_url = DATA_SERVER_COMPOSITES_BASE_URL,
    details_map_path = composite_details_map_path,
    itd_map_path = itd_map_path,
    details_datatable_path = composite_details_datatable_path,
    next_dataset_page_path = next_dataset_page_path,
    previous_dataset_page_path = previous_dataset_page_path,
    website_repo_content_path = website_repo_content_path,
    details_page_dir = composite_details_page_dir,
    display_data = display_data,
    additional_jinjar_vars = list(
      composite_id = composite_id,
      mission_id_a = mission_id_a,
      mission_id_b = mission_id_b,
      individual_mission_page_path_a = individual_mission_page_path_a,
      individual_mission_page_path_b = individual_mission_page_path_b,
      composite_details_map_path = composite_details_map_path,
      composite_details_datatable_path = composite_details_datatable_path
    )
  )
}


# Top-level orchestrator for a single composite's detail page
make_composite_details_page = function(composite_id_foc,
                                       all_composite_ids,
                                       composite_summaries,
                                       composite_points,
                                       s3_file_listing,
                                       website_static_path,
                                       website_content_path,
                                       datatable_header_files_dir,
                                       leaflet_header_files_dir,
                                       composite_details_template_filepath,
                                       composite_details_datatable_dir,
                                       composite_details_map_dir,
                                       composite_details_page_dir,
                                       composite_itd_map_dir) {

  # Filter to focal composite
  composite_summary_foc = composite_summaries |>
    filter(composite_id == composite_id_foc) |>
    arrange(mission_type)  # "hn" is higher and comes first alphabetically

  composite_points_foc = composite_points |>
    filter(composite_id == composite_id_foc)

  s3_file_listing_foc = s3_file_listing |>
    filter(composite_id == composite_id_foc)

  # Extract mission IDs
  mission_id_a = composite_summary_foc$mission_id[1]  # higher
  mission_id_b = composite_summary_foc$mission_id[2]  # lower

  # Make composite details map
  composite_details_map_path = make_composite_details_map(
    composite_summary_foc = composite_summary_foc,
    composite_points_foc = composite_points_foc,
    website_static_path = website_static_path,
    leaflet_header_files_dir = leaflet_header_files_dir,
    composite_details_map_dir = composite_details_map_dir
  )

  # Make composite details datatable
  composite_details_datatable_path = make_composite_details_datatable(
    composite_summary_foc = composite_summary_foc,
    website_static_path = website_static_path,
    datatable_header_files_dir = datatable_header_files_dir,
    composite_details_datatable_dir = composite_details_datatable_dir
  )

  # Check for ITD data and create ITD map if exists
  itd_map_path = NA
  product_urls = compute_s3_product_urls(composite_id_foc, s3_file_listing_foc, DATA_SERVER_COMPOSITES_BASE_URL)
  if (product_urls$ttops_exists) {
    # Download ITD data
    itd_path_mostrecent = product_urls$itd_path_mostrecent
    ttops_url = product_urls$ttops_url
    ttops_tempfile = tempfile(fileext = ".gpkg")

    download_result = tryCatch({
      download.file(ttops_url, ttops_tempfile, quiet = TRUE, method = "wget")
      TRUE
    }, error = function(e) {
      FALSE
    })

    if (download_result && file.exists(ttops_tempfile)) {
      ttops = st_read(ttops_tempfile, quiet = TRUE)

      # Create ITD map (simplified version without attribute switching)
      m_itd = leaflet() |>
        addProviderTiles(providers$Esri.WorldTopoMap, group = "Topo", options = providerTileOptions(maxZoom = 22)) |>
        addProviderTiles(providers$Esri.WorldImagery, group = "Imagery", options = providerTileOptions(maxZoom = 22)) |>
        addCircleMarkers(data = ttops,
                         radius = 3,
                         fillColor = "green",
                         fillOpacity = 0.6,
                         color = "darkgreen",
                         weight = 1) |>
        addLayersControl(baseGroups = c("Topo", "Imagery"),
                         options = layersControlOptions(collapsed = FALSE))

      itd_map_filename = paste0(composite_id_foc, ".html")
      save_widget_html(m_itd,
                       website_static_path = website_static_path,
                       header_files_dir = leaflet_header_files_dir,
                       html_dir = composite_itd_map_dir,
                       html_filename = itd_map_filename,
                       delete_folder_first = FALSE)

      itd_map_path = paste(composite_itd_map_dir, itd_map_filename, sep = "/")

      unlink(ttops_tempfile)
    }
  }

  # Compute prev/next composite navigation (wraparound)
  current_index = which(all_composite_ids == composite_id_foc)
  next_index = ifelse(current_index == length(all_composite_ids), 1, current_index + 1)
  previous_index = ifelse(current_index == 1, length(all_composite_ids), current_index - 1)
  next_dataset_page_path = paste0("/", composite_details_page_dir, all_composite_ids[next_index], "/")
  previous_dataset_page_path = paste0("/", composite_details_page_dir, all_composite_ids[previous_index], "/")

  # Render composite details page
  render_composite_details_page(
    template_filepath = composite_details_template_filepath,
    composite_summary_foc = composite_summary_foc,
    s3_file_listing = s3_file_listing_foc,
    composite_details_map_path = composite_details_map_path,
    itd_map_path = itd_map_path,
    composite_details_datatable_path = composite_details_datatable_path,
    next_dataset_page_path = next_dataset_page_path,
    previous_dataset_page_path = previous_dataset_page_path,
    website_repo_content_path = website_content_path,
    composite_details_page_dir = composite_details_page_dir,
    display_data = TRUE
  )

  gc()
}
