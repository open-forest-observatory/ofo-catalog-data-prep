# Purpose: Take the raw imagery folders and compute the "deliverable" versions of the outputs by
# postprocessing. Outputs include example images and a zipped folder of all images. Perform at the
# mission level.

# TODO: This relies on an attribute of the image-level metadata called "received_image_path" that
# really represents the path of the image in the raw, organized imagery folder (i.e., likely
# different than as it was "received"). Upstream, we should change the name of this attribute (or
# add a new one that we use here instead) that is more descriptive of this.

library(tidyverse)
library(sf)
library(magick)
library(furrr)

source("src/utils.R")

magick:::magick_threads(1)


# Function to do all the imagery prep for a given mission, with pre-subsetted metadata and footprint
make_raw_imagery_thumbnails_and_zip = function(mission_id_foc, use_post_curation = TRUE) {

  cat("\n **** Making thumbnails and zip for mission", mission_id_foc, "**** \n")

  # Skip if the mission already has all outputs, asuming that if the zip file exists, the entire
  # mission was processed to completion
  zip_outpath = file.path(IMAGERY_ZIP_AND_EXAMPLES_PATH, mission_id_foc, "images", paste0(mission_id_foc, "_images.zip"))

  if (SKIP_EXISTING && file.exists(zip_outpath)) {
    cat("Already exists. Skipping.\n")
    return()
  }

  if(!dir.exists(IN_PROCESS_PATH)) {
    dir.create(IN_PROCESS_PATH, recursive = TRUE)
  }

  # Save a file that indicates the mission is being processed
  processing_file = file.path(IN_PROCESS_PATH, paste0(mission_id_foc, ".csv"))
  fake_df = data.frame(a = 1, b = 1)
  write.csv(fake_df, processing_file)

  # Select appropriate metadata path
  if (use_post_curation) {
    points_filepath = file.path(POST_CURATION_PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
                                paste0(mission_id_foc, "_image-metadata.gpkg"))
    footprint_filepath = file.path(POST_CURATION_FULL_METADATA_PER_MISSION_PATH,
                                    paste0(mission_id_foc, "_mission-metadata.gpkg"))
  } else {
    points_filepath = file.path(PARSED_EXIF_FOR_RETAINED_IMAGES_PATH,
                                paste0(mission_id_foc, "_image-metadata.gpkg"))
    footprint_filepath = file.path(FULL_METADATA_PER_MISSION_PATH,
                                    paste0(mission_id_foc, "_mission-metadata.gpkg"))
  }

  # Get the mission images metadata (points gpkg)
  mission_images_metadata = st_read(points_filepath)

  # Get the mission footprint (polygon gpkg)
  mission_footprint = st_read(footprint_filepath)

  # Project mission image locs and footprint to the local UTM zone
  mission_images_metadata = transform_to_local_utm(mission_images_metadata)
  mission_footprint = transform_to_local_utm(mission_footprint)

  # In case one was in a different UTM zone than the other, reproject to the same zone
  mission_footprint = st_transform(mission_footprint, st_crs(mission_images_metadata))

  # Get the footprint area
  mission_footprint_area = mission_footprint |> st_area()

  # Buffer in the polygon by 100 m
  mission_footprint_buffered = mission_footprint |> st_buffer(-100)

  # Get the new area
  mission_footprint_buffered_area = mission_footprint_buffered |> st_area()

  # If the new area is < 40% of the original, use a smaller buffer, first trying 50 m, then 25 m,
  # then 0
  if (as.vector(mission_footprint_buffered_area / mission_footprint_area) < 0.4) {
    mission_footprint_buffered = mission_footprint |> st_buffer(-50)
    mission_footprint_buffered_area = mission_footprint_buffered |> st_area()
  }

  if (as.vector(mission_footprint_buffered_area / mission_footprint_area) < 0.4) {
    mission_footprint_buffered = mission_footprint |> st_buffer(-25)
    mission_footprint_buffered_area = mission_footprint_buffered |> st_area()
  }

  if (as.vector(mission_footprint_buffered_area / mission_footprint_area) < 0.4) {
    mission_footprint_buffered = mission_footprint
  }

  # Get the image locations within the buffered footprint
  mission_images_interior = mission_images_metadata |> st_intersection(mission_footprint_buffered)

  # Group points into 4 clusters
  coords = st_coordinates(mission_images_interior)
  clusts = kmeans(coords, centers = N_EXAMPLE_IMAGES)
  mission_images_interior$cluster = clusts$cluster

  # Get the cluster centroids, then the image point nearest each, to be the 4 representative images
  # from the mission
  selected_images = data.frame()
  for (i in 1:N_EXAMPLE_IMAGES) {
    cluster_pts = mission_images_interior |> filter(cluster == i)
    cluster_centroid = cluster_pts |> st_union() |> st_centroid()
    nearest_idx = st_nearest_feature(cluster_centroid, cluster_pts)
    nearest_pt = cluster_pts[nearest_idx, ]
    selected_images = rbind(selected_images, nearest_pt)
  }

  # Symlink the selected images to the publishable folder
  inpaths = file.path(SORTED_IMAGERY_PATH, selected_images$image_path_ofo)
  extensions = tools::file_ext(inpaths)
  outpaths = file.path(IMAGERY_ZIP_AND_EXAMPLES_PATH, mission_id_foc, "images", "examples", "fullsize", paste0("example_", 1:N_EXAMPLE_IMAGES, ".", extensions))

  outdirs = unique(dirname(outpaths))
  walk(outdirs, dir.create, recursive = TRUE)

  # Remove the file if it already exists
  exists = file.exists(outpaths)
  outpaths_remove = outpaths[exists]
  file.remove(outpaths_remove)

  file.symlink(inpaths, outpaths)

  # Create thumbnails of the images using magick package and save to the publishable folder
  for (i in 1:N_EXAMPLE_IMAGES) {
    img = image_read(inpaths[i])
    dim_string = paste0(THUMBNAIL_SIZE, "x", THUMBNAIL_SIZE)
    img = image_resize(img, dim_string)

    thumb_outpath = file.path(IMAGERY_ZIP_AND_EXAMPLES_PATH, mission_id_foc, "images", "examples", "thumbnails", paste0("example_", i, ".", extensions[i]))

    outdirs = unique(dirname(thumb_outpath))
    walk(outdirs, dir.create, recursive = TRUE, showWarnings = FALSE)

    image_write(img, thumb_outpath)
  }

  # Zip the entirety of the raw images folder and save to the publishable folder
  # Save to a tempfile while creating the zip, then move to the final location, because if the
  # process is terminated we don't want to leave a partial temp file in the file tree
  inpath = file.path(SORTED_IMAGERY_PATH, mission_id_foc)

  tempfile = file.path(TEMPDIR, "catalog-prep_tempzip", paste0("tempzip_", mission_id_foc, ".zip"))
  # Delete if exists
  if (file.exists(tempfile)) {
    file.remove(tempfile)
  }
  # Create dir
  dir.create(dirname(tempfile), recursive = TRUE, showWarnings = FALSE)

  command = paste("cd", shQuote(inpath), "; zip -r -0", shQuote(tempfile), ".")
  system(command, ignore.stdout = TRUE)
  file.rename(tempfile, zip_outpath)

  # Remove the file indicating the mission is being processed
  file.remove(processing_file)

  gc()

}
