# Purpose: For the specified mission ID, download the images from CyVerse and unzip them into a
# staging folder for photogrammetry processing

source("deploy/drone-imagery-ingestion/00_set-constants.R")


download_unzip_images = function(mission_id_foc) {

  # Destination folder:
  zip_filepath = file.path(
    PHOTOGRAMMETRY_DIR, DOWNLOADED_IMAGERY_ZIP_SUBDIR,
    paste0(mission_id_foc, "_images.zip")
  )

  # Make sure the destination exists
  dir.create(dirname(zip_filepath), recursive = TRUE, showWarnings = FALSE)


  # Option 1: Download over HTTP

  cyverse_url = paste0(
    "https://data.cyverse.org/dav-anon/",
    CYVERSE_MISSIONS_DIR, mission_id_foc,
    "/images/",
    mission_id_foc, "_images.zip"
  )

  download.file(
    url = cyverse_url,
    destfile = zip_filepath,
    method = "curl",
    extra = "-L" # Follow redirects
  )

  # # Option 2: Download with iRODS (2-3x faster, but still just a tiny fraction of the pipeline time,
  # and more robust because it's not restricted to 10 simultaneous connections). But it will fail if
  # there are more than 10 simultaneous requests to irods, which is the case e.g. at the very start of
  # running the parallelized pipeline. So we are not using it until that can be fixed.

  # irods_path = paste0(
  #   CYVERSE_MISSIONS_DIR, mission_id_foc,
  #   "/images/",
  #   mission_id_foc, "_images.zip"
  # )

  # command = paste("iget -P -f -K -T", shQuote(irods_path), shQuote(zip_filepath))
  # system(command)

  # Confirm that download worked

  if (file.exists(zip_filepath)) {
    size = file.info(zip_filepath)$size

    if (size < 1000000) {
      warning("Downloaded zip file is too small (is it an error message?) for mission: ", mission_id_foc)
      return(FALSE)
    }
  } else {
    warning("Zip download failed (no server response or no remote file) for mission: ", mission_id_foc)
    return(FALSE)
  }

  # Unzip it

  unzip_dir = file.path(
    PHOTOGRAMMETRY_DIR, INPUT_IMAGES_SUBDIR,
    mission_id_foc
  )

  if (!dir.exists(unzip_dir)) {
    dir.create(unzip_dir, recursive = TRUE)
  }

  command = paste("unzip -o -q", shQuote(zip_filepath), "-d", shQuote(unzip_dir))
  unzip_result = system(command)

  if (unzip_result != 0) {
    warning("Unzip failed for mission: ", mission_id_foc)
    return(FALSE)
  }

  return(TRUE)
}

# mission_id = "000574"
# download_unzip_images(mission_id)
