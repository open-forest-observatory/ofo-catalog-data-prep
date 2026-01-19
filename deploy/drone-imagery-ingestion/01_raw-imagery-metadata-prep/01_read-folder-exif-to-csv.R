# Purpose: Read in the manually cleaned drone images in mission folders and extract EXIF
# data. Save to a file for further processing in next script (and others).

library(tidyverse)
library(exifr)
library(furrr)

source("deploy/drone-imagery-ingestion/00_set-constants.R")
source("src/metadata-extraction_imagery_general.R")


get_image_data = function(dataset_folder) {

  # Get the full path to the dataset folder
  base_folder = basename(dataset_folder)

  # Get list of image files
  image_filepaths = list.files(dataset_folder, full.names = TRUE, recursive = TRUE, pattern = "(.jpg$)|(.jpeg$)|(.JPG$)|(.JPEG$)")

  # Get exif data
  exif = read_exif_drop_thumbnails(image_filepaths)

  # Compile relevant per-image info (including potential ways to distinguish two drones) into data
  # frame

  # Check if required columns exist at all (would indicate all images are corrupted)
  if (is.null(exif$DateTimeOriginal) || is.null(exif$Model)) {
    warning("Required EXIF columns missing entirely in ", base_folder, ". Skipping folder.")
    return(NULL)
  }

  # Filter out individual corrupted images (missing date or model)
  corrupted_mask = is.na(exif$DateTimeOriginal) | is.na(exif$Model)
  n_corrupted = sum(corrupted_mask)
  if (n_corrupted > 0) {
    warning(n_corrupted, " corrupted image(s) skipped in ", base_folder)
    exif = exif[!corrupted_mask, ]
    image_filepaths = image_filepaths[!corrupted_mask]
  }

  # If all images were corrupted, skip the folder
  if (nrow(exif) == 0) {
    warning("All images corrupted in ", base_folder, ". Skipping folder.")
    return(NULL)
  }

  serialnumber_null = is.null(exif$SerialNumber)

  # If the serial number is missing (as in the case of the Matrice 100) replace it with the model
  if (serialnumber_null) {
    exif$SerialNumber = exif$Model
  }
  exif = exif %>%
    mutate(SerialNumber = ifelse(is.na(SerialNumber), Model, SerialNumber))

  exif$folder_in = base_folder
  exif$image_path_in = image_filepaths

  gc()

  return(exif)
}

imagery_input_path = file.path(CONTRIBUTED_IMAGERY_PATH, PROJECT_NAME_TO_PROCESS_RAW_IMAGERY_METADATA)

# All folders
folders = list.dirs(imagery_input_path, full.names = TRUE, recursive = FALSE)

plan = future::plan(multicore)
exif_list = future_map(folders, get_image_data, .progress = TRUE, .options = furrr_options(seed = TRUE, 
                                                                                   scheduling = Inf))

exif = bind_rows(exif_list)

exif_output_path = file.path(RAW_EXIF_PATH, paste0(PROJECT_NAME_TO_PROCESS_RAW_IMAGERY_METADATA, ".csv"))
if(!dir.exists(dirname(exif_output_path))) {
  dir.create(dirname(exif_output_path), recursive = TRUE)
}
write_csv(exif, file.path(exif_output_path))
