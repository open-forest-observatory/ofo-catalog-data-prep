# Purpose: For the provided mission_id and run_id, upload the photogrammetry outputs to the object
# store

# IMPORTANT NOTE: You must have already configured an rclone remote on this machine for the object
# store. This is already done automatically on ofo dev images. The config should look like this:
# https://github.com/open-forest-observatory/ofo-ansible/blob/main/roles/ofo/files/rclone.conf In
# addition, you must have your S3 credentials set in the environment variables RCLONE_S3_ACCESS_KEY_ID and
# RCLONE_S3_SECRET_ACCESS_KEY. You can do this by running the following command (change the values to
# yours): 'export RCLONE_S3_ACCESS_KEY_ID=your_access_key_id; export
# RCLONE_S3_SECRET_ACCESS_KEY=your_secret_access_key' ad then reboot.

source("deploy/drone-imagery-ingestion/00_set-constants.R")

mission_id_foc = "000000"
config_id_foc = "01"

## Workflow

upload_photogrammetry_outputs_to_object_store = function(mission_id_foc, config_id_foc) {

  cat(
    "\n **** Uploading photogrammetry outputs to object store for mission",
    mission_id_foc,
    "config",
    config_id_foc,
    "**** \n"
  )

  # Change ownership of all these files to the current user
  command = paste0(
    "sudo chown ", Sys.getenv("USER"), " ",
    file.path(PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_SUBDIR, paste0(config_id_foc, "_", mission_id_foc, "_*"))
  )
  system(command)

  # And change ownership of the containing folder for metashape outputs and projects
  command = paste0(
    "sudo chown ", Sys.getenv("USER"), " ",
    file.path(PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_SUBDIR)
  )
  system(command)

  command = paste0(
    "sudo chown ", Sys.getenv("USER"), " ",
    file.path(PHOTOGRAMMETRY_DIR, METASHAPE_PROJECT_SUBDIR)
  )
  system(command)

  # Construct data transfer command line call
  local_dir = file.path(PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_SUBDIR)
  remote_dir = paste0(RCLONE_REMOTE, ":", REMOTE_PHOTOGRAMMETRY_DIR, paste0("config_", config_id_foc))
  filter_phrase = paste0("--include ", config_id_foc, "_", mission_id_foc, "_*")
  command = paste("rclone copy", local_dir, remote_dir, filter_phrase, "--max-depth 1 --progress --transfers 32 --checkers 32 --stats 1s --retries 5 --retries-sleep=15s --s3-upload-cutoff 100Mi --s3-chunk-size 100Mi --s3-upload-concurrency 16", sep = " ")
  result = system(command)

  # Check if we got an error, and if so, retry up to 10 times, waiting 1 minute between each
  tries = 0
  while (result != 0 && tries < 10) {
    tries = tries + 1
    Sys.sleep(60) # Wait 1 minute before retrying
    cat("\n **** Uploading photogrammetry outputs to object store for mission", mission_id_foc, "failed. Retrying... (attempt", tries, ") **** \n")
    result = system(command)
  }

  # If we still got an upload error, print a warning and save to log and return false (don't do the
  # next step which is deleting the local directory)
  if (result != 0) {
    toprint = (paste(Sys.time(), "- Error uploading photogrammetry outputs to object store for mission", mission_id_foc, "(all 10 tries failed).\n"))
    warning(toprint)
    write(toprint, file = UPLOAD_ERROR_LOG, append = TRUE)
    return(FALSE)
  }

  # TODO: consider if there is a better way to check if the upload was successful, other than the
  # current approach of simply checking the return code for no errors.

  # If upload was successful (no errors), delete the just-uploaded files for the processing run. TODO: Consider whether this should be a
  # separate function. If it's moved, only call it if the upload was successful.

  uploaded_files = list.files(local_dir, pattern = paste0(config_id_foc, "_", mission_id_foc, "_*"), full.names = TRUE)

  unlink(uploaded_files, recursive = TRUE)
  # If the containing folder (photogrammetry-outputs) is now empty, delete that too
  mission_dir = file.path(PHOTOGRAMMETRY_DIR, METASHAPE_OUTPUT_SUBDIR)
  files_in_mission_dir = list.files(mission_dir)
  if (length(files_in_mission_dir) == 0) {
    unlink(mission_dir, recursive = TRUE)
  }

  return(TRUE)
}
