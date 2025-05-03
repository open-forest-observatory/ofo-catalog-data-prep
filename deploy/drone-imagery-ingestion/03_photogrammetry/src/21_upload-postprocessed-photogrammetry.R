# Purpose: For the provided mission_id and run_id, upload the post-processed files to CyVerse

# IMPORTANT NOTE: You must have already authenticated irods with your CyVerse account on this
# machine using iinit. Note that by default, OFO dev instances come pre-authenticated for an
# anonymous (read-only) user, so you will need to run the following lines (change the username to
# yours) to authenticate as yourself and then type in your password when prompted
## echo '{"irods_host": "data.cyverse.org", "irods_port": 1247, "irods_user_name": "djyoung", "irods_zone_name": "iplant"}' > /home/exouser/.irods/irods_environment.json; iinit


mission_id_foc = "000545"
run_id_foc = "20250502T0351"

source("deploy/drone-imagery-ingestion/00_set-constants.R")

## Workflow

upload_postprocessed_photogrammetry_to_cyverse = function(mission_id_foc, run_id_foc) {

  cat(
    "\n **** Uploading post-processed photogrammetry outputs to CyVerse for mission",
    mission_id_foc,
    "run",
    run_id_foc,
    "**** \n"
  )

  local_dir = file.path(
    PHOTOGRAMMETRY_DIR,
    PHOTOGRAMMETRY_POSTPROCESSED_SUBDIR,
    mission_id_foc,
    paste0("processed_", run_id_foc)
  )

  remote_dir = file.path(
    CYVERSE_MISSIONS_DIR,
    mission_id_foc
  )

  # Make sure the remote dir ends with a trailing slash. It's OK if it ends in 2, so we don't need to
  # check if it already does and can just add one regardless.
  remote_dir = paste0(remote_dir, "/")

  command = paste("iput -P -f -T -K -r", local_dir, remote_dir)

  result = system(command)
  # With intern = FALSE (default), the result is the return code (0 for success, other val for failure)
    

  # Check if we got an error, and if so, retry up to 10 times, waiting 1 minute between each
  tries = 0
  while (result != 0 && tries < 10) {
    tries = tries + 1
    Sys.sleep(60) # Wait 1 minute before retrying
    cat("\n **** Uploading postprocessed photogrammetry to CyVerse for mission", mission_id_foc, "failed. Retrying... (attempt", tries, ") **** \n")
    result = system(command)
  }

  # If we still got an upload error, print a warning and save to log and return false (don't do the
  # next step which is deleting the local directory)
  if (result != 0) {
    toprint = (paste(Sys.time(), "- Error uploading postprocessed photogrammetry to CyVerse for mission", mission_id_foc, "(all 10 tries failed).\n"))
    warning(toprint)
    write(toprint, file = UPLOAD_ERROR_LOG, append = TRUE)
    return(FALSE)
  }

  # TODO: consider if there is a better way to check if the upload was successful, other than the
  # current approach of simply checking the return code for no errors.

  # If upload was successful (no errors), delete the local directory. TODO: Consider whether this should be a
  # separate function. If it's moved, only call it if the upload was successful.

  unlink(local_dir, recursive = TRUE)
  # If the containing folder (mission_id_foc) is now empty, delete that too
  mission_dir = file.path(PHOTOGRAMMETRY_DIR, PHOTOGRAMMETRY_POSTPROCESSED_SUBDIR, mission_id_foc)
  files_in_mission_dir = list.files(mission_dir) 
  if (length(files_in_mission_dir) == 0) {
    unlink(mission_dir, recursive = TRUE)
  }

  return(TRUE)
}
