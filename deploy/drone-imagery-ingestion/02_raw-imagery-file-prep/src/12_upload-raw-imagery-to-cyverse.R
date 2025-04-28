## Purpose: Upload the data in the upload staging directory created in previous step to CyVerse.

# IMPORTANT NOTE: You must have already authenticated irods with your CyVerse account on this
# machine using iinit

upload_raw_imagery_to_cyverse = function(mission_id_foc) {

  # Construct data transfer command line call
  local_dir = file.path(UPLOAD_STAGING_DIR_PATH, mission_id_foc)

  # When the remote dir ends in a trailing slash, make sure you do not specify the base folder
  # name of the local dir being transferred (e.g. `000350`), as it will be created within the remote
  # dir.
  
  # Make sure the remote dir ends with a trailing slash. It's OK if it ends in 2, so we don't need to
  # check if it already does
  CYVERSE_MISSIONS_DIR = paste0(CYVERSE_MISSIONS_DIR, "/")

  command = paste("iput -P -f -T -K -r", local_dir, CYVERSE_MISSIONS_DIR)

  system(command)

  return(TRUE)
}
