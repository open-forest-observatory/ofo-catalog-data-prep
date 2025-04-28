## Purpose: Upload the data in the upload staging directory created in previous step to CyVerse.

# IMPORTANT NOTE: You must have already authenticated irods with your CyVerse account on this
# machine using iinit

upload_raw_imagery_to_cyverse = function(mission_id_foc) {

  cat("\n **** Uploading raw imagery to CyVerse for mission", mission_id_foc, "**** \n")

  # Construct data transfer command line call
  local_dir = file.path(UPLOAD_STAGING_DIR_PATH, mission_id_foc)

  # When the remote dir ends in a trailing slash, make sure you do not specify the base folder
  # name of the local dir being transferred (e.g. `000350`), as it will be created within the remote
  # dir.
  
  # Make sure the remote dir ends with a trailing slash. It's OK if it ends in 2, so we don't need to
  # check if it already does
  CYVERSE_MISSIONS_DIR = paste0(CYVERSE_MISSIONS_DIR, "/")

  command = paste("iput -P -f -T -K -r", local_dir, CYVERSE_MISSIONS_DIR)

  # result_file = "/ofo-share/drone-imagery-organization/temp/cyverse-upload-log.txt"
  # to_write = paste0("\n **** Uploading raw imagery to CyVerse for mission", mission_id_foc, "**** \n",
  #     "Command: ", command, "\n")
  # write(to_write, file = result_file, append = TRUE)

  result = system(command, intern = TRUE)
  result = paste(result, collapse = "\n")

  # Check if we got an error
  if (any(grepl("ERROR", result))) {
    toprint = (paste("*#*#*#*#*# Error uploading raw imagery to CyVerse for mission", mission_id_foc, ". Result was:\n", result))
    cat(toprint)
  } else {
    toprint = paste0("\n **** Successfully uploaded raw imagery to CyVerse for mission", mission_id_foc,
        ". Result was :\n", result, "\n")
    cat(toprint)
  }

  # Append result to a file
  result_file = "/ofo-share/drone-imagery-organization/temp/cyverse-upload-log.txt"
  write(toprint, file = result_file, append = TRUE)

  return(TRUE)
}
