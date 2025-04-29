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
  # check if it already does and can just add one regardless.
  CYVERSE_MISSIONS_DIR = paste0(CYVERSE_MISSIONS_DIR, "/")

  command = paste("iput -P -f -T -K -r", local_dir, CYVERSE_MISSIONS_DIR)

  result = system(command)
  # With intern = FALSE (default), the result is the return code (0 for success, 2 for failure)
  
  # Check if we got an error, and if so, retry up to 3 times
  tries = 0
  while (result == 2 && tries < 3) {
    tries = tries + 1
    cat("\n **** Uploading raw imagery to CyVerse for mission", mission_id_foc, "failed. Retrying... (attempt", tries, ") **** \n")
    result = system(command)
  }

  # If we still got an upload error, print a warning and save to log and return false
  if (result == 2) {
    toprint = (paste(Sys.time(), "- Error uploading raw imagery to CyVerse for mission", mission_id_foc, "(all 3 tries failed).\n"))
    warning(toprint)
    write(toprint, file = UPLOAD_ERROR_LOG, append = TRUE)
    return(FALSE)
  }

  return(TRUE)
}


# Function for checking whether the zip file upload exists on cyverse
check_if_zip_file_exists = function(mission_id_foc) {
  # Check if the zip file exists on CyVerse
  coll_name = shQuote(paste0(CYVERSE_MISSIONS_DIR, mission_id_foc, "/images"))
  data_name = shQuote(paste0(mission_id_foc, "_images.zip"))
  call = paste0("iquest --no-page '%s/%s' \"select COLL_NAME, DATA_NAME where COLL_NAME like ", coll_name, " and DATA_NAME like ", data_name, "\"")
  # Make any double-slashes into single slashes
  call = gsub("//", "/", call)
  uploaded = system(call, intern = TRUE)
  uploaded_relevant_files = any(str_detect(uploaded, fixed("/iplant/home/shared/ofo/public/mission")))

  return(uploaded_relevant_files)
}
