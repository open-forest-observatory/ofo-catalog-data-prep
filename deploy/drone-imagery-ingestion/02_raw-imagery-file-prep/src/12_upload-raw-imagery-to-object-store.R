## Purpose: Upload the data in the upload staging directory created in previous step to CyVerse.

# IMPORTANT NOTE: You must have already configured an rclone remote on this machine for the object
# store. This is done entirely througn env vars:
# export RCLONE_CONFIG_JS2S3_TYPE=s3
# export RCLONE_CONFIG_JS2S3_PROVIDER=Other
# export RCLONE_CONFIG_JS2S3_ENDPOINT=https://js2.jetstream-cloud.org:8001/
# export RCLONE_CONFIG_JS2S3_ENV_AUTH=true
#  export AWS_ACCESS_KEY_ID=<access-key-id>
#  export AWS_SECRET_ACCESS_KEY=<secret-access-key>

upload_raw_imagery_to_object_store = function(mission_id_foc) {

  cat("\n **** Uploading raw imagery to remote object storage for mission", mission_id_foc, "**** \n")

  # Construct data transfer command line call
  local_dir = file.path(UPLOAD_STAGING_DIR_PATH, mission_id_foc)

  remote_dir = paste0(RCLONE_REMOTE, ":", REMOTE_MISSIONS_DIR, mission_id_foc)

  command = paste("rclone copy", local_dir, remote_dir, "--progress --transfers 32 --checkers 32 --stats 1s --retries 5 --retries-sleep=15s --s3-upload-cutoff 100Mi --s3-chunk-size 100Mi --s3-upload-concurrency 16 --config /dev/null", sep = " ")
  # There are only 11 files so as long as there are at least 11 tranfers we are maxed on throughput

  # # OLD METHOD FOR CYVERSE IRODS

  # # When the remote dir ends in a trailing slash, make sure you do not specify the base folder
  # # name of the local dir being transferred (e.g. `000350`), as it will be created within the remote
  # # dir.
  
  # # Make sure the remote dir ends with a trailing slash. It's OK if it ends in 2, so we don't need to
  # # check if it already does and can just add one regardless.
  # CYVERSE_MISSIONS_DIR = paste0(CYVERSE_MISSIONS_DIR, "/")

  # command = paste("iput -P -f -T -K -r", local_dir, CYVERSE_MISSIONS_DIR)

  result = system(command)
  # With intern = FALSE (default), the result is the return code (0 for success, 2 for failure)
  
  # Check if we got an error, and if so, retry up to 3 times, waiting 1 minute between each
  tries = 0
  while (result != 0 && tries < 3) {
    tries = tries + 1
    Sys.sleep(60) # Wait 1 minute before retrying
    cat("\n **** Uploading raw imagery to remote object storage for mission", mission_id_foc, "failed. Retrying... (attempt", tries, ") **** \n")
    result = system(command)
  }

  # If we still got an upload error, print a warning and save to log and return false
  if (result != 0) {
    toprint = (paste(Sys.time(), "- Error uploading raw imagery to remote object storage for mission", mission_id_foc, "(all 3 tries failed).\n"))
    warning(toprint)
    write(toprint, file = UPLOAD_ERROR_LOG, append = TRUE)
    return(FALSE)
  }

  return(TRUE)
}


# # Function for checking whether the zip file upload exists on cyverse
# check_if_zip_file_exists_cyverse = function(mission_id_foc) {
#   # Check if the zip file exists on CyVerse
#   coll_name = shQuote(paste0(CYVERSE_MISSIONS_DIR, mission_id_foc, "/images"))
#   data_name = shQuote(paste0(mission_id_foc, "_images.zip"))
#   call = paste0("iquest --no-page '%s/%s' \"select COLL_NAME, DATA_NAME where COLL_NAME like ", coll_name, " and DATA_NAME like ", data_name, "\"")
#   # Make any double-slashes into single slashes
#   call = gsub("//", "/", call)
#   uploaded = system(call, intern = TRUE)
#   uploaded_relevant_files = any(str_detect(uploaded, fixed("/iplant/home/shared/ofo/public/mission")))

#   return(uploaded_relevant_files)
# }

# # Function for checking whether the zip file upload exists on osn
# check_if_zip_file_exists = function(mission_id_foc) {
#   # Check if the zip file exists on CyVerse

#   remote_zip_path = paste0(OSN_RCLONE_REMOTE, ":", OSN_MISSIONS_DIR, mission_id_foc, "/images/", mission_id_foc, "_images.zip")

#   call = paste0("rclone lsf ", remote_zip_path)

#   uploaded = system(call, intern = TRUE)
#   uploaded_success = length(uploaded) > 0

#   return(uploaded_success)
# }
