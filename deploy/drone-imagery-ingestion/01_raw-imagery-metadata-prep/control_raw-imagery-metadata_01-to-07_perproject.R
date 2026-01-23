# Purpose: a controller script that runs the raw imagery ingestion pipeline from script 01 to 11.

repo_root = "/ofo-share/repos/derek/ofo-catalog-data-prep"

cat("\n\n****Starting script 1****\n\n")
source(file.path(repo_root,"deploy/drone-imagery-ingestion/01_raw-imagery-metadata-prep/01_read-folder-exif-to-csv.R"))
cat("\n\n****Starting script 2****\n\n")
source(file.path(repo_root,"deploy/drone-imagery-ingestion/01_raw-imagery-metadata-prep/02_plan-imagery-folder-standardization.R"))
cat("\n\n****Starting script 3****\n\n")
source(file.path(repo_root,"deploy/drone-imagery-ingestion/01_raw-imagery-metadata-prep/03_compile-baserow-metadata.R"))
cat("\n\n****Starting script 4****\n\n")
source(file.path(repo_root,"deploy/drone-imagery-ingestion/01_raw-imagery-metadata-prep/04_determine-missions-to-process.R"))
cat("\n\n****Starting script 5****\n\n")
source(file.path(repo_root,"deploy/drone-imagery-ingestion/01_raw-imagery-metadata-prep/05_parse-exif-metadata-per-image.R"))
cat("\n\n****Starting script 6****\n\n")
source(file.path(repo_root,"deploy/drone-imagery-ingestion/01_raw-imagery-metadata-prep/06_summarize-exif-metadata-per-mission.R"))
cat("\n\n****Starting script 7****\n\n")
source(file.path(repo_root,"deploy/drone-imagery-ingestion/01_raw-imagery-metadata-prep/07_merge-exif-and-baserow.R"))
