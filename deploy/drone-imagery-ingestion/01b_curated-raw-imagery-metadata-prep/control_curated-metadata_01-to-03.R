# deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/control_curated-metadata_01-to-03.R
# Purpose: Run the post-curation metadata preparation pipeline
#
# This control script runs the three post-curation scripts in sequence:
# 1. Apply curation filters (remove extraneous images, update IDs for formerly curated missions)
# 2. Re-summarize metadata at mission/sub-mission level with curated images
# 3. Merge metadata with contributed metadata and add curation anomaly notes

repo_root = "/ofo-share/repos/derek/ofo-catalog-data-prep"

cat("\n\n**** Starting post-curation script 01: Apply curation filters ****\n\n")
source(file.path(repo_root, "deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/01_apply-curation-filters.R"))

cat("\n\n**** Starting post-curation script 02: Summarize curated metadata ****\n\n")
source(file.path(repo_root, "deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/02_summarize-curated-metadata-per-mission.R"))

cat("\n\n**** Starting post-curation script 03: Merge metadata and curation notes ****\n\n")
source(file.path(repo_root, "deploy/drone-imagery-ingestion/01b_curated-raw-imagery-metadata-prep/03_merge-metadata-and-curation-notes.R"))

cat("\n\n**** Post-curation metadata preparation complete! ****\n\n")
