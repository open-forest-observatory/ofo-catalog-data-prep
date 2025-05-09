# Purpose: Prepare the metashape config and shell files for processing all image sets

library(dplyr)

# Handle difference in how the current directory is set between debugging and command line call
if (file.exists("deploy/drone-imagery-ingestion/imagery_project_name.txt")) {
  IMAGERY_PROJECT_NAME_FILE = "deploy/drone-imagery-ingestion/imagery_project_name.txt"
} else {
  IMAGERY_PROJECT_NAME_FILE = "imagery_project_name.txt"
}
IMAGERY_PROJECT_NAME = readr::read_lines(IMAGERY_PROJECT_NAME_FILE)

DATASET_DIR = "/ofo-share/drone-imagery-organization/2_sorted"
BASE_YAML_FILEPATH = "/ofo-share/repos-derek/ofo-catalog-data-prep/deploy/drone-imagery-ingestion/full-run-configs/base.yml"
DERIVED_YAML_OUTFOLDER = "/ofo-share/repos-derek/ofo-catalog-data-prep/deploy/drone-imagery-ingestion/full-run-configs/02"
AUTOMATE_METASHAPE_PATH = "/ofo-share/repos-derek/automate-metashape"
METASHAPE_OUTPUT_PATH = "/ofo-share/drone-imagery-processed/01/metashape-outputs"
METASHAPE_PROJECT_PATH = "/ofo-share/drone-imagery-processed/01/metashape-projects"
N_SHELL_SPLITS = 8

dataset_dir = file.path(DATASET_DIR, IMAGERY_PROJECT_NAME)
derived_yaml_out_folder = file.path(DERIVED_YAML_OUTFOLDER, IMAGERY_PROJECT_NAME)

devtools::load_all()

# Processing

# Get the list of datasets
datasets = list.dirs(dataset_dir, recursive = FALSE)

scenarios = data.frame()
for (dataset in datasets) {

  subdirs = list.dirs(dataset, recursive = FALSE)
  config_filename = basename(dataset)

  scenario = data.frame(photo_path = I(list(subdirs)),
                   config_filename = config_filename)

  scenarios = bind_rows(scenarios, scenario)
}

scenarios$output_path = METASHAPE_OUTPUT_PATH
scenarios$project_path = METASHAPE_PROJECT_PATH

make_derived_configs(
  base_yaml_filepath = BASE_YAML_FILEPATH,
  scenarios,
  derived_yaml_out_folder = derived_yaml_out_folder,
  automate_metashape_path = AUTOMATE_METASHAPE_PATH,
  n_shell_splits = min(N_SHELL_SPLITS, nrow(scenarios)) # Make sure we don't request more splits than config files
)
