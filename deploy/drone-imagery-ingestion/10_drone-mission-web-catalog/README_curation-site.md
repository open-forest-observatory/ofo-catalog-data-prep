# Drone Data Curation Site Generation

This directory contains scripts for generating both the public drone data catalog and an internal curation site for data analysts.

## Scripts Overview

- **10_create-drone-data-catalog-webpages.R**: Generates the public-facing drone data catalog with full mission details, processed data products, and download links
- **11_create-drone-data-curation-webpages.R**: Generates the internal curation site with simplified pages showing only mission metadata and image location maps

## Curation Site Features

The curation site (`11_create-drone-data-curation-webpages.R`) is designed for internal data quality review and provides:

1. **Mission-level metadata table**: All metadata attributes for each mission
2. **Interactive leaflet map**: Shows drone photo locations, flight path, and image attributes
3. **Secondary metadata support**: Allows using alternative image metadata for specific missions
4. **Mission override list**: Specify which missions should use secondary metadata

## Using the Curation Site Generator

### Basic Usage

To generate the curation site with default settings:

```r
source("deploy/drone-imagery-ingestion/10_drone-mission-web-catalog/11_create-drone-data-curation-webpages.R")
```

This will:
- Use mission metadata from: `/ofo-share/project-data/catalog-data-prep/05_drone-imagery-web-catalog/01_metadata/mission-metadata.gpkg`
- Use primary image metadata from: `/ofo-share/project-data/catalog-data-prep/05_drone-imagery-web-catalog/01_metadata/image-metadata.gpkg`
- Generate pages in the website repo at: `ofo-website-3/content/data/drone-curation/mission-details/`

### Using Secondary Image Metadata

To use alternative image metadata for specific missions:

1. **Prepare your secondary metadata file** (`.gpkg` format with same structure as primary)

2. **Update the override list CSV** with missions that should use secondary metadata (single column named `mission_id`)

3. **Edit the constants file** (`deploy/drone-imagery-ingestion/00_set-constants.R`) to specify both paths:
   ```r
   # Set both to actual file paths:
   SECONDARY_IMAGE_METADATA_FILEPATH = "/path/to/your/secondary-image-metadata.gpkg"
   IMAGERY_METADATA_MISSION_OVERRIDE_LIST_FILEPATH = "/path/to/your/override-list.csv"

   # Or set both to empty strings ("") to disable override functionality:
   SECONDARY_IMAGE_METADATA_FILEPATH = ""
   IMAGERY_METADATA_MISSION_OVERRIDE_LIST_FILEPATH = ""
   ```

4. **Run the script**:
   ```r
   source("deploy/drone-imagery-ingestion/10_drone-mission-web-catalog/11_create-drone-data-curation-webpages.R")
   ```

**Important:** Both `SECONDARY_IMAGE_METADATA_FILEPATH` and `IMAGERY_METADATA_MISSION_OVERRIDE_LIST_FILEPATH` must be either both populated (with valid paths) or both empty strings. Don't mix populated and empty.

The script will automatically:
- Use secondary metadata for missions listed in the override list
- Use primary metadata for all other missions
- Report how many missions use each metadata source

### Example Override List

The override list CSV should look like:

```csv
mission_id
000013
000019
000020
000022
```

## Configuration

Key constants are defined in `deploy/drone-imagery-ingestion/00_set-constants.R`:

- `MISSION_METADATA_FILEPATH`: Path to primary mission metadata (.gpkg)
- `IMAGE_METADATA_FILEPATH`: Path to primary image metadata (.gpkg)
- `SECONDARY_IMAGE_METADATA_FILEPATH`: Path to secondary image metadata (.gpkg), or empty string
- `IMAGERY_METADATA_MISSION_OVERRIDE_LIST_FILEPATH`: Path to override list CSV, or empty string
- `CURATION_MISSION_DETAILS_TEMPLATE_FILEPATH`: Template for curation pages
- `CURATION_MISSION_DETAILS_PAGE_DIR`: Output directory for curation pages
- `CURATION_MISSION_DETAILS_DATATABLE_DIR`: Output directory for datatables
- `CURATION_MISSION_DETAILS_MAP_DIR`: Output directory for maps

## Output Structure

The script generates:

1. **Mission detail pages**: One markdown file per mission in `data/drone-curation/mission-details/`
2. **Leaflet maps**: Interactive HTML maps in `static/drone-curation-mission-details-maps/`
3. **Data tables**: HTML tables in `static/drone-curation-mission-details-datatables/`
4. **Header files**: Shared JavaScript/CSS libraries for datatables and leaflet

## Differences from Public Catalog

The curation site differs from the public catalog in these ways:

| Feature | Public Catalog | Curation Site |
|---------|---------------|---------------|
| Processed products | ✓ (orthomosaics, CHMs, etc.) | ✗ |
| Example images | ✓ | ✗ |
| Download links | ✓ | ✗ |
| Mission metadata | ✓ | ✓ |
| Image location map | ✓ | ✓ |
| Secondary metadata | ✗ | ✓ |
| Override list | ✗ | ✓ |

## Code Reuse

Both scripts reuse the same core functions from:
- `src/web-catalog-creation_shared-functions.R`: Common widget saving functions
- `src/web-catalog-creation_drone-imagery-catalog.R`: Mission-specific catalog functions

This ensures consistency between public and internal sites while minimizing code duplication.

## Troubleshooting

**Issue**: Script can't find metadata files
- **Solution**: Verify paths in constants file and ensure `.gpkg` files exist

**Issue**: Override list not working
- **Solution**: Verify CSV has `mission_id` column header and mission IDs match exactly
- **Solution**: Ensure both `SECONDARY_IMAGE_METADATA_FILEPATH` and `IMAGERY_METADATA_MISSION_OVERRIDE_LIST_FILEPATH` are set to valid file paths (not empty strings)

**Issue**: Secondary metadata file not being used
- **Solution**: Ensure both `SECONDARY_IMAGE_METADATA_FILEPATH` and `IMAGERY_METADATA_MISSION_OVERRIDE_LIST_FILEPATH` are populated with valid file paths
- **Solution**: Both must be either populated or both empty strings - don't mix

**Issue**: Website repo not found
- **Solution**: Verify `WEBSITE_REPO_PATH` in `00_set-constants.R` points to the `ofo-website-3` submodule
