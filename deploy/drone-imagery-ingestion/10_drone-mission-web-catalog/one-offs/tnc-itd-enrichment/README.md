The scripts in this folder are for adding additional detail to the photogrammetry and ITD products to support deliverables due to TNC. Tasks that these scripts perform:

- Composite missions
  - Download all composite mission polygons, take the intersection of the 2 layers of each (to convert the 2 layers into 1), and bind into one GPKG
    - [polygons are on S3 in the ofo-public bucket with the following file path format: drone/mission-composites_01/000001_000002/metadata-mission/000001_000002_mission-metadata.gpkg]
    - [for guidance on intersecting, see what is done for displaying the mission polygon in the composite mission map in src/web-catalog-creation_drone-imagery-catalog.R]
    - [save to /ofo-share/project-data/tnc-yuba-deliverables/composite-missions/overall/composite-drone-plot-summaries.gpkg]
- Nadir missions
  - Download all nadir mission polygons and bind into one GPKG (after filtering as described below)
    - [save to (deliverables-dir)/individual-missions/overall/individual-drone-plot-summaries.gpkg]
    - [polygons are on S3 in the ofo-public bucket with the following file path format: drone/missions_03/000001/metadata-mission/000001_mission-metadata.gpkg]
    - These include both nadir and oblique missions. Need to filter to nadir only based on the attribute camera_pitch_nominal or camera_pitch_derived (at least one must be < 10) and altitude_agl_nominal (must be > 90)
- CHMs for nadir missions
  - Download all CHM data from *nadir* missions only (using the nadir mission drone plot summaries file to know what's nadir, and using the S3 file listing to know what's available)
    - [save individual CHMs to /ofo-share/project-data/tnc-yuba-deliverables/individual-missions/canopy-height-models/]
    - [CHMs are on S3 in the ofo-public bucket with the following file structure: drone/missions_03/000001/photogrammetry_03/full/000001_chm-mesh.tif (all CHMs are in missions_03)]
  - Compute overall canopy cover for each and add as an attribute to the nadir mission gpkg. Canopy is defined as anything over 5 m height. Also compute average height of everything over 5 m height.
    - [same drone-plot-summaries.gpkg file]
    - [attributes named 'canopy_cover' and 'canopy_height']
  - Compute a raster for each nadir mission, based on the CHM, that depicts canopy cover in a 100 m square (moving window pixel by pixel)
    - [save to /ofo-share/project-data/tnc-yuba-deliverables/individual-missions/canopy-cover-rasters/]
- ITD data
  - Download all *composite* ITD data and compile it into a single gpkg
    - [save to (deliverables-dir)/composite-missions/all-detected-trees.gpkg]
  - Use the cloud2trees package to estimate the DBH of each tree based on its height and save as an attribute to all-detected-trees.gpkg
  - For each composite drone plot, compute the following summary statistics and save as attributes to the drone-plot-summaries.gpkg file: tree density (trees per ha and acre), basal area (sq m per ha, sq ft per acre), density of trees with DBH > 30 cm (ha and acre), and proportion pine (species code starting with 'PI') (by basal area). Also make a raster spanning the bounds of the plot that is a "heatmap" of proportion pine (at the 1 ha scale, moving window, 10 m resolution) and another that is for density of trees with DBH > 30 cm. These rasters should go in (deliverables-dir)/composite-missions/heatmaps-pine-proportion/ and (deliverables-dir)/composite-missions/heatmaps-large-trees/


  General guidance: Include an initial script that pulls and saves a S3 file listing to a temp dir in /ofo-share/project-data/tnc-yuba-deliverables/tmp/s3-file-listing.csv following the approach in deploy/drone-imagery-ingestion/10_drone-mission-web-catalog/05b_get-all-composite-metadata.R. Other scripts that need to find and download S3 data should reference this. Download the data (in the other scripts) using HTTP so no auth is required, following the approach used for accessing the ITD treetop points in src/web-catalog-creation_drone-imagery-catalog.R.