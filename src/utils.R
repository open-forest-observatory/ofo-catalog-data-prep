# Create a directory if it doesn't exist
create_dir <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

drop_units_if_present = function(x) {
  if (inherits(x, "units")) {
    return(x |> units::drop_units())
  } else {
    return(x)
  }
}

### **** UTILS from below this line were copied from the ofo-r repo ****

# Reproject a sf object into the CRS representing its local UTM zone
transform_to_local_utm = function(sf) {
  geo = sf::st_transform(sf, 4326)
  geo_noz = sf::st_zm(geo, drop = TRUE)
  lonlat = sf::st_centroid(geo_noz) |> sf::st_coordinates()
  utm = lonlat_to_utm_epsg(lonlat)

  sf_transf = sf::st_transform(sf, utm)

  return(sf_transf)
}

# Take a lon/lat coordinates dataframe and convert to the local UTM zone EPSG code
lonlat_to_utm_epsg = function(lonlat) {
  utm = (floor((lonlat[, 1] + 180) / 6) %% 60) + 1
  utms = ifelse(lonlat[, 2] > 0, utm + 32600, utm + 32700)

  utms_unique = unique(utms)

  if (length(utms_unique) > 2) {
    stop("The geometry spans 3 or more UTM zones")
  } else if (length(utms_unique) > 1) {
    if (abs(diff(utms_unique)) > 1) {
      stop("The geometry spans 2 non-adjacent UTM zones.")
    }
  }

  return(utms_unique[1])
}

# Strip double slashes from paths (except preserve the double slash after "https:")
strip_double_slashes = function(path) {
  path = gsub("([^:])//+", "\\1/", path)
  return(path)
}
