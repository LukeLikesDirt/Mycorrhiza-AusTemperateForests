# Load libraries
require(ggtext)
require(terra)
require(data.table)
require(tidyverse)

# (1) Read in the data ---------------------------------------------------------

# Temperate forest vector
temperate_forest_vect <- vect("data/temperate_forest_ecoregions.shp")

# Temperate forest raster
temperate_forest_rast <- rast("data/temperate_forest_raster.tif")

# Sample metadata
metadata <- fread("data/sample_metadata.csv")

# Extract coordinates
my_sites <- metadata %>%
  select(sample_id, longitude, latitude) %>%
  distinct() %>%
  glimpse(.)

# Create base coordinate vector in longlat
coords_longlat <- my_sites %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(., crs = '+proj=longlat')

# Project to forest raster CRS (for forest extraction)
coords_forest <- project(coords_longlat, temperate_forest_rast)

# (2) Read in the covariate rasters --------------------------------------------

# Define bioclim variables
bioclim_vars <- c("bio1.tif", "bio4.tif", "bio7.tif", "bio12.tif", "bio15.tif")

# Create bioclim raster stack
bioclim_stack <- rast(paste0("data/covariates_georeferenced/bioclim/", bioclim_vars))

# Aridity index
aridity_rast <- rast("data/covariates_georeferenced/aridity_index/ADI.flt") %>%
  tidyterra::rename(aridity_index = ADI)

# Define soil variables
soil_vars <- c("N_total_5cm.tif", "P_available_5cm.tif", 
               "pH_CaCl2_5cm.tif", "SOC_5cm.tif")
# Create soilgrid raster stack
soilgrid_stack <- rast(paste0("data/covariates_georeferenced/soil_grid/", soil_vars))

# Name the layers
names(soilgrid_stack) <- c("nitrogen_georef", "phosphorus_georef", 
                           "ph_georef", "carbon_georef")

# Combine all layers into one raster stack
myco_tree_rast <- c(
  rast("data/covariates_georeferenced/basal_area_AM_raster.tif") %>%
    tidyterra::select(am_tree_basal_area = predicted, am_tree_basal_area_pctl = predicted_percentile),
  rast("data/covariates_georeferenced/basal_area_EcM_raster.tif") %>%
    tidyterra::select(ecm_tree_basal_area = predicted, ecm_basal_area_percentile_pctl = predicted_percentile),
  rast("data/covariates_georeferenced/richness_AM_raster.tif") %>%
    tidyterra::select(am_tree_richness = predicted, am_tree_richness_pctl = predicted_percentile),
  rast("data/covariates_georeferenced/richness_EcM_raster.tif") %>%
    tidyterra::select(ecm_tree_richness = predicted, ecm_richness_percentile_pctl = predicted_percentile)
)

# Check the names
names(myco_tree_rast)

# (3) Extract bioclim and soilgrid values --------------------------------------

# Project coordinates to bioclim CRS
coords_bioclim <- project(coords_longlat, bioclim_stack)

# Project coordinates to aridity index CRS
coords_aridity <- project(coords_longlat, aridity_rast)

# Project coordinates to soilgrid CRS
coords_soil <- project(coords_longlat, soilgrid_stack)

# Project coordinates to myco_tree_rast CRS
coords_myco <- project(coords_longlat, myco_tree_rast)

# Extract bioclim values at site coordinates
bioclim_values <- terra::extract(bioclim_stack, coords_bioclim) %>%
  as.data.frame() %>%
  select(-ID)

# Extract aridity index values at site coordinates
aridity_values <- terra::extract(aridity_rast, coords_aridity) %>%
  as.data.frame() %>%
  select(-ID)

# Extract soilgrid values at site coordinates
soilgrid_values <- terra::extract(soilgrid_stack, coords_soil) %>%
  as.data.frame() %>%
  select(-ID)

# Extract myco_tree values at site coordinates
myco_tree_values <- terra::extract(myco_tree_rast, coords_myco) %>%
  as.data.frame() %>%
  select(-ID)

# Check for NAs
sum(is.na(bioclim_values))
sum(is.na(aridity_values))
sum(is.na(soilgrid_values))
sum(is.na(myco_tree_values))

# Fill NAs with nearest neighbour values: bioclim
if (any(is.na(bioclim_values))) {
  na_rows <- which(rowSums(is.na(bioclim_values)) > 0)
  
  for (i in na_rows) {
    bioclim_values[i, ] <- terra::extract(
      bioclim_stack, 
      coords_bioclim[i],
      method = "simple",
      cells = TRUE
    ) %>%
      select(-ID, -cell) %>%
      as.numeric()
    
    # If still NA, use nearest neighbor
    if (any(is.na(bioclim_values[i, ]))) {
      bioclim_values[i, ] <- terra::extract(
        bioclim_stack, 
        coords_bioclim[i],
        method = "bilinear"
      ) %>%
        select(-ID) %>%
        as.numeric()
    }
  }
}

# Fill NAs with nearest neighbour values: aridity index
if (any(is.na(aridity_values))) {
  na_rows <- which(rowSums(is.na(aridity_values)) > 0)
  
  for (i in na_rows) {
    aridity_values[i, ] <- terra::extract(
      aridity_rast, 
      coords_aridity[i],
      method = "simple",
      cells = TRUE
    ) %>%
      select(-ID, -cell) %>%
      as.numeric()
    
    # If still NA, use bilinear interpolation
    if (any(is.na(aridity_values[i, ]))) {
      aridity_values[i, ] <- terra::extract(
        aridity_rast, 
        coords_aridity[i],
        method = "bilinear"
      ) %>%
        select(-ID) %>%
        as.numeric()
    }
  }
}

# Fill NAs with nearest neighbor values: soilgrid
if (any(is.na(soilgrid_values))) {
  na_rows <- which(rowSums(is.na(soilgrid_values)) > 0)
  
  for (i in na_rows) {
    soilgrid_values[i, ] <- terra::extract(
      soilgrid_stack, 
      coords_soil[i],
      method = "simple",
      cells = TRUE
    ) %>%
      select(-ID, -cell) %>%
      as.numeric()
    
    # If still NA, use bilinear interpolation
    if (any(is.na(soilgrid_values[i, ]))) {
      soilgrid_values[i, ] <- terra::extract(
        soilgrid_stack, 
        coords_soil[i],
        method = "bilinear"
      ) %>%
        select(-ID) %>%
        as.numeric()
    }
  }
}

# Fill NAs with nearest neighbour values: myco_tree
if (any(is.na(myco_tree_values))) {
  na_rows <- which(rowSums(is.na(myco_tree_values)) > 0)
  
  for (i in na_rows) {
    myco_tree_values[i, ] <- terra::extract(
      myco_tree_rast, 
      coords_myco[i],
      method = "simple",
      cells = TRUE
    ) %>%
      select(-ID, -cell) %>%
      as.numeric()
    
    # If still NA, use bilinear interpolation
    if (any(is.na(myco_tree_values[i, ]))) {
      myco_tree_values[i, ] <- terra::extract(
        myco_tree_rast, 
        coords_myco[i],
        method = "bilinear"
      ) %>%
        select(-ID) %>%
        as.numeric()
    }
  }
}

# Check for remaining NAs
sum(is.na(bioclim_values))
sum(is.na(aridity_values))
sum(is.na(soilgrid_values))
sum(is.na(myco_tree_values))

# Fill remaining NAs with nearest non-NA cell values: soilgrid
if (any(is.na(soilgrid_values))) {
  na_rows <- which(rowSums(is.na(soilgrid_values)) > 0)
  
  for (i in na_rows) {
    # Get coordinates of site with NA
    site_xy <- crds(coords_soil[i])
    
    # Create a point from the coordinates
    site_point <- vect(site_xy, crs = crs(soilgrid_stack))
    
    # Try increasing search distances until we find non-NA values
    found <- FALSE
    for (search_dist in seq(100, 10000, by = 100)) {
      # Buffer the point
      buffered_point <- buffer(site_point, width = search_dist)
      
      # Crop raster to buffer (much smaller area)
      soil_crop <- crop(soilgrid_stack, buffered_point)
      
      # Get cells as dataframe
      soil_cells_buffered <- as.data.frame(soil_crop, xy = TRUE, na.rm = TRUE)
      
      if (nrow(soil_cells_buffered) > 0) {
        # Calculate distance to all cells in buffer
        distances <- sqrt((soil_cells_buffered$x - site_xy[1])^2 + 
                            (soil_cells_buffered$y - site_xy[2])^2)
        
        # Find the closest one
        min_idx <- which.min(distances)
        
        # Extract values
        soilgrid_values[i, ] <- soil_cells_buffered[min_idx, names(soilgrid_values)]
        
        cat("Row", i, "- filled from nearest cell at distance:", 
            round(distances[min_idx], 2), "m\n")
        found <- TRUE
        break
      }
    }
    
    if (!found) {
      warning("Could not find non-NA values within 100km for row ", i)
    }
  }
}

# Fill remaining NAs with nearest non-NA cell values: myco_tree
if (any(is.na(myco_tree_values))) {
  na_rows <- which(rowSums(is.na(myco_tree_values)) > 0)
  
  for (i in na_rows) {
    # Get coordinates of site with NA
    site_xy <- crds(coords_myco[i])
    
    # Create a point from the coordinates
    site_point <- vect(site_xy, crs = crs(myco_tree_rast))
    
    # Try increasing search distances until we find non-NA values
    found <- FALSE
    for (search_dist in seq(100, 10000, by = 100)) {
      # Buffer the point
      buffered_point <- buffer(site_point, width = search_dist)
      
      # Crop raster to buffer (much smaller area)
      myco_crop <- crop(myco_tree_rast, buffered_point)
      
      # Get cells as dataframe
      myco_cells_buffered <- as.data.frame(myco_crop, xy = TRUE, na.rm = TRUE)
      
      if (nrow(myco_cells_buffered) > 0) {
        # Calculate distance to all cells in buffer
        distances <- sqrt((myco_cells_buffered$x - site_xy[1])^2 + 
                            (myco_cells_buffered$y - site_xy[2])^2)
        
        # Find the closest one
        min_idx <- which.min(distances)
        
        # Extract values
        myco_tree_values[i, ] <- myco_cells_buffered[min_idx, names(myco_tree_values)]
        
        cat("Row", i, "- filled from nearest cell at distance:", 
            round(distances[min_idx], 2), "m\n")
        found <- TRUE
        break
      }
    }
    
    if (!found) {
      warning("Could not find non-NA values within 100km for row ", i)
    }
  }
}

# Check for remaining NAs
sum(is.na(soilgrid_values))
sum(is.na(myco_tree_values))

# (4) Save the metadata with extracted covariates ------------------------------

# Combine extracted values with site metadata
covariate_data <- my_sites %>%
  bind_cols(bioclim_values) %>%
  bind_cols(aridity_values) %>%
  bind_cols(soilgrid_values) %>%
  bind_cols(myco_tree_values)

# Check for NAs
sum(is.na(covariate_data))

# Save to txt file
inner_join(
  metadata, 
  covariate_data, 
  by = c("sample_id", "longitude", "latitude")
  ) %>%
  select(-project, -flow_id) %>%
  fwrite(
    .,
    "data/sample_covariates.txt", 
    sep = "\t"
    )
