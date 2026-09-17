# Load libraries
require(ggtext)
require(terra)
require(data.table)
require(tidyverse)

# (1) Trim forest raster -------------------------------------------------------

# Forest raster
forest_rast <- rast(
  "data/georeferenced_covariates/native_forest.tif"
)

# Temperate forest ecoregions
temperate_forest_vect <- vect("data/georeferenced_covariates/Ecoregions2017/Ecoregions2017.shp") %>%
  tidyterra::filter(BIOME_NAME == "Temperate Broadleaf & Mixed Forests") %>%
  # Crop to Australia extent
  terra::crop(ext(110, 155, -45, -10))

# Ensure CRS match
temperate_forest_vect <- project(temperate_forest_vect, forest_rast)

# Crop and mask the raster to the ecoregion
forest_rast_trimmed <- crop(forest_rast, temperate_forest_vect, mask = TRUE)

# Check the result
plot(forest_rast_trimmed)

# (2) Trim tree inventory plots ------------------------------------------------

# Mycorrhizal types
mycorrhizal_types <- fread("data/tree_mycorrhizal_types.txt") %>%
  rename(species = scientific_name)

# Load tree inventory plots
trees <- fread("data/tree_induvidual_measurments_dbh.txt") %>%
  rename(species = scientific_name) %>%
  # Join with mycorrhizal types
  left_join(
    mycorrhizal_types %>%
      select(
        genus, mycorrhizal_type
      ),
    by = c("genus"),
    relationship = "many-to-many"
  ) %>%
  # Keep only certain mycorrhizal types
  filter(
    mycorrhizal_type %in% c("AM", "NM-AM", "EcM", "EcM-AM", "ErM", "NM")
  ) %>%
  # Create a new AM column where "AM" and "NM-AM" have a value of 1, and 
  # "EcM-AM" have a values of 0.5, and all others have a value of 0
  mutate(
    AM = case_when(
      mycorrhizal_type %in% c("AM", "NM-AM") ~ 1,
      mycorrhizal_type == "EcM-AM" ~ 0.5,
      TRUE ~ 0
    )
  ) %>%
  # Create a new EcM column where "EcM" have a value of 1, and "EcM-AM" have a
  # values of 0.5, and all others have a value of 0
  mutate(
    EcM = case_when(
      mycorrhizal_type == "EcM" ~ 1,
      mycorrhizal_type == "EcM-AM" ~ 0.5,
      TRUE ~ 0
    )
  ) %>%
  glimpse(.)

# Extract unique sites with coordinates
all_sites <- trees %>%
  select(site, longitude, latitude) %>%
  distinct() %>%
  glimpse(.)

# Project longlat according to the forest_rast
coords <- all_sites %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  # Longitude and latitude as a spatial vector
  vect(., crs = '+proj=longlat') %>%
  # Project coordinates according to forest_rast
  project(., forest_rast_trimmed)

# Forest structure values for each plot
structure_values <- terra::extract(forest_rast_trimmed, coords)

# Forest sites
forest_sites <- all_sites %>%
  bind_cols(structure_values) %>%
  filter(
    !is.na(FOR_CATEGO)
  ) %>%
  select(-FOR_CATEGO)

# Forest trees
forest_trees <- trees %>%
  filter(
    site %in% forest_sites$site
  ) %>%
  distinct() %>%
  glimpse(.)

# Number of forest tree species
forest_trees %>%
  select(species) %>%
  distinct() %>%
  nrow() %>%
  message("Number of tree species: ", .)

# Number of forest sites
forest_sites %>%
  select(site) %>%
  distinct() %>%
  nrow() %>%
  message("Number of sites: ", .)

# (3) Agglomerate forest sites to 1 km grids ----------------------------------

# Create a spatial vector from forest sites
forest_sites_vect <- forest_sites %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  vect(crs = '+proj=longlat') %>%
  project(forest_rast_trimmed)

# Extract cell numbers for each site
forest_sites <- forest_sites %>%
  mutate(cell = cells(forest_rast_trimmed, forest_sites_vect)[, "cell"])

# Get unique cells
unique_cells <- unique(forest_sites$cell)

# Get cell center coordinates
cell_centers <- xyFromCell(forest_rast_trimmed, unique_cells) %>%
  as.data.frame()

# Convert cell centers back to longitude/latitude
cell_centers_longlat <- cell_centers %>%
  as.matrix() %>%
  vect(crs = crs(forest_rast_trimmed)) %>%
  project('+proj=longlat') %>%
  crds() %>%
  as.data.frame() %>%
  rename(longitude = x, latitude = y) %>%
  mutate(cell = unique_cells)

# Aggregate sites by grid cell with cell center coordinates
sites_per_cell <- forest_sites %>%
  group_by(cell) %>%
  summarise(
    n_sites = n(),
    sites = list(site)
  ) %>%
  ungroup() %>%
  left_join(cell_centers_longlat, by = "cell")

# Aggregate tree data by grid cell WITH CELL CENTER COORDINATES
forest_trees_by_cell <- forest_trees %>%
  left_join(
    forest_sites %>% select(site, cell),
    by = "site"
  ) %>%
  select(-longitude, -latitude) %>%  # Remove original coordinates
  left_join(cell_centers_longlat, by = "cell")  # Add cell center coordinates

# Compute mycorrhizal proportions per grid cell
mycorrhizal_proportions_by_cell <- forest_trees_by_cell %>%
  group_by(cell) %>%
  summarise(
    richness = n_distinct(species),
    total_abundance = n(),
    total_diameter = sum(diameter_cm, na.rm = TRUE),
    AM_abundance = sum(AM, na.rm = TRUE),
    EcM_abundance = sum(EcM, na.rm = TRUE),
    AM_abundance_prop = AM_abundance / total_abundance,
    EcM_abundance_prop = EcM_abundance / total_abundance,
    AM_diameter = sum(diameter_cm * AM, na.rm = TRUE),
    EcM_diameter = sum(diameter_cm * EcM, na.rm = TRUE),
    AM_diameter_prop = AM_diameter / total_diameter,
    EcM_diameter_prop = EcM_diameter / total_diameter
  ) %>%
  ungroup()

# Check results
message("Number of unique grid cells: ", nrow(mycorrhizal_proportions_by_cell))
glimpse(mycorrhizal_proportions_by_cell)

# Check the ralationship between abundance proportions and diameter proportions
mycorrhizal_proportions_by_cell %>%
  select(ends_with("prop")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    mycorrhizal_type = if_else(
      str_detect(variable, "AM"), "Arbuscular mycorrhizal", "Ectomycorrhizal"
    ),
    metric = if_else(
      str_detect(variable, "abundance"), "Abundance", "Basal area"
    ),
    value = case_when(
      metric == "Abundance" & mycorrhizal_type == "Arbuscular mycorrhizal" ~ value^1.5,
      metric == "Abundance" & mycorrhizal_type == "Ectomycorrhizal" ~ sqrt(value),
      metric == "Basal area" ~ value
    )
  ) %>%
  select(-variable) %>%  # Remove the original variable column
  pivot_wider(
    names_from = metric,
    values_from = value
  ) %>%
  unnest(cols = c(Abundance, `Basal area`)) %>%  # Unnest the list columns
  ggplot(aes(x = Abundance, y = `Basal area`)) +
  geom_point(color = "black", alpha = 0.5) +
  # Add a 1:1 line
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggpubr::stat_cor(
    aes(label = after_stat(rr.label)),
  ) +
  facet_wrap(~mycorrhizal_type) +
  theme_minimal() +
  theme(
    plot.tag = element_markdown(size = 14),
    axis.title = element_markdown(size = 10),
    axis.text = element_markdown(size = 9),
    aspect.ratio = 1
  )

# Check how many site of abundance only data we can potentially use to increase 
# sample size
abundance_sites <- fread("data/tree_site_measurement_abundance.txt") %>%
  select(site, longitude, latitude) %>%
  unique()

# Trim abundance sites to temperate forests
abundance_sites_coords <- abundance_sites %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  # Longitude and latitude as a spatial vector
  vect(., crs = '+proj=longlat') %>%
  # Project coordinates according to forest_rast
  project(., forest_rast_trimmed)

# Extract forest structure values for abundance sites
abundance_structure_values <- terra::extract(forest_rast_trimmed, abundance_sites_coords)

# Filter to forest sites only
abundance_forest_sites <- abundance_sites %>%
  bind_cols(abundance_structure_values) %>%
  filter(
    !is.na(FOR_CATEGO)
  ) %>%
  select(-FOR_CATEGO)

# Check how many sites are in temperate forests
message("Total abundance sites: ", nrow(abundance_sites))
message("Abundance sites in temperate forests: ", nrow(abundance_forest_sites))

# Check if these are additional sites (not already in forest_sites)
additional_sites <- abundance_forest_sites %>%
  filter(!site %in% forest_sites$site)

message("Additional sites from abundance data: ", nrow(additional_sites))

# (4) Read in the covariate rasters --------------------------------------------

# Define bioclim variables
bioclim_vars <- c("bio1.tif", "bio4.tif", "bio5.tif", "bio6.tif", 
                  "bio12.tif", "bio15.tif", "bio16.tif", "bio17.tif")

# Create bioclim raster stack
bioclim_stack <- rast(paste0("data/georeferenced_covariates/bioclim/", bioclim_vars))

# Define soil variables
soil_vars <- c("N_total_5cm.tif", "P_available_5cm.tif", 
               "pH_CaCl2_5cm.tif", "SOC_5cm.tif", "CEC_5cm.tif")

# Create soilgrid raster stack
soilgrid_stack <- rast(paste0("data/georeferenced_covariates/soil_grid/", soil_vars))

# Name the layers
names(soilgrid_stack) <- c("N_total", "P_available", 
                            "pH_CaCl2", "SOC", "CEC")

# Process bioclim stack
# Project temperate forest vector to bioclim CRS
temperate_forest_bioclim <- project(temperate_forest_vect, bioclim_stack)

# Crop to temperate forest ecoregions
bioclim_cropped <- crop(bioclim_stack, temperate_forest_bioclim)

# Project to forest_rast_trimmed CRS and mask
bioclim_trimmed <- project(bioclim_cropped, forest_rast_trimmed) %>%
  mask(forest_rast_trimmed)

# Process soilgrid stack
# Project temperate forest vector to soilgrid CRS
temperate_forest_soil <- project(temperate_forest_vect, soilgrid_stack)

# Crop to temperate forest ecoregions
soilgrid_cropped <- crop(soilgrid_stack, temperate_forest_soil)

# Project to forest_rast_trimmed CRS and mask
soilgrid_trimmed <- project(soilgrid_cropped, forest_rast_trimmed) %>%
  mask(forest_rast_trimmed)

# Get all cell numbers from forest_rast_trimmed
all_forest_cells <- cells(forest_rast_trimmed)

# Get cell center coordinates for all forest cells
all_cell_centers <- xyFromCell(forest_rast_trimmed, all_forest_cells) %>%
  as.data.frame() %>%
  mutate(cell = all_forest_cells)

# Convert cell centers back to longitude/latitude
all_cell_centers_longlat <- all_cell_centers %>%
  select(x, y) %>%
  as.matrix() %>%
  vect(crs = crs(forest_rast_trimmed)) %>%
  project('+proj=longlat') %>%
  crds() %>%
  as.data.frame() %>%
  rename(longitude = x, latitude = y) %>%
  bind_cols(cell = all_cell_centers$cell)

# Extract bioclim variables for all forest cells
bioclim_values_all <- terra::extract(
  bioclim_trimmed,
  all_cell_centers %>% select(x, y),
  ID = FALSE
)

# Extract soilgrid variables for all forest cells
soilgrid_values_all <- terra::extract(
  soilgrid_trimmed,
  all_cell_centers %>% select(x, y),
  ID = FALSE
)

# Create data frame with all forest cells and covariates
all_cells_with_covariates <- all_cell_centers_longlat %>%
  bind_cols(bioclim_values_all) %>%
  bind_cols(soilgrid_values_all) %>%
  # Drop NAs
  drop_na()

# Identify cells that are in mycorrhizal_proportions_by_cell (estimation data)
est_cells <- mycorrhizal_proportions_by_cell$cell

# Create the master tibble
master_tibble <- all_cells_with_covariates %>%
  # Add model_type column
  mutate(
    model_type = if_else(cell %in% est_cells, "est", "pred")
  ) %>%
  # Left join with mycorrhizal data (will add NAs for pred cells)
  left_join(
    mycorrhizal_proportions_by_cell,
    by = "cell",
    suffix = c("", "_myc")
  ) %>%
  # Rearrange columns to put model_type and mycorrhizal variables first
  select(
    cell, longitude, latitude, model_type,
    richness, total_abundance, total_diameter,
    AM_abundance, EcM_abundance, AM_abundance_prop, EcM_abundance_prop,
    AM_diameter, EcM_diameter, AM_diameter_prop, EcM_diameter_prop,
    everything()
  )

# Check results
message("Total number of forest cells: ", nrow(master_tibble))
message("Estimation cells (est): ", sum(master_tibble$model_type == "est"))
message("Prediction cells (pred): ", sum(master_tibble$model_type == "pred"))
glimpse(master_tibble)
master_tibble %>%
  filter(model_type == "est") %>%
  glimpse(.)

# (5) Compute PCs for climate and soil variables -------------------------------

# Check the distribution of bioclim variables
all_cells_with_covariates %>%
  select(-longitude, -latitude, -cell) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  ggplot(aes(x = value)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

# Variables to not transform
no_transform_vars <- c(
  "bio1", "bio15", "bio4", "bio5", "bio6"
)

all_cells_with_covariates %>%
  select(-longitude, -latitude, -cell, -all_of(no_transform_vars)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  ggplot(aes(x = value)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal() +
  scale_x_log10()

# Variables to log transform
log_transform_vars <- c(
  "bio12", "bio16", "bio17", "CEC", "N_total", "P_available", "SOC", 
  "pH_CaCl2", "SOC")

# Compute the PCA
pca <- all_cells_with_covariates %>%
  # Select only covariate columns
  select(-longitude, -latitude) %>%
  column_to_rownames("cell") %>%
  mutate(
    # Log transform certain variables
    across(all_of(log_transform_vars), ~ log10(.)),
    # Scale all variables
    across(everything(), ~ scale(.))
    ) %>%
  psych::principal(
    .,
    nfactors = 7,
    rotate = "none"
  )

# Grab the first four components
pca_scores <- as.data.frame(pca$scores) %>%
  rownames_to_column("cell") %>%
  select(cell, PC1, PC2, PC3, PC4)

# Join PCA scores to master tibble
master_tibble_final <- master_tibble %>%
  # Make "cell" column character to match pca_scores
  mutate(
    cell = as.character(cell)
  ) %>%
  left_join(
    pca_scores,
    by = "cell"
  )

# (5) Save files ---------------------------------------------------------------

# Save the master tibble
fwrite(
  master_tibble_final,
  "data/tree_mycorrhizal_dominance.txt",
  sep = "\t"
)

# Save temperate forest raster
writeRaster(
  forest_rast_trimmed,
  "data/temperate_forest_raster.tif",
  overwrite = TRUE
)

# Save the tempersate forest vector
# Dissolve into a single polygon
temperate_forest_vect <- aggregate(temperate_forest_vect, dissolve = TRUE)

# Check the result
plot(temperate_forest_vect)

# Save as GeoPackage (more efficient than shapefile)
writeVector(
  temperate_forest_vect,
  "data/temperate_forest_ecoregions.gpkg",
  overwrite = TRUE
)
writeVector(
  temperate_forest_vect,
  "data/temperate_forest_ecoregions.shp",
  overwrite = TRUE
)
