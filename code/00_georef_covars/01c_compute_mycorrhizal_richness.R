
# Load libraries
require(ggtext)
require(terra)
require(data.table)
require(tidyverse)

# (1) Trim forest raster -------------------------------------------------------

# Forest raster
forest_rast <- rast(
  "data/covariates_georeferenced/native_forest.tif"
)

# Temperate forest ecoregions
temperate_forest_vect <- vect("data/covariates_georeferenced/Ecoregions2017/Ecoregions2017.shp") %>%
  tidyterra::filter(BIOME_NAME == "Temperate Broadleaf & Mixed Forests") %>%
  # Crop to Australia extent
  terra::crop(ext(110, 155, -45, -10))

# Ensure CRS match
temperate_forest_vect <- project(temperate_forest_vect, forest_rast)

# Crop and mask the raster to the ecoregion
forest_rast_trimmed <- crop(forest_rast, temperate_forest_vect, mask = TRUE)

# Check the result
plot(forest_rast_trimmed)

# (2) Trim tree presence plots ------------------------------------------------

# Mycorrhizal types
mycorrhizal_types <- fread("data/tree_mycorrhizal_types.txt") %>%
  rename(species = scientific_name)

# Extract unique sites with coordinates
all_sites <- fread("data/tree_presence.txt") %>%
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
forest_trees <- fread("data/tree_presence.txt") %>%
  rename(species = scientific_name) %>%
  filter(
    site %in% forest_sites$site
  ) %>%
  distinct() %>%
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

# Compute sample effort per grid cell
sample_effort <- forest_trees %>%
  left_join(
    forest_sites %>% select(site, cell),
    by = "site"
  ) %>%
  group_by(cell) %>%
  summarise(
    n_obs = n()
  ) %>%
  ungroup()

# Aggregate tree data by grid cell WITH CELL CENTER COORDINATES
forest_trees_by_cell <- forest_trees %>%
  left_join(
    forest_sites %>% select(site, cell),
    by = "site"
  ) %>%
  select(-longitude, -latitude) %>%  # Remove original coordinates
  left_join(cell_centers_longlat, by = "cell")  # Add cell center coordinates

# Compute mycorrhizal proportions per grid cell
mycorrhizal_richness_by_cell <- forest_trees_by_cell %>%
  select(cell, longitude, latitude, species, AM, EcM) %>%
  distinct() %>%
  group_by(cell, longitude, latitude) %>%
  summarise(
    absolute_tree_richness = n_distinct(species),
    am_tree_richness = sum(AM),
    ecm_tree_richness = sum(EcM)
  ) %>%
  ungroup() %>%
  left_join(sample_effort, by = "cell")

# Check results
message("Number of unique grid cells: ", nrow(mycorrhizal_richness_by_cell))
glimpse(mycorrhizal_richness_by_cell)

# (4) Read in the covariate rasters --------------------------------------------

# Covariates were generated in "01a.Compute_mycorrhizal_dominance.R", so I will
# use those but ensure they align with the current data first

# Create data frame with all forest cells and covariates
all_cells_with_covariates <- fread("data/tree_mycorrhizal_dominance.txt") %>%
  select("cell", "longitude",  "latitude", "PC1", "PC2", "PC3", "PC4")

# Check cells and long lat match - JOIN ONLY BY CELL
all_cells_with_covariates %>%
  inner_join(
    mycorrhizal_richness_by_cell,
    by = "cell",
    suffix = c("_saved", "_new")
  )

# Join by cell and check coordinate differences
comparison <- all_cells_with_covariates %>%
  inner_join(
    mycorrhizal_richness_by_cell,
    by = "cell",
    suffix = c("_saved", "_new")
  )

# Check if coordinates match within tolerance
comparison %>%
  mutate(
    lon_diff = abs(longitude_saved - longitude_new),
    lat_diff = abs(latitude_saved - latitude_new)
  ) %>%
  select(cell, lon_diff, lat_diff) %>%
  summary()

# Looks to be tiny floating-point precision differences, so will round to 4 
# decimal places (~11m accuracy at the equator)

# Round coordinates in both datasets before joining
all_cells_with_covariates_rounded <- all_cells_with_covariates %>%
  mutate(
    longitude = round(longitude, 4),
    latitude = round(latitude, 4)
  )

mycorrhizal_richness_by_cell_rounded <- mycorrhizal_richness_by_cell %>%
  mutate(
    longitude = round(longitude, 4),
    latitude = round(latitude, 4)
  )

# Check cells and long lat match: Sould be around 90k cells for estimation
all_cells_with_covariates_rounded %>%
  inner_join(
    mycorrhizal_richness_by_cell_rounded,
    by = c("cell", "longitude", "latitude")
  )

# Now join by all three
mycorrhizal_richness_with_covariates <- all_cells_with_covariates_rounded %>%
  full_join(
    mycorrhizal_richness_by_cell_rounded,
    by = c("cell", "longitude", "latitude")
  ) %>%
  # Add the "est" and "pred" model type
  mutate(
    model_type = if_else(
      is.na(absolute_tree_richness), "pred", "est"
    )
  )

# Check results
message("Total number of forest cells: ", nrow(mycorrhizal_richness_with_covariates))
message("Estimation cells (est): ", sum(mycorrhizal_richness_with_covariates$model_type == "est"))
message("Prediction cells (pred): ", sum(mycorrhizal_richness_with_covariates$model_type == "pred"))

# (5) Save files ---------------------------------------------------------------

# Save the master tibble
fwrite(
  mycorrhizal_richness_with_covariates,
  "data/tree_mycorrhizal_richness.txt",
  sep = "\t"
)
