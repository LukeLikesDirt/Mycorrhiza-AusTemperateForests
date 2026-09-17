# Required packages and functions
require(statmod)
require(tweedie)
require(INLA) # inla.list.models()
require(fmesher)
require(terra)
require(ggtext)
require(tidyverse)
source("code/functions_inla.R")

# Set theme
tag_size <- 14
strip_size <- 12
title_size <- 10
text_size <- 9
common_theme <- theme_minimal() +
  theme(
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    panel.grid = element_line(colour = "grey90", linewidth = 0.1),
    axis.ticks = element_blank(),
    axis.text = element_text(size = text_size),
    axis.title = element_text(size = title_size),
    plot.title = element_text(face = "bold", size = title_size, hjust = 0.5),
    plot.tag = element_markdown(size = tag_size),
    strip.text = element_text(face = "bold", size = strip_size),
    plot.margin = margin(1, 1, 1, 1, "pt")
  )

# (1) Read data ################################################################

# Forest raster
forest_raster <- rast("data/temperate_forest_raster.tif")
forest_vector <- vect("data/temperate_forest_ecoregions.shp")

# Read the estimation data
data_est <- data.table::fread(
  "data/tree_mycorrhizal_richness.txt",
  stringsAsFactors = TRUE
) %>%
  na.omit() %>%
  filter(
    model_type == "est",
    n_obs >= 5,
    n_obs <= 500
  ) %>%
  mutate(
    row_id = row_number(),
    n_obs_std = (log(n_obs) - mean(log(n_obs)) / sd(log(n_obs))) # Standardise n_obs
  )

# Check the distribution of sampling effort
data_est %>%
  ggplot(aes(n_obs_std)) +
  geom_density(fill = "lightblue")

# Get the data with too low or too high sampling effort to add to prediction data
data_est_outliers_1 <- data.table::fread(
  "data/tree_mycorrhizal_richness.txt",
  stringsAsFactors = TRUE
) %>%
  na.omit() %>%
  filter(
    n_obs < 5 | n_obs > 500
  ) %>%
  mutate(
    model_type = "pred"
  ) %>%
  select(cell, longitude,  latitude, PC1, PC2, PC3, PC4, model_type)

# Get the outliers detected during model fitting
source("code/richness_outlier_ids.R")
data_est_outliers_2 <- data_est %>%
  filter(row_id %in% outlier_ids) %>%
  mutate(
    model_type = "pred"
  ) %>%
  select(cell, longitude,  latitude, PC1, PC2, PC3, PC4, model_type)

# Update the estimation data to remove outliers
data_est <- data_est %>%
  filter(!row_id %in% outlier_ids)

# Count the number of cells for estimation
n_est <- nrow(data_est) %>%
  print()

# Create estimation covariates with n_obs fixed to median to predict richness
# across prediction sites while accounting for sampling effort
data_est_marg <- data_est %>%
  mutate(n_obs = median(data_est$n_obs_std)) %>%
  select(n_obs, PC1, PC2, PC3, PC4)

# Read in the covariates for prediction
data_pred <- data.table::fread(
  "data/tree_mycorrhizal_dominance.txt",
  stringsAsFactors = TRUE
) %>%
  filter(model_type == "pred") %>%
  select(cell, longitude,  latitude, PC1, PC2, PC3, PC4, model_type) %>%
  # Add the outliers from estimation data
  bind_rows(data_est_outliers_1) %>%
  bind_rows(data_est_outliers_2)

# Count the number of cells for prediction
n_pred <- nrow(data_pred) %>%
  print()

# NOTE: I should have saved Albers coords in the data for INLA modelling

# Get x and y coordinates in the raster's CRS from cell IDs
cell_coords <- xyFromCell(forest_raster, data_est$cell) %>%
  as.data.frame() %>%
  rename(x_albers = x, y_albers = y)

# Add to your estimation data
data_est <- data_est %>%
  bind_cols(cell_coords)

# Same for prediction data
cell_coords_pred <- xyFromCell(forest_raster, data_pred$cell) %>%
  as.data.frame() %>%
  rename(x_albers = x, y_albers = y)

data_pred <- data_pred %>%
  bind_cols(cell_coords_pred)

# Mutate to KMs
data_est <- data_est %>%
  mutate(
    x_km = x_albers / 1000,
    y_km = y_albers / 1000
  )
data_pred <- data_pred %>%
  mutate(
    x_km = x_albers / 1000,
    y_km = y_albers / 1000
  )

# (2) Model formulas ###########################################################

# Ensure relationships are linear
data_est %>%
  select(PC1, PC2, PC3, PC4, am_tree_richness) %>%
  # Add 1 to richness to avoid log(0)
  mutate(am_tree_richness = am_tree_richness + 1) %>%
  # Sample 10000 points
  slice_sample(n = 10000) %>%
  pivot_longer(cols = -am_tree_richness) %>%
  ggplot(aes(value, am_tree_richness)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_log10() +
  facet_wrap(~ name, scales = "free")
data_est %>%
  select(PC1, PC2, PC3, PC4, ecm_tree_richness) %>%
  # Add 1 to richness to avoid log(0)
  mutate(ecm_tree_richness = ecm_tree_richness + 1) %>%
  # Sample 1000 points
  slice_sample(n = 10000) %>%
  select(PC1, PC2, PC3, PC4, ecm_tree_richness) %>%
  pivot_longer(cols = -ecm_tree_richness) %>%
  ggplot(aes(value, ecm_tree_richness)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_y_log10() +
  facet_wrap(~ name, scales = "free")

# Define the formulas
formula_AM <- as.formula(paste0(
  "y_AM ~ -1 + intercept + n_obs + ",
  "PC1 + PC2 + PC3 + PC4 +",
  "f(spatial, model = spde)"
))
formula_EcM <- as.formula(paste0(
  "y_EcM ~ -1 + intercept + n_obs + ",
  "PC1 + PC2 + PC3 + PC4 +",
  "f(spatial, model = spde)"
))

# (3) Prepare the mesh #########################################################

# Main tools to control the mesh:
#   (1) 'max.edge'    The largest allowed triangle length.
#                     The lower the value for 'max.edge' the higher the resolution. 
#   (2) 'quantileoff' Sites with a distance value smaller than the quantileoff 
#                     are modelled by a single vertex in the mesh.
#   (3) 'boundary'    Useful for islands and fjords.
#   (4) 'offset'      Determines how far the inner and outer boundaries extend.

# Site locations for the estimation and prediction data
loc_est <- data_est %>%
  select(x_km, y_km) %>%
  as.matrix()
loc_pred <- data_pred %>%
  select(x_km, y_km) %>%
  as.matrix()

# Distances between estimation site
small_scale <- 200

# Maximun edge for the triangulation (about 1/3 of the spatial correlation range)
max_edge <- small_scale / 5

# Create convex hull from estimation locations
nonconvex_hull <- fm_nonconvex_hull(
  loc_est %>%
    as.data.frame() %>%
    # Sample 1000 points to speed up
    slice_sample(n = 1000) %>%
    as.matrix(),
  convex = -0.1
)

# Prepare the mesh
mesh <- fm_mesh_2d_inla(
  loc = loc_est,
  max.edge = c(1, 5) * max_edge,
  boundary = nonconvex_hull,
  cutoff = max_edge / 5
)

# Define the mesh CSR
fm_crs(mesh) <- fm_CRS(crs(forest_vector))

# Extract geometry and scale to km
geom_data <- geom(forest_vector)
geom_data[, c("x", "y")] <- geom_data[, c("x", "y")] / 1000

# Rebuild the vector
forest_vector_km <- vect(geom_data, type = "polygons", crs = crs(forest_vector))

# Visualise the mesh using ggplot
ggplot() +
  theme_minimal() +
  fmesher::geom_fm(data = mesh) +
  geom_sf(
    data = forest_vector_km,
    fill = "transparent",
    col = "red"
  ) +
  geom_point(
    data = data_est %>%
      # Sample 10% to speed up plotting
      slice_sample(prop = 0.1),
    aes(x_km, y_km), size = 0.1, alpha = 0.1
  ) +
  labs(
    x = NULL,
    y = NULL
  )

# Define projector matrix for the estimation and prediction data
A_est <- inla.spde.make.A(
  mesh = mesh,
  loc = loc_est
)
dim(A_est)
A_pred <- inla.spde.make.A(
  mesh = mesh,
  loc = loc_pred
)
dim(A_pred)

# (4) Define the SPDE ##########################################################

# Range and sigma estimations for the priors
range_est <- 60
sigma_est <- 0.75

# Define the SPDE for the estimation
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  # Range of the spatial correlation is unlikely to be less range_est at p=0.01
  prior.range = c(range_est, 0.01),
  # Prior for the standard deviation is unlikely to be larger than sigma_est at p=0.01
  prior.sigma = c(sigma_est, 0.01)
)

# Define the SPDE index
index_spatial <- inla.spde.make.index(
  name = "spatial",
  n.spde = spde$n.spde
)

# The size of the mesh and n.spde should match
mesh$n
spde$n.spde

# (5) Stack the data ###########################################################

# Covariates for the estimation and prediction (spatial and marginal effects)
covariates <- bind_rows(
  tibble(
    # Observed n_obs for estimation
    n_obs = c(data_est$n_obs_std, rep(median(data_est$n_obs_std), n_pred)),
    PC1 = c(data_est$PC1, data_pred$PC1),
    PC2 = c(data_est$PC2, data_pred$PC2),
    PC3 = c(data_est$PC3, data_pred$PC3),
    PC4 = c(data_est$PC4, data_pred$PC4),
    tag = c(rep("est", n_est), rep("pred", n_pred))
  ),
  # Marginalised estimation (estimation sites with n_obs = median)
  tibble(
    n_obs = rep(median(data_est$n_obs_std), n_est),
    PC1 = data_est$PC1,
    PC2 = data_est$PC2,
    PC3 = data_est$PC3,
    PC4 = data_est$PC4,
    tag = "est_marg"
  ),
  tibble(
    n_obs = seq(min(data_est$n_obs_std), max(data_est$n_obs_std), length.out = 100),
    PC1 = median(data_est$PC1),
    PC2 = median(data_est$PC2),
    PC3 = median(data_est$PC3),
    PC4 = median(data_est$PC4),
    tag = "pred_n_obs"
  ),
  tibble(
    n_obs = median(data_est$n_obs_std),
    PC1 = seq(
      min(data_est$PC1), max(data_est$PC1), length.out = 100
    ),
    PC2 = median(data_est$PC2),
    PC3 = median(data_est$PC3),
    PC4 = median(data_est$PC4),
    tag = "pred_PC1"
  ),
  tibble(
    n_obs = median(data_est$n_obs_std),
    PC1 = median(data_est$PC1),
    PC2 = seq(
      min(data_est$PC2), max(data_est$PC2), length.out = 100
    ),
    PC3 = median(data_est$PC3),
    PC4 = median(data_est$PC4),
    tag = "pred_PC2"
  ),
  tibble(
    n_obs = median(data_est$n_obs_std),
    PC1 = median(data_est$PC1),
    PC2 = median(data_est$PC2),
    PC3 = seq(
      min(data_est$PC3), max(data_est$PC3), length.out = 100
    ),
    PC4 = median(data_est$PC4),
    tag = "pred_PC3"
  ),
  tibble(
    n_obs = median(data_est$n_obs_std),
    PC1 = median(data_est$PC1),
    PC2 = median(data_est$PC2),
    PC3 = median(data_est$PC3),
    PC4 = seq(
      min(data_est$PC4), max(data_est$PC4), length.out = 100
    ),
    tag = "pred_PC4"
  )
)

# Stack the estimation data
stack_est <- inla.stack(
  tag = "est",
  data = list(
    #y_all = data_est$absolute_tree_richness,
    y_AM = data_est$am_tree_richness,
    y_EcM = data_est$ecm_tree_richness
  ),
  A = list(1, 1, A_est),
  effects = list(
    intercept = rep(1, n_est),
    X = covariates %>% filter(tag == "est") %>% select(-tag),
    spatial = index_spatial
  )
)

stack_est_marg <- inla.stack(
  tag = "est_marg",
  data = list(
    y_AM = rep(NA_real_, n_est),
    y_EcM = rep(NA_real_, n_est)
  ),
  A = list(1, 1, A_est),
  effects = list(
    intercept = rep(1, n_est),
    X = covariates %>% filter(tag == "est_marg") %>% select(-tag),
    spatial = index_spatial
  )
)

# Prediction data (new cells)
stack_pred <- inla.stack(
  tag = "pred",
  data = list(
    y_AM = rep(NA_real_, n_pred),
    y_EcM = rep(NA_real_, n_pred)
  ),
  A = list(1, 1, A_pred),
  effects = list(
    intercept = rep(1, n_pred),
    X = covariates %>% filter(tag == "pred") %>% select(-tag),
    spatial = index_spatial
  )
)

# Stack the data for the marginal effects
stack_PC1 <- inla.stack(
  tag = "pred_PC1",
  data = list(
    y_AM = rep(NA_real_, 100),
    y_EcM = rep(NA_real_, 100)
  ),
  # The order needs to matches the order in effects
  A = list(1, 1),
  effects = list(
    intercept = rep(1, 100),
    X = covariates %>% filter(tag == "pred_PC1") %>% select(-tag)
  )
)
stack_PC2 <- inla.stack(
  tag = "pred_PC2",
  data = list(
    y_AM = rep(NA_real_, 100),
    y_EcM = rep(NA_real_, 100)
  ),
  # The order needs to matches the order in effects
  A = list(1, 1),
  effects = list(
    intercept = rep(1, 100),
    X = covariates %>% filter(tag == "pred_PC2") %>% select(-tag)
  )
)
stack_PC3 <- inla.stack(
  tag = "pred_PC3",
  data = list(
    y_AM = rep(NA_real_, 100),
    y_EcM = rep(NA_real_, 100)
  ),
  # The order needs to matches the order in effects
  A = list(1, 1),
  effects = list(
    intercept = rep(1, 100),
    X = covariates %>% filter(tag == "pred_PC3") %>% select(-tag)
  )
)
stack_PC4 <- inla.stack(
  tag = "pred_PC4",
  data = list(
    y_AM = rep(NA_real_, 100),
    y_EcM = rep(NA_real_, 100)
  ),
  # The order needs to matches the order in effects
  A = list(1, 1),
  effects = list(
    intercept = rep(1, 100),
    X = covariates %>% filter(tag == "pred_PC4") %>% select(-tag)
  )
)
stack_n_obs <- inla.stack(
  tag = "pred_n_obs",
  data = list(
    y_AM = rep(NA_real_, 100),
    y_EcM = rep(NA_real_, 100)
  ),
  # The order needs to matches the order in effects
  A = list(1, 1),
  effects = list(
    intercept = rep(1, 100),
    X = covariates %>% filter(tag == "pred_n_obs") %>% select(-tag)
  )
)

# Join the stacks
stack <- inla.stack.join(
  stack_est,
  stack_est_marg,
  stack_pred,
  stack_PC1,
  stack_PC2,
  stack_PC3,
  stack_PC4,
  stack_n_obs
)

# Retrieve the stack indexes
index_est <- inla.stack.index(stack, "est")$data
index_est_marg <- inla.stack.index(stack, "est_marg")$data
index_pred <- inla.stack.index(stack, "pred")$data
index_PC1 <- inla.stack.index(stack, "pred_PC1")$data
index_PC2 <- inla.stack.index(stack, "pred_PC2")$data
index_PC3 <- inla.stack.index(stack, "pred_PC3")$data
index_PC4 <- inla.stack.index(stack, "pred_PC4")$data
index_n_obs <- inla.stack.index(stack, "pred_n_obs")$data

# (6) Fit the models ###########################################################

# AM
model_AM <- inla(
  formula_AM,
  family = 'tweedie',
  data = inla.stack.data(stack),
  control.predictor = list(
    compute = TRUE, link = 1,
    A = inla.stack.A(stack)
  ),
  control.compute = control_compute
)

# EcM
model_EcM <- inla(
  formula_EcM,
  family = 'tweedie',
  data = inla.stack.data(stack),
  control.predictor = list(
    compute = TRUE, link = 1,
    A = inla.stack.A(stack)
  ),
  control.compute = control_compute
)

# Check the summaries
summary(model_AM)
summary(model_EcM)

# (7) Model diagnostics ######################################################

# Compare sigma
posterior_marginals(model_AM, model_EcM)

#### * Compute diagnostic data * ####

# Diagnostic data: Compute fitted values, residuals, and posterior predictive 
# checks
data_diagnostics <- inla_spat_diagnostic_data(
  model_names = c("model_AM", "model_EcM"),
  stack_names = c("stack_est", "stack_est"),
  observed_values = list(
    data_est$am_tree_richness,
    data_est$ecm_tree_richness
  ),
  posterior_sampling = FALSE,
  observation_indices = list(
    data_est$row_id,
    data_est$row_id
  )
)

# Extract parameters from the Tweedie models
# Dispersion parameter
disp_AM  <- model_AM$summary.hyperpar["Dispersion parameter for Tweedie", "mean"]
disp_EcM <- model_EcM$summary.hyperpar["Dispersion parameter for Tweedie", "mean"]

# Power parameter (usually called 'p' or 'power')
# Check the exact name in your model
power_AM  <- model_AM$summary.hyperpar["p parameter for Tweedie", "mean"]
power_EcM <- model_EcM$summary.hyperpar["p parameter for Tweedie", "mean"]

# If power is fixed (not estimated), you need to extract it from your model formula
# For example, if you specified family = 'tweedie' with a fixed power:
# power_AM <- 1.5  # or whatever you specified

# Add parameters to your diagnostic data
data_diagnostics <- data_diagnostics %>%
  mutate(
    disp = case_when(
      model == "Model AM" ~ disp_AM,
      model == "Model EcM" ~ disp_EcM,
      TRUE ~ NA_real_
    ),
    power = case_when(
      model == "Model AM" ~ power_AM,
      model == "Model EcM" ~ power_EcM,
      TRUE ~ NA_real_
    )
  )

# Now calculate randomised quantile residuals
data_diagnostics <- data_diagnostics %>%
  rowwise() %>%
  mutate(
    qresiduals = qnorm(
      ptweedie(
        q = observed,
        mu = fitted,
        phi = disp,
        power = power
      )
    )
  ) %>%
  ungroup()

#### * Q-Q plot * ####
qq_plot <- ggplot(
  data_diagnostics %>%
    # Select 1000 points to speed up plotting
    group_by(model) %>%
    slice_sample(n = 1000) %>%
    ungroup(), aes(sample = qresiduals)) +
  stat_qq(alpha = 0.3, shape = 1) +
  stat_qq_line(colour = "red") +
  facet_wrap(~ model) +
  labs(x = "Theoretical Quantiles", y = "Quantile Residuals") +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    aspect.ratio = 1
  )

# Display the plot
print(qq_plot)

#### * Homogeneity of variance * ####

homogeneity_plot <- ggplot(
  data_diagnostics %>%
    filter(fitted < 40) %>%
    # Select 1000 points to speed up plotting
    group_by(model) %>%
    slice_sample(n = 1000),
  aes(x = fitted, y = qresiduals)) +
  geom_point(alpha = 0.3, shape = 1) +
  geom_smooth(method = "loess", colour = "red") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  facet_wrap(~ model, nrow = 1, scales = "free_x") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    aspect.ratio = 1,
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, "pt"),
    strip.text = element_blank()
  ) +
  labs(tag = "b", x = "Fitted values", y = "Residuals")

# Display the plot
print(homogeneity_plot)

#### * Observed vs predicted * ####
# limits_list <- list(
#   "Model AM" = c(
#     min(data_diagnostics %>% filter(model == "Model AM") %>% pull(observed), 
#         data_diagnostics %>% filter(model == "Model AM") %>% pull(fitted)),
#     max(data_diagnostics %>% filter(model == "Model AM") %>% pull(observed), 
#         data_diagnostics %>% filter(model == "Model AM") %>% pull(fitted))
#   ),
#   "Model EcM" = c(
#     min(data_diagnostics %>% filter(model == "Model EcM") %>% pull(observed), 
#         data_diagnostics %>% filter(model == "Model EcM") %>% pull(fitted)),
#     max(data_diagnostics %>% filter(model == "Model EcM") %>% pull(observed), 
#         data_diagnostics %>% filter(model == "Model EcM") %>% pull(fitted))
#   )
# )

# Observed vs predicted plot
obs_fit_plot <- ggplot(
  data_diagnostics %>%
    # Select 1000 points to speed up plotting
    group_by(model) %>%
    slice_sample(n = 10000) %>%
    ungroup(),
  aes(x = fitted, y = observed)
) +
  geom_point(alpha = 0.3, shape = 1, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, colour = "red") +
  ggpubr::stat_cor(aes(label = after_stat(rr.label)), colour = "red", size = 3) +
  facet_wrap(~ model, nrow = 1, scales = "free") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    aspect.ratio = 1,
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, "pt")
  ) +
  labs(x = "Posterior mean predictions", y = "Observed values")
# +
#   ggh4x::facetted_pos_scales(
#     x = list(
#       "Model AM" = scale_x_continuous(limits = limits_list[["Model AM"]]),
#       "Model EcM" = scale_x_continuous(limits = limits_list[["Model EcM"]])
#     ),
#     y = list(
#       "Model AM" = scale_y_continuous(limits = limits_list[["Model AM"]]),
#       "Model EcM" = scale_y_continuous(limits = limits_list[["Model EcM"]])
#     )
#   )

# Display the plot
print(obs_fit_plot)

#### * Calibration plot * ####

data_diagnostics %>%
  group_by(model) %>%
  mutate(fitted_bin = cut(fitted, breaks = 20)) %>%
  group_by(model, fitted_bin) %>%
  summarise(
    mean_observed = mean(observed),
    mean_fitted = mean(fitted),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = mean_fitted, y = mean_observed)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  facet_wrap(~ model, scales = "free") +
  labs(x = "Mean Fitted", y = "Mean Observed") +
  theme_minimal() +
  theme(aspect.ratio = 1)

#### * Residuals vs predictors * ####

# AM
tibble(
  residuals = data_diagnostics %>% filter(model == "Model AM") %>% pull("qresiduals"),
  PC1 = data_est[["PC1"]],
  PC2 = data_est[["PC2"]],
  PC3 = data_est[["PC3"]],
  PC4 = data_est[["PC4"]]
) %>%
  # Sample 1000 points to speed up plotting
  sample_n(1000) %>%
  pivot_longer(
    -residuals
  ) %>%
  ggplot(aes(x = value, y = residuals)) +
  geom_point() +
  stat_smooth(method = "loess") +
  facet_wrap(~ name, scales = "free_x") +
  theme_minimal() +
  theme(aspect.ratio = 1)

# EcM
tibble(
  residuals = data_diagnostics %>% filter(model == "Model EcM") %>% pull("qresiduals"),
  PC1 = data_est[["PC1"]],
  PC2 = data_est[["PC2"]],
  PC3 = data_est[["PC3"]],
  PC4 = data_est[["PC4"]]
) %>%
  sample_n(1000) %>%
  pivot_longer(
    -residuals
  ) %>%
  ggplot(aes(x = value, y = residuals)) +
  geom_point() +
  stat_smooth(method = "loess") +
  facet_wrap(~ name, scales = "free_x") +
  theme_minimal() +
  theme(aspect.ratio = 1)

#### * Check outliers * ####
# Identify outliers based on quantile residuals
# Values beyond ±3 are typically considered outliers
data_diagnostics <- data_diagnostics %>%
  mutate(
    outlier_qres = abs(qresiduals) > 3,
    outlier_severity = case_when(
      abs(qresiduals) > 4 ~ "Extreme",
      abs(qresiduals) > 3 ~ "Moderate",
      TRUE ~ "Normal"
    )
  )

# Summary of outliers by model
data_diagnostics %>%
  group_by(model) %>%
  summarise(
    n_total = n(),
    n_moderate_outliers = sum(abs(qresiduals) > 3, na.rm = TRUE),
    n_extreme_outliers = sum(abs(qresiduals) > 4, na.rm = TRUE),
    pct_outliers = round(100 * n_moderate_outliers / n_total, 2)
  )

# !!! Run when removing outliers from the full dataset !!!

# # View the actual outliers
# outliers <- data_diagnostics %>%
#   filter(abs(qresiduals) > 3) %>%
#   arrange(desc(abs(qresiduals))) %>%
#   print()
# 
# # Get unique outlier row_ids (since each appears twice - once per model)
# outlier_indices <- outliers %>%
#   pull(obs_index) %>%
#   unique() %>%
#   sort()
# 
# # Print in a format you can copy-paste
# cat("outlier_ids <- c(", paste(outlier_indices, collapse = ", "), ")\n")
# 
# # Save to a file for future reference
# writeLines(
#   c("# Outlier row_ids identified on [DATE]",
#     paste0("outlier_ids <- c(", paste(outlier_indices, collapse = ", "), ")")),
#   "code/richness_outlier_ids.R"
# )

# (8) Predicted values #########################################################

#### * Predicted values * ####

# Define the model names and get observed values
model_info <- list(
  AM = list(
    name = "model_AM",
    observed = data_est$am_tree_richness
  ),
  EcM = list(
    name = "model_EcM",
    observed = data_est$ecm_tree_richness
  )
)

# Create an empty list to store the predicted values for each model
fitted_values_list <- list()

# Loop through each model
for (model_type in names(model_info)) {
  # Get the current model object
  current_model <- get(model_info[[model_type]]$name)
  
  # 1. Predictions for residuals (observed n_obs)
  predicted_resid <- current_model$summary.fitted$mean[index_est]
  
  # 2. Predictions for mapping (marginalised)
  predicted_map <- c(
    current_model$summary.fitted$mean[index_est_marg],
    current_model$summary.fitted$mean[index_pred]
  )
  
  current_fitted <- tibble(
    model = model_type,
    predicted = predicted_map,  # <-- mapping
    upper = c(
      current_model$summary.fitted$`0.975quant`[index_est_marg],
      current_model$summary.fitted$`0.975quant`[index_pred]
    ),
    lower = c(
      current_model$summary.fitted$`0.025quant`[index_est_marg],
      current_model$summary.fitted$`0.025quant`[index_pred]
    ),
    sd = c(
      current_model$summary.fitted$sd[index_est_marg],
      current_model$summary.fitted$sd[index_pred]
    ),
    observed = c(
      model_info[[model_type]]$observed,
      rep(NA_integer_, n_pred)
    ),
    # Residuals should come from observed-fit, not observed-marginal
    residuals = c(
      model_info[[model_type]]$observed - predicted_resid,
      rep(NA_real_, n_pred)
    ),
    pearson_residuals = c(
      (model_info[[model_type]]$observed - predicted_resid) /
        current_model$summary.fitted$sd[index_est],
      rep(NA_real_, n_pred)
    ),
    latitude = c(data_est$latitude, data_pred$latitude),
    longitude = c(data_est$longitude, data_pred$longitude),
    x_albers = c(data_est$x_albers, data_pred$x_albers),
    y_albers = c(data_est$y_albers, data_pred$y_albers)
  ) %>%
    mutate(
      uncertainty = upper - lower,
      x_albers = as.numeric(x_albers),
      y_albers = as.numeric(y_albers)
    )
  
  # Spatial field (same as before)
  A_posterior <- fmesher::fm_basis(
    mesh,
    loc = matrix(
      c(current_fitted$x_albers / 1000, current_fitted$y_albers / 1000),
      ncol = 2
    )
  )
  
  spatial_field <- tibble(
    spatial_field = as.vector(A_posterior %*% current_model$summary.random$spatial[, "mean"]),
    x_albers = current_fitted$x_albers,
    y_albers = current_fitted$y_albers
  ) %>%
    unique()
  
  fitted_values_list[[model_type]] <- current_fitted %>%
    left_join(spatial_field, by = c("x_albers", "y_albers"))
}

# Combine all models into one dataframe
fitted_values <- bind_rows(fitted_values_list)

#### * Rasterise predicted values * ####

# AM raster
fitted_rast_AM <- rasterize(
  vect(
    fitted_values %>% filter(model == "AM"),
    geom = c("x_albers", "y_albers"),
    crs = crs(forest_raster)
  ),
  forest_raster,
  field = names(fitted_values %>% select(-latitude, -longitude, -x_albers, -y_albers)),
  background = NA
)

# EcM raster
fitted_rast_EcM <- rasterize(
  vect(
    fitted_values %>% filter(model == "EcM"),
    geom = c("x_albers", "y_albers"),
    crs = crs(forest_raster)
  ),
  forest_raster,
  field = names(fitted_values %>% select(-latitude, -longitude, -x_albers, -y_albers)),
  background = NA
)

# Standardised richness to percentiles
standardised_rast_AM <- standardise_raster_percentiles(fitted_rast_AM[["predicted"]])
standardised_rast_EcM <- standardise_raster_percentiles(fitted_rast_EcM[["predicted"]])

#### * Plot predicted values * ####

# Plot the quantile rasters
myco_richness_plot_quant <- ggplot() +
  # Plot the forest vector
  tidyterra::geom_spatvector(
    data = forest_vector,
    fill = "grey90",
    col = "grey90",
    linewidth = 0.25
  ) +
  tidyterra::geom_spatraster(
    data = terra::rast(list(
      AM = standardised_rast_AM,
      EcM = standardised_rast_EcM
    )),
    na.rm = TRUE
  ) +
  scale_fill_stepsn(
    colors = rev(paletteer::paletteer_c("grDevices::Spectral", 100)),
    breaks = 1:100,
    labels = c(rep("", 4), "Low", rep("", 90), "High", rep("", 4)),
    na.value = "transparent"
  ) +
  facet_wrap(~ lyr, nrow = 1) +
  labs(
    tag = "(**a**)",
    fill = "Tree richness"
  ) +
  common_theme +
  # Adjust the tag to account fror the strip text
  theme(
    plot.tag.position = c(0.02, 0.95),
    axis.text.x = element_blank()
  )

#### * Hot spot plot * ####

myco_richness_plot_hotspot <- ggplot() +
  # Plot the forest vector
  tidyterra::geom_spatvector(
    data = forest_vector,
    fill = "grey90",
    col = "grey90",
    linewidth = 0.25
  ) +
  # Add the hotspot rasters
  tidyterra::geom_spatraster(
    data = terra::rast(list(
      AM = standardised_rast_AM %>%
        tidyterra::mutate(predicted = ifelse(predicted >= 95, "1", "0")),
      EcM = standardised_rast_EcM %>%
        tidyterra::mutate(predicted = ifelse(predicted >= 95, "1", "0"))
    )),
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = c("1" = "#cb181d", "0" = "#4A4A4A"),
    breaks = c("1", "0"),
    labels = c("Hotspot", "Non-hotspot"),
    name = "Richness\nhotspots",
    na.value = "transparent"
  ) +
  facet_wrap(~lyr, nrow = 1) +
  labs(
    tag = "(**b**)"
  ) +
  common_theme +
  theme(
    strip.text = element_blank()
  )

# Join the two plots
myco_richness_combined_plot <- patchwork::wrap_plots(
  myco_richness_plot_quant,
  myco_richness_plot_hotspot,
  nrow = 2
)

# Save the plots
ggsave(
  filename = "output/supplementary/myco_tree_richness.png",
  plot = myco_richness_combined_plot,
  width = 16,
  height = 20,
  units = "cm",
  bg = "white",
  dpi = 300
)

# (9) Save the rasters ##########################################################

# Add predicted percentiles to the rasters
fitted_rast_AM$predicted_percentile <- standardised_rast_AM
fitted_rast_EcM$predicted_percentile <- standardised_rast_EcM

# Check the names
names(fitted_rast_AM)
names(fitted_rast_EcM)

# Create the data/covariates_georeferenced/ directory if it doesn't exist
if (!dir.exists("data/covariates_georeferenced/")) {
  dir.create("data/covariates_georeferenced/")
}

# Save the rasters
writeRaster(
  fitted_rast_AM,
  overwrite = TRUE,
  filename = "data/covariates_georeferenced/richness_AM_raster.tif"
)
writeRaster(
  fitted_rast_EcM,
  overwrite = TRUE,
  filename = "data/covariates_georeferenced/richness_EcM_raster.tif"
)
