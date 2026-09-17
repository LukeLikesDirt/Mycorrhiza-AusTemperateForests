# INLA model prams #############################################################

control_compute <- list(
  config = TRUE,
  waic = TRUE,
  dic = TRUE,
  residuals = TRUE
)

# INLA diagnostics #############################################################

inla_spat_diagnostic_data <- function(
    model_names, observed_values, stack_names, backtransform = "none",
    n_samples = 1000, posterior_sampling = TRUE, observation_indices = NULL) {
  
  # Input validation
  if (!is.character(model_names) || !is.character(stack_names)) {
    stop("model_names and stack_names must be character vectors")
  }
  if (length(model_names) != length(stack_names)) {
    stop("model_names and stack_names must have the same length")
  }
  if (!backtransform %in% c("none", "exp", "plogis")) {
    stop("backtransform must be one of: 'none', 'exp', 'plogis'")
  }
  if (length(observed_values) != length(model_names)) {
    stop("observed_values must have the same length as model_names")
  }
  if (!is.logical(posterior_sampling)) {
    stop("posterior_sampling must be TRUE or FALSE")
  }
  
  # Initialize vectors
  observed <- c()
  fitted <- c()
  upper <- c()
  lower <- c()
  residuals <- c()
  pearson_residuals <- c()
  posterior_mean <- c()
  posterior_lower <- c()
  posterior_upper <- c()
  models <- c()
  obs_index <- c() 
  
  # Loop through each model
  for (i in seq_along(model_names)) {
    # Get model object from its name with error handling
    model <- try(get(model_names[i]), silent = TRUE)
    if (inherits(model, "try-error")) {
      stop(paste("Could not find model:", model_names[i]))
    }
    
    # Extract model number from the suffix
    model_suffix <- sub(".*?_", "", model_names[i])
    
    # Define the stack name
    stack_name <- try(get(stack_names[i]), silent = TRUE)
    if (inherits(stack_name, "try-error")) {
      stop(paste("Could not find stack:", stack_names[i]))
    }
    
    # Estimation data indexes and number of observations
    index_est <- inla.stack.index(stack_name, "est")$data
    n_observations <- length(index_est)
    
    # Extract fitted values
    fitted_values <- model$summary.fitted.values[index_est, "mean"]
    lower_values <- model$summary.fitted.values[index_est, "0.025quant"]
    upper_values <- model$summary.fitted.values[index_est, "0.975quant"]
    
    # Get the corresponding observed values for this model
    current_observed <- observed_values[[i]][index_est]
    
    # Estimation data indexes and number of observations
    index_est <- inla.stack.index(stack_name, "est")$data
    n_observations <- length(index_est)
    
    # Calculate residuals
    residuals_i <- current_observed - fitted_values
    
    # Calculate pearson residuals
    pearson_residuals_i <- residuals_i / sqrt(fitted_values)
    
    # Generate posterior predictive samples (conditional)
    if (posterior_sampling) {
      tryCatch({
        set.seed(1986)
        posterior_samples <- inla.posterior.sample(n_samples, model)
        predictive_values <- sapply(posterior_samples, function(sample) {
          sample$latent[grep("Predictor", rownames(sample$latent))]
        })
        predictive_mean <- apply(predictive_values, 1, mean)
        predictive_quantiles <- apply(predictive_values, 1, quantile, probs = c(0.025, 0.5, 0.975))
      }, error = function(e) {
        stop(paste("Error in posterior sampling for model", model_names[i], ":", e$message))
      })
      
      # Apply backtransformation
      if (backtransform == "exp") {
        posterior_pred_df <- data.frame(
          posterior_mean = exp(predictive_mean[index_est]),
          lower = exp(predictive_quantiles[1, index_est]),
          upper = exp(predictive_quantiles[3, index_est])
        )
      } else if (backtransform == "plogis") {
        posterior_pred_df <- data.frame(
          posterior_mean = plogis(predictive_mean[index_est]),
          lower = plogis(predictive_quantiles[1, index_est]),
          upper = plogis(predictive_quantiles[3, index_est])
        )
      } else {
        posterior_pred_df <- data.frame(
          posterior_mean = predictive_mean[index_est],
          lower = predictive_quantiles[1, index_est],
          upper = predictive_quantiles[3, index_est]
        )
      }
    } else {
      # Set posterior columns to NA when posterior_sampling is FALSE
      posterior_pred_df <- data.frame(
        posterior_mean = rep(NA_real_, n_observations),
        lower = rep(NA_real_, n_observations),
        upper = rep(NA_real_, n_observations)
      )
    }
    
    # Append values to vectors
    observed <- c(observed, current_observed)
    fitted <- c(fitted, fitted_values)
    upper <- c(upper, upper_values)
    lower <- c(lower, lower_values)
    residuals <- c(residuals, residuals_i)
    pearson_residuals <- c(pearson_residuals, pearson_residuals_i)
    posterior_mean <- c(posterior_mean, posterior_pred_df$posterior_mean)
    posterior_lower <- c(posterior_lower, posterior_pred_df$lower)
    posterior_upper <- c(posterior_upper, posterior_pred_df$upper)
    models <- c(models, rep(paste0("Model ", model_suffix), n_observations))
    
    # Add observation indices
    if (!is.null(observation_indices)) {
      obs_index <- c(obs_index, observation_indices[[i]][index_est])
    } else {
      obs_index <- c(obs_index, index_est)
    }
  }
  
  # Create the tibble
  data_diagnostics <- tibble(
    obs_index = obs_index,  # Add this
    observed = observed,
    fitted = fitted,
    upper = upper,
    lower = lower,
    residuals = residuals,
    pearson_residuals = pearson_residuals,
    posterior_mean = posterior_mean,
    posterior_lower = posterior_lower,
    posterior_upper = posterior_upper,
    model = models
  )
  
  return(data_diagnostics)
}
# Posterior Marginal Plot ######################################################

# Define the custom theme
my_theme <- function() {
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_line(colour = 'black', linewidth = 0.25),
    plot.title = element_text(hjust = 0.5, vjust = 1, size = rel(0.9)),
    legend.position = c(1, 1),  # Place legend at the top-right corner
    legend.justification = c(1.1, 1.1),  # Justify legend to the top-right corner
    legend.box.just = "right",  # Align legend box to the right
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Set margin to 0
    legend.background = element_blank()
  )
}

# Support function for the labels
sigma_dist_labs <- function() {
  labs(
    x = expression(sigma),
    y = expression(paste("P(", sigma, " | Data)"))
  )
}

# Plotting function for posterior marginals of sigma

posterior_marginals <- function(...) {
  models <- list(...)
  
  marginals_list <- lapply(seq_along(models), function(i) {
    model <- models[[i]]
    model_name <- paste0("Model ", i)
    
    inla.tmarginal(
      fun = function(x) exp(-0.5 * x),
      marg = model$internal.marginals.hyperpar[[1]]
    ) %>%
      as_tibble() %>%
      mutate(
        Parameter = model_name
      )
  })
  
  marginals_df <- bind_rows(marginals_list)
  
  ggplot(marginals_df, aes(x, y, group = Parameter, color = Parameter)) +
    geom_line() +
    ggtitle("Posterior marginal distribution of sigma") +
    my_theme() +
    sigma_dist_labs()
}

# Standardies raster percentiles ###############################################

# Function to standardise richness values scale for each map using 1% brackets
standardise_raster_percentiles <- function(raster_data) {
  # Create a copy of the raster
  result <- raster_data
  
  # Get the original values and remove NAs for quantile calculation
  values <- terra::values(raster_data)
  non_na_values <- values[!is.na(values)]
  
  # Calculate quantile breaks for 1% intervals (100 brackets)
  breaks <- quantile(non_na_values, probs = seq(0, 1, by = 0.01))
  
  # Create a new vector to hold the quantized values
  new_values <- values
  
  # Assign percentile categories (1-100)
  for (i in 1:100) {
    if (i == 1) {
      # First percentile includes the minimum value
      mask <- values <= breaks[2] & !is.na(values)
    } else if (i == 100) {
      # Last percentile includes the maximum value
      mask <- values > breaks[100] & !is.na(values)
    } else {
      # Middle percentiles
      mask <- values > breaks[i] & values <= breaks[i+1] & !is.na(values)
    }
    
    # Assign the percentile number to matching cells
    new_values[mask] <- i
  }
  
  # Update the raster with the new values
  terra::values(result) <- new_values
  return(result)
}