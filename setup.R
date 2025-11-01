#!/usr/bin/env Rscript
# Setup script for Mycorrhiza-AusTemperateForests project
# Purpose: Install required R packages and check dependencies
# Author: Luke Florence
# Date: 2025

cat("===== Setting up Mycorrhiza-AusTemperateForests Project =====\n\n")

# List of required packages
required_packages <- c(
  "vegan",        # Community ecology analyses
  "dplyr",        # Data manipulation
  "tidyr",        # Data tidying
  "readr",        # Fast CSV reading
  "ggplot2",      # Data visualization
  "patchwork",    # Combining plots
  "MASS"          # Statistical functions (glm.nb)
)

# Optional packages (enhance functionality but not required)
optional_packages <- c(
  "phyloseq",     # Microbiome data handling
  "lme4",         # Mixed effects models
  "indicspecies", # Indicator species analysis
  "BiocManager"   # For installing Bioconductor packages
)

cat("Checking required packages...\n")

# Function to check and install packages
install_if_missing <- function(packages, is_optional = FALSE) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  
  if(length(new_packages)) {
    cat("\nInstalling", length(new_packages), "package(s):", 
        paste(new_packages, collapse = ", "), "\n")
    
    if (is_optional) {
      cat("(These are optional packages)\n")
      user_response <- readline(prompt = "Install? (y/n): ")
      if (tolower(user_response) != "y") {
        cat("Skipping optional packages\n")
        return()
      }
    }
    
    tryCatch({
      install.packages(new_packages, dependencies = TRUE, 
                      repos = "https://cran.rstudio.com/")
      cat("Successfully installed packages\n")
    }, error = function(e) {
      cat("Error installing packages:", e$message, "\n")
      cat("Please install manually using: install.packages(c('", 
          paste(new_packages, collapse = "', '"), "'))\n", sep = "")
    })
  } else {
    if (is_optional) {
      cat("All optional packages already installed\n")
    } else {
      cat("All required packages already installed\n")
    }
  }
}

# Install required packages
install_if_missing(required_packages, is_optional = FALSE)

cat("\n")

# Check optional packages
cat("Checking optional packages...\n")
install_if_missing(optional_packages, is_optional = TRUE)

cat("\n")

# Verify installation
cat("Verifying package installation...\n")

all_packages <- c(required_packages, optional_packages)
installation_status <- data.frame(
  Package = all_packages,
  Installed = all_packages %in% installed.packages()[,"Package"],
  Required = c(rep("Yes", length(required_packages)), 
               rep("No", length(optional_packages)))
)

print(installation_status)

cat("\n")

# Check R version
r_version <- R.version$version.string
cat("R version:", r_version, "\n")

if (as.numeric(R.version$major) < 4) {
  cat("WARNING: R version 4.0.0 or higher is recommended\n")
  cat("Current version may work but is not tested\n")
}

cat("\n")

# Create directory structure if it doesn't exist
cat("Checking project directory structure...\n")

dirs_to_create <- c(
  "data/raw",
  "data/processed",
  "scripts",
  "outputs/figures",
  "outputs/tables",
  "docs"
)

for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cat("  Created:", dir, "\n")
  } else {
    cat("  Exists:", dir, "\n")
  }
}

cat("\n")

# Summary
cat("===== Setup Complete =====\n\n")

failed_packages <- installation_status %>%
  filter(Required == "Yes" & Installed == FALSE)

if (nrow(failed_packages) > 0) {
  cat("WARNING: Some required packages failed to install:\n")
  print(failed_packages$Package)
  cat("\nPlease install these packages manually before running analyses\n")
} else {
  cat("All required packages are installed!\n")
  cat("You can now run the analysis scripts:\n")
  cat("  1. source('scripts/01_data_processing.R')\n")
  cat("  2. source('scripts/02_diversity_analysis.R')\n")
  cat("  3. source('scripts/03_nutrient_analysis.R')\n")
  cat("  4. source('scripts/04_visualization.R')\n\n")
  cat("Make sure to place your data files in data/raw/ first!\n")
}

cat("\nFor more information, see README.md and docs/methods.md\n")

# Save session info
cat("\nSaving session info...\n")
sink("setup_session_info.txt")
sessionInfo()
sink()
cat("Session info saved to setup_session_info.txt\n")
