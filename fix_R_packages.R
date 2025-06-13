#!/usr/bin/env Rscript

# Automated R Package Installer for Viral Mutation Visualizer
# Run this script ONCE to install all required packages
# Simply double-click this file or run: Rscript fix_R_packages.R

cat("==========================================================\n")
cat("Installing R Packages for Viral Mutation Visualizer\n")
cat("This may take 5-10 minutes. Please wait...\n")
cat("==========================================================\n\n")

# Function to safely install packages with better error handling
safe_install <- function(packages, repo_type = "CRAN") {
  for (pkg in packages) {
    cat(paste("Checking package:", pkg, "..."))
    
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(" INSTALLING\n")
      tryCatch({
        if (repo_type == "CRAN") {
          install.packages(pkg, repos = "https://cloud.r-project.org/", 
                         dependencies = TRUE, quiet = FALSE)
        } else if (repo_type == "Bioconductor") {
          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", repos = "https://cloud.r-project.org/")
          }
          BiocManager::install(pkg, quiet = FALSE)
        }
        cat(paste("   ✓ Successfully installed", pkg, "\n"))
      }, error = function(e) {
        cat(paste("   ✗ Failed to install", pkg, ":", e$message, "\n"))
      })
    } else {
      cat(" ALREADY INSTALLED ✓\n")
    }
  }
}

# Install CRAN packages (removing cowplot as it's not actually used)
cat("Installing CRAN packages...\n")
cran_packages <- c(
  "optparse",      # Command line parsing
  "ggplot2",       # Plotting
  "dplyr",         # Data manipulation
  "tidyr",         # Data tidying
  "readr",         # Reading data files
  "stringr",       # String manipulation
  "patchwork",     # Combining plots
  "ggrepel",       # Text labels that don't overlap
  "RColorBrewer",  # Color palettes
  "gridExtra",     # Grid graphics
  "grid"           # Base grid graphics
)

safe_install(cran_packages, "CRAN")

# Install Bioconductor packages
cat("\nInstalling Bioconductor packages...\n")
bioc_packages <- c(
  "rentrez",       # NCBI data retrieval
  "seqinr"         # Sequence analysis
)

safe_install(bioc_packages, "Bioconductor")

# Test that everything works
cat("\n==========================================================\n")
cat("Testing package installation...\n")
cat("==========================================================\n")

all_packages <- c(cran_packages, bioc_packages)
failed_packages <- c()

for (pkg in all_packages) {
  tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    cat(paste("✓", pkg, "loaded successfully\n"))
  }, error = function(e) {
    cat(paste("✗", pkg, "failed to load\n"))
    failed_packages <<- c(failed_packages, pkg)
  })
}

# Final report
cat("\n==========================================================\n")
cat("INSTALLATION COMPLETE\n")
cat("==========================================================\n")

if (length(failed_packages) == 0) {
  cat("🎉 SUCCESS! All packages installed correctly.\n")
  cat("You can now run the viral mutation visualization scripts.\n\n")
  cat("To visualize mutations, run:\n")
  cat("Rscript visualize_mutations.R --input your_data.tsv --output plot.pdf\n")
} else {
  cat("⚠️  Some packages failed to install:\n")
  for (pkg in failed_packages) {
    cat(paste("   -", pkg, "\n"))
  }
  cat("\nPlease contact your bioinformatics support for help.\n")
}

cat("\n==========================================================\n")