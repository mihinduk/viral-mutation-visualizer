# Setup renv for package management
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# Initialize renv
renv::init()

# Install required packages
renv::install(c(
  "optparse", 
  "ggplot2", 
  "dplyr", 
  "readr", 
  "stringr", 
  "patchwork", 
  "ggrepel", 
  "RColorBrewer"
))

# For Bioconductor packages
renv::install("BiocManager")
BiocManager::install(c("rentrez", "seqinr"))

# Save the state of the project library
renv::snapshot()