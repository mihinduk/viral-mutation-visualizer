#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# visualize_mutations_v5.R - Non-synonymous mutations only with separate output table
# -------------------------------------------------------------------------

# Function to install and load required packages
install_and_load_packages <- function() {
  # Define required packages
  cran_packages <- c("optparse", "ggplot2", "dplyr", "tidyr", "readr", "stringr", 
                     "patchwork", "ggrepel", "RColorBrewer", "gridExtra", "grid", "cowplot")
  bioc_packages <- c("rentrez", "seqinr")
  
  # Install missing CRAN packages
  missing_cran <- cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_cran) > 0) {
    cat("Installing missing CRAN packages:", paste(missing_cran, collapse = ", "), "\n")
    install.packages(missing_cran, repos = "https://cloud.r-project.org/", quiet = TRUE)
  }
  
  # Install BiocManager if needed for Bioconductor packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager", repos = "https://cloud.r-project.org/", quiet = TRUE)
  }
  
  # Install missing Bioconductor packages
  missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    cat("Installing missing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
    BiocManager::install(missing_bioc, quiet = TRUE)
  }
  
  # Load all packages
  all_packages <- c(cran_packages, bioc_packages)
  suppressPackageStartupMessages({
    for (pkg in all_packages) {
      library(pkg, character.only = TRUE)
    }
  })
  
  cat("All required packages loaded successfully.\n")
}

# Install and load required packages
install_and_load_packages()

# Parse command line arguments
option_list <- list(
  make_option(c("--input"), type="character", default=NULL, 
              help="Input TSV file with SnpEff-annotated LoFreq output", metavar="FILE"),
  make_option(c("--output"), type="character", default=NULL, 
              help="Output file path (PDF or PNG)", metavar="FILE"),
  make_option(c("--cutoff"), type="numeric", default=0.05, 
              help="Allele frequency cutoff [default= %default]", metavar="NUMBER"),
  make_option(c("--mutation-genes"), type="character", default="all", 
              help="Genes to display mutations for: 'all', 'structural', 'non-structural', or specific gene names (comma-separated)"),
  make_option(c("--colors"), type="character", default=NULL, 
              help="Custom color scheme, comma-separated hex values [optional]"),
  make_option(c("--width"), type="numeric", default=16, 
              help="Plot width in inches [default= %default]"),
  make_option(c("--height"), type="numeric", default=10, 
              help="Plot height in inches [default= %default]"),
  make_option(c("--accession"), type="character", default=NULL, 
              help="Manually specify GenBank accession (overrides accession in CHROM column)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$input)) {
  stop("Input file path is required. Use --input to specify the path.")
}
if (is.null(opt$output)) {
  stop("Output file path is required. Use --output to specify the path.")
}

# Function to read SnpEff-annotated LoFreq file
read_mutations <- function(file_path) {
  # First check if file has header
  first_line <- readLines(file_path, n = 1)
  has_header <- grepl("CHROM", first_line, ignore.case = TRUE)
  
  # Define column names based on SnpEff-annotated LoFreq output
  if (has_header) {
    mutations <- read_tsv(file_path, col_types = cols(.default = col_character()))
  } else {
    mutations <- read_tsv(file_path, col_names = c(
      "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", 
      "Total_Depth", "Allele_Frequency", "strand_bias", "DP4",
      "EFFECT", "PUTATIVE_IMPACT", "GENE_NAME", "GENE_ID", "FEATURE_TYPE", "FEATURE_ID", "TRANSCRIPT_TYPE",
      "HGVSc", "HGVSp", "cDNA_POSITION_AND_LENGTH", "CDS_POSITION_AND_LENGTH", "PROTEIN_POSITION_AND_LENGTH", "ERROR"
    ), col_types = cols(.default = col_character()))
  }
  
  # Make sure all expected columns exist
  required_cols <- c("CHROM", "POS", "REF", "ALT", "GENE_NAME")
  for (col in required_cols) {
    if (!(col %in% colnames(mutations))) {
      stop(paste("Required column", col, "not found in input file"))
    }
  }
  
  # Handle Allele_Frequency
  if ("Allele_Frequency" %in% colnames(mutations)) {
    mutations$Allele_Frequency <- as.numeric(mutations$Allele_Frequency)
  } else if ("AF" %in% colnames(mutations)) {
    mutations$Allele_Frequency <- as.numeric(mutations$AF)
  } else {
    # Try to extract from INFO field
    if ("INFO" %in% colnames(mutations)) {
      mutations$Allele_Frequency <- as.numeric(sapply(mutations$INFO, function(info) {
        af_match <- regexpr("AF=([0-9.]+)", info)
        if (af_match > 0) {
          af_value <- regmatches(info, af_match)
          return(as.numeric(gsub("AF=", "", af_value)))
        } else {
          return(NA)
        }
      }))
    } else {
      stop("Could not find allele frequency information")
    }
  }
  
  # Convert other numeric columns
  mutations$POS <- as.numeric(mutations$POS)
  
  # Extract protein positions if available
  if ("HGVSp" %in% colnames(mutations)) {
    mutations$protein_pos <- as.numeric(str_extract(mutations$HGVSp, "(?<=p\\.[A-Za-z]{1,3})\\d+"))
  }
  
  return(mutations)
}

# Function to get genome features
get_genome_features <- function(accession) {
  cat("Fetching genome features for accession:", accession, "\n")
  
  # Try to fetch from NCBI
  genome_info <- tryCatch({
    # Fetch the GenBank record
    gb_record <- entrez_fetch(db = "nucleotide", id = accession, rettype = "gb", retmode = "text")
    
    # Parse the GenBank record
    lines <- strsplit(gb_record, "\n")[[1]]
    
    # Extract definition
    definition_line <- grep("^DEFINITION", lines, value = TRUE)[1]
    definition <- gsub("^DEFINITION\\s+", "", definition_line)
    definition <- gsub("\\.$", "", definition)
    
    # Extract genome length from LOCUS line
    locus_line <- grep("^LOCUS", lines, value = TRUE)[1]
    genome_length <- as.numeric(gsub(".*\\s+(\\d+)\\s+bp.*", "\\1", locus_line))
    
    # Extract CDS and mat_peptide features
    cds_starts <- c()
    cds_ends <- c()
    gene_names <- c()
    products <- c()
    
    i <- 1
    while (i <= length(lines)) {
      if (grepl("^\\s+(CDS|mat_peptide)\\s+", lines[i])) {
        # Extract position - it's on the same line after the feature type
        pos_line <- lines[i]
        # The position is everything after CDS or mat_peptide
        pos_str <- trimws(gsub("^\\s+(CDS|mat_peptide)\\s+", "", pos_line))
        
        # Handle complement and join notations
        pos_str <- gsub("complement\\(|\\)|join\\(|\\)", "", pos_str)
        
        # Extract simple start..end positions
        if (grepl("\\d+\\.\\.\\d+", pos_str)) {
          # Extract numbers from position string
          matches <- gregexpr("\\d+", pos_str)
          positions <- as.numeric(regmatches(pos_str, matches)[[1]])
          if (length(positions) >= 2) {
            start_pos <- positions[1]
            end_pos <- positions[length(positions)]
            
            # Look for gene and product information
            gene_name <- NA
            product_name <- NA
            j <- i + 1
            while (j <= length(lines) && grepl("^\\s+/", lines[j])) {
              if (grepl("/gene=", lines[j])) {
                gene_name <- gsub('.*gene="(.*)".*', '\\1', lines[j])
              }
              if (grepl("/product=", lines[j])) {
                product_name <- gsub('.*product="(.*)".*', '\\1', lines[j])
                # Handle multi-line products
                k <- j + 1
                while (k <= length(lines) && !grepl("^\\s+/", lines[k]) && grepl('^\\s+"', lines[k])) {
                  product_name <- paste(product_name, gsub('\\s*"(.*)".*', '\\1', lines[k]))
                  k <- k + 1
                }
              }
              j <- j + 1
            }
            
            # For mat_peptide features without gene names, infer from product
            if (is.na(gene_name) && !is.na(product_name)) {
              # Map common product names to gene names
              if (grepl("capsid|core", product_name, ignore.case = TRUE)) gene_name <- "C"
              else if (grepl("premembrane|preM|prM", product_name, ignore.case = TRUE)) gene_name <- "prM"
              else if (grepl("envelope|^E protein", product_name, ignore.case = TRUE)) gene_name <- "Env"
              else if (grepl("NS1|nonstructural protein 1", product_name, ignore.case = TRUE)) gene_name <- "NS1"
              else if (grepl("NS2A|nonstructural protein 2A", product_name, ignore.case = TRUE)) gene_name <- "NS2a"
              else if (grepl("NS2B|nonstructural protein 2B", product_name, ignore.case = TRUE)) gene_name <- "NS2b"
              else if (grepl("NS3|nonstructural protein 3", product_name, ignore.case = TRUE)) gene_name <- "NS3"
              else if (grepl("NS4A|nonstructural protein 4A", product_name, ignore.case = TRUE)) gene_name <- "NS4a"
              else if (grepl("NS4B|nonstructural protein 4B", product_name, ignore.case = TRUE)) gene_name <- "NS4b"
              else if (grepl("NS5|nonstructural protein 5", product_name, ignore.case = TRUE)) gene_name <- "NS5"
            }
            
            # Store if we have a gene name or product
            if (!is.na(gene_name) || !is.na(product_name)) {
              cds_starts <- c(cds_starts, start_pos)
              cds_ends <- c(cds_ends, end_pos)
              gene_names <- c(gene_names, ifelse(is.na(gene_name), "", gene_name))
              products <- c(products, ifelse(is.na(product_name), "", product_name))
            }
          }
        }
      }
      i <- i + 1
    }
    
    # Process gene names based on virus type
    if (length(gene_names) > 0 || length(products) > 0) {
      # Check virus type based on products
      is_flavivirus <- any(grepl("NS[1-5]|nonstructural protein [1-5]", products, ignore.case = TRUE))
      is_alphavirus <- any(grepl("nsp[1-4]|E[1-3]|6K", products, ignore.case = TRUE))
      
      if (is_alphavirus) {
        # This is likely an alphavirus
        cat("Detected alphavirus genome structure\n")
        
        # Clean up gene names based on products if needed
        if (all(gene_names == "" | is.na(gene_names))) {
          # Map products to standard gene names for alphaviruses
          for (i in 1:length(products)) {
            if (grepl("nsp1", products[i], ignore.case = TRUE)) gene_names[i] <- "nsp1"
            else if (grepl("nsp2", products[i], ignore.case = TRUE)) gene_names[i] <- "nsp2"
            else if (grepl("nsp3", products[i], ignore.case = TRUE)) gene_names[i] <- "nsp3"
            else if (grepl("nsp4", products[i], ignore.case = TRUE)) gene_names[i] <- "nsp4"
            else if (grepl("^C$|capsid", products[i], ignore.case = TRUE)) gene_names[i] <- "C"
            else if (grepl("E3", products[i], ignore.case = TRUE)) gene_names[i] <- "E3"
            else if (grepl("E2", products[i], ignore.case = TRUE)) gene_names[i] <- "E2"
            else if (grepl("6K", products[i], ignore.case = TRUE)) gene_names[i] <- "6K"
            else if (grepl("E1", products[i], ignore.case = TRUE)) gene_names[i] <- "E1"
          }
        }
      } else if (is_flavivirus && length(cds_starts) >= 10) {
        # This is likely a flavivirus with mat_peptides
        cat("Detected flavivirus genome structure\n")
        
        # Clean up gene names based on products if needed
        if (all(gene_names == "" | is.na(gene_names))) {
          # Map products to standard gene names
          for (i in 1:length(products)) {
            if (grepl("capsid", products[i], ignore.case = TRUE)) gene_names[i] <- "C"
            else if (grepl("membrane-associated glycoprotein precursor|preM", products[i], ignore.case = TRUE)) gene_names[i] <- "prM"
            else if (grepl("envelope", products[i], ignore.case = TRUE)) gene_names[i] <- "Env"
            else if (grepl("NS1", products[i], ignore.case = TRUE)) gene_names[i] <- "NS1"
            else if (grepl("NS2a", products[i], ignore.case = TRUE)) gene_names[i] <- "NS2a"
            else if (grepl("NS2b", products[i], ignore.case = TRUE)) gene_names[i] <- "NS2b"
            else if (grepl("NS3", products[i], ignore.case = TRUE)) gene_names[i] <- "NS3"
            else if (grepl("NS4a", products[i], ignore.case = TRUE)) gene_names[i] <- "NS4a"
            else if (grepl("NS4b", products[i], ignore.case = TRUE)) gene_names[i] <- "NS4b"
            else if (grepl("NS5", products[i], ignore.case = TRUE)) gene_names[i] <- "NS5"
          }
        }
      } else if (is_flavivirus && length(cds_starts) == 1) {
        # Single polyprotein without mat_peptides - need to divide it
        polyprotein_start <- cds_starts[1]
        polyprotein_end <- cds_ends[1]
        polyprotein_length <- polyprotein_end - polyprotein_start + 1
        
        # Use approximate proportions based on typical flavivirus structure
        proportions <- c(0.037, 0.051, 0.146, 0.098, 0.063, 0.037, 0.175, 0.043, 0.072, 0.278)
        
        cds_starts <- numeric(10)
        cds_ends <- numeric(10)
        gene_names <- c("C", "prM", "Env", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5")
        current_pos <- polyprotein_start
        
        for (i in 1:10) {
          cds_starts[i] <- current_pos
          gene_length <- round(polyprotein_length * proportions[i])
          cds_ends[i] <- current_pos + gene_length - 1
          current_pos <- cds_ends[i] + 1
        }
        # Adjust last gene to end at polyprotein end
        cds_ends[10] <- polyprotein_end
        
        products <- c("capsid protein", "premembrane protein", "envelope protein E",
                     "nonstructural protein 1", "nonstructural protein 2A", "nonstructural protein 2B",
                     "nonstructural protein 3", "nonstructural protein 4A", "nonstructural protein 4B",
                     "nonstructural protein 5")
      }
      
      # Filter out empty entries and create genes data frame
      valid_indices <- which(gene_names != "" & !is.na(gene_names))
      if (length(valid_indices) > 0) {
        genes_data <- data.frame(
          gene = gene_names[valid_indices],
          start = cds_starts[valid_indices],
          end = cds_ends[valid_indices],
          product = products[valid_indices],
          strand = rep("+", length(valid_indices)),
          stringsAsFactors = FALSE
        )
        
        # Add category based on virus type
        if (any(genes_data$gene %in% c("nsp1", "nsp2", "nsp3", "nsp4"))) {
          # Alphavirus categories
          genes_data$category <- ifelse(genes_data$gene %in% c("C", "E3", "E2", "6K", "E1"), "structural", "non-structural")
        } else {
          # Flavivirus categories
          genes_data$category <- ifelse(genes_data$gene %in% c("C", "prM", "Env"), "structural", "non-structural")
        }
      } else {
        # If no valid gene names, return NULL to trigger default
        return(NULL)
      }
      
      return(list(
        genes = genes_data,
        genome_length = genome_length,
        definition = definition,
        accession = accession
      ))
    }
    
    # If no genes found, return NULL to trigger default
    return(NULL)
    
  }, error = function(e) {
    cat("Error fetching from NCBI:", e$message, "\n")
    return(NULL)
  })
  
  # If NCBI fetch failed or no genes found, use defaults
  if (is.null(genome_info)) {
    cat("Using default WNV genome structure\n")
    # Default WNV gene structure
    genes_data <- data.frame(
      gene = c("C", "prM", "Env", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5"),
      start = c(97, 466, 967, 2470, 3526, 4219, 4612, 6469, 6916, 7672),
      end = c(465, 966, 2469, 3525, 4218, 4611, 6468, 6915, 7671, 10395),
      product = c("capsid protein", "premembrane protein", "envelope protein E", 
                  "nonstructural protein 1", "nonstructural protein 2A", "nonstructural protein 2B",
                  "nonstructural protein 3", "nonstructural protein 4A", "nonstructural protein 4B", 
                  "nonstructural protein 5"),
      strand = rep("+", 10),
      category = c(rep("structural", 3), rep("non-structural", 7)),
      stringsAsFactors = FALSE
    )
    
    genome_length <- 11000
    definition <- "West Nile virus, complete genome (default)"
    
    genome_info <- list(
      genes = genes_data,
      genome_length = genome_length,
      definition = definition,
      accession = accession
    )
  }
  
  # Calculate UTR regions
  if (!is.null(genome_info$genes) && nrow(genome_info$genes) > 0) {
    first_cds_start <- min(genome_info$genes$start)
    last_cds_end <- max(genome_info$genes$end)
    
    genome_info$utrs <- data.frame(
      region = c("5'UTR", "3'UTR"),
      start = c(1, last_cds_end + 1),
      end = c(first_cds_start - 1, genome_info$genome_length),
      stringsAsFactors = FALSE
    )
  } else {
    # Create empty UTR data frame
    genome_info$utrs <- data.frame(
      region = character(0),
      start = numeric(0),
      end = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  return(genome_info)
}

# Function to format amino acid change
format_aa_change <- function(hgvsp) {
  # Vectorized function to handle multiple values
  sapply(hgvsp, function(x) {
    if (is.na(x) || x == "") return("")
    
    # Remove "p." prefix
    aa_change <- gsub("^p\\.", "", x)
    
    # Extract position number if present
    position_match <- regexpr("\\d+", aa_change)
    if (position_match > 0) {
      position <- regmatches(aa_change, position_match)
      # Get the part before and after the position
      before_pos <- substr(aa_change, 1, position_match - 1)
      after_pos <- substr(aa_change, position_match + attr(position_match, "match.length"), nchar(aa_change))
      
      # If we have amino acids before position, add space
      if (nchar(before_pos) > 0) {
        aa_change <- paste0(before_pos, " ", position, after_pos)
      }
    }
    
    # Keep full 3-letter codes
    return(aa_change)
  })
}

# Function to create per-gene mutation tables (NON-SYNONYMOUS ONLY)
create_gene_mutation_tables <- function(mutations, gene_colors) {
  # Filter for non-synonymous mutations only - INCLUDING stop mutations
  nonsyn_mutations <- mutations %>%
    filter(
      # Include any mutation with these effects
      grepl("missense|nonsense|stop_gained|stop_lost", EFFECT, ignore.case = TRUE) |
      # Include mutations with protein changes that aren't synonymous
      (!is.na(HGVSp) & HGVSp != "" & !grepl("=$", HGVSp)) |
      # Include mutations marked as HIGH impact
      grepl("HIGH", PUTATIVE_IMPACT, ignore.case = TRUE)
    ) %>%
    # Exclude only pure synonymous variants (not stop_gained)
    filter(!grepl("^synonymous_variant$", EFFECT, ignore.case = TRUE))
  
  # Get unique genes with mutations
  genes_with_mutations <- unique(nonsyn_mutations$GENE_NAME)
  genes_with_mutations <- genes_with_mutations[genes_with_mutations %in% names(gene_colors)]
  
  # Create table for each gene
  gene_tables <- list()
  
  for (gene in genes_with_mutations) {
    gene_muts <- nonsyn_mutations %>%
      filter(GENE_NAME == gene) %>%
      arrange(POS) %>%
      mutate(
        nt_change = paste0(POS, " ", REF, ">", ALT),
        aa_change = format_aa_change(HGVSp),
        is_stop = grepl("\\*", aa_change) | grepl("stop_gained|nonsense", EFFECT, ignore.case = TRUE)
      ) %>%
      # Remove any remaining synonymous mutations (where AA doesn't change)
      filter(!grepl("=$", aa_change))
    
    # Create formatted entries
    table_entries <- character()
    for (i in 1:nrow(gene_muts)) {
      entry <- gene_muts$nt_change[i]
      if (gene_muts$aa_change[i] != "") {
        if (gene_muts$is_stop[i]) {
          # Replace asterisk with TER for terminator in figure
          aa_display <- gsub("\\*", "TER", gene_muts$aa_change[i])
          entry <- paste0(entry, "\t", aa_display)
        } else {
          entry <- paste0(entry, "\t", gene_muts$aa_change[i])
        }
      }
      table_entries <- c(table_entries, entry)
    }
    
    # Only create table if there are mutations
    if (length(table_entries) > 0) {
      # Create data frame for table
      table_df <- data.frame(
        Mutations = table_entries,
        stringsAsFactors = FALSE
      )
      
      # Store with gene name as title
      gene_tables[[gene]] <- list(
        title = gene,
        data = table_df,
        color = gene_colors[gene]
      )
    }
  }
  
  return(gene_tables)
}

# Function to write all mutations to a separate table file
write_mutation_table <- function(mutations, output_file, gene_selection) {
  # Create a comprehensive mutation table
  mutation_table <- mutations %>%
    arrange(POS) %>%
    mutate(
      Nucleotide_Change = paste0(POS, " ", REF, ">", ALT),
      Amino_Acid_Change = format_aa_change(HGVSp),
      Mutation_Type = case_when(
        grepl("synonymous_variant", EFFECT, ignore.case = TRUE) & 
          !grepl("non_synonymous", EFFECT, ignore.case = TRUE) ~ "Synonymous",
        grepl("missense|nonsense|stop_gained", EFFECT, ignore.case = TRUE) ~ "Non-synonymous",
        grepl("upstream|downstream|intergenic", EFFECT, ignore.case = TRUE) ~ "Non-coding",
        TRUE ~ "Other"
      ),
      Frequency_Percent = round(Allele_Frequency * 100, 2)
    ) %>%
    select(
      Gene = GENE_NAME,
      Position = POS,
      Nucleotide_Change,
      Amino_Acid_Change,
      Mutation_Type,
      Frequency_Percent,
      Effect = EFFECT
    )
  
  # Create output filename
  base_name <- tools::file_path_sans_ext(output_file)
  table_file <- paste0(base_name, "_mutations_table.tsv")
  
  # Write the table
  write_tsv(mutation_table, table_file)
  cat("Mutation table written to:", table_file, "\n")
  
  # Also create a summary
  summary <- mutation_table %>%
    group_by(Gene, Mutation_Type) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = Mutation_Type, values_from = Count, values_fill = 0)
  
  summary_file <- paste0(base_name, "_mutations_summary.tsv")
  write_tsv(summary, summary_file)
  cat("Summary table written to:", summary_file, "\n")
}

# Function to create genome visualization with multiple tables
create_genome_visualization <- function(mutations, genome_features, cutoff, gene_selection, custom_colors = NULL) {
  # Filter mutations by allele frequency
  filtered_mutations <- mutations %>% filter(Allele_Frequency >= cutoff)
  
  # Map gene names to standard names and handle polyprotein
  if (nrow(filtered_mutations) == 0) {
    cat("No mutations to process after filtering\n")
    return(filtered_mutations)
  }
  
  for (i in 1:nrow(filtered_mutations)) {
    # Debug: Check what we're dealing with
    if (!("GENE_NAME" %in% names(filtered_mutations))) {
      cat("WARNING: GENE_NAME column not found in filtered_mutations\n")
      filtered_mutations$GENE_NAME <- NA
    }
    
    # Handle NULL, empty, or missing gene names
    if (is.null(filtered_mutations$GENE_NAME) || 
        length(filtered_mutations$GENE_NAME) < i ||
        is.null(filtered_mutations$GENE_NAME[i]) || 
        length(filtered_mutations$GENE_NAME[i]) == 0 ||
        filtered_mutations$GENE_NAME[i] == "") {
      filtered_mutations$GENE_NAME[i] <- "Intergenic"
      next
    }
    
    gene_name <- tolower(as.character(filtered_mutations$GENE_NAME[i]))
    
    if (is.na(gene_name)) {
      filtered_mutations$GENE_NAME[i] <- "Intergenic"
      next
    }
    
    # Check if it's polyprotein and map by position
    if (grepl("polyprotein", gene_name)) {
      pos <- filtered_mutations$POS[i]
      mapped <- FALSE
      
      # Find which gene this position falls in
      if (!is.null(genome_features$genes) && nrow(genome_features$genes) > 0) {
        for (j in 1:nrow(genome_features$genes)) {
          if (pos >= genome_features$genes$start[j] && pos <= genome_features$genes$end[j]) {
            filtered_mutations$GENE_NAME[i] <- genome_features$genes$gene[j]
            mapped <- TRUE
            break
          }
        }
      }
      
      if (!mapped) {
        # Check if in UTR regions
        if (!is.null(genome_features$genes) && nrow(genome_features$genes) > 0) {
          if (pos < min(genome_features$genes$start)) {
            filtered_mutations$GENE_NAME[i] <- "5'UTR"
          } else if (pos > max(genome_features$genes$end)) {
            filtered_mutations$GENE_NAME[i] <- "3'UTR"
          } else {
            filtered_mutations$GENE_NAME[i] <- "Intergenic"
          }
        } else {
          filtered_mutations$GENE_NAME[i] <- "Intergenic"
        }
      }
    } else {
      # Try to match specific gene names
      if (grepl("capsid|^c$", gene_name)) filtered_mutations$GENE_NAME[i] <- "C"
      else if (grepl("prm|premembrane", gene_name)) filtered_mutations$GENE_NAME[i] <- "prM"
      else if (grepl("env|envelope", gene_name)) filtered_mutations$GENE_NAME[i] <- "Env"
      else if (grepl("ns1", gene_name)) filtered_mutations$GENE_NAME[i] <- "NS1"
      else if (grepl("ns2a", gene_name)) filtered_mutations$GENE_NAME[i] <- "NS2a"
      else if (grepl("ns2b", gene_name)) filtered_mutations$GENE_NAME[i] <- "NS2b"
      else if (grepl("ns3", gene_name)) filtered_mutations$GENE_NAME[i] <- "NS3"
      else if (grepl("ns4a", gene_name)) filtered_mutations$GENE_NAME[i] <- "NS4a"
      else if (grepl("ns4b", gene_name)) filtered_mutations$GENE_NAME[i] <- "NS4b"
      else if (grepl("ns5", gene_name)) filtered_mutations$GENE_NAME[i] <- "NS5"
    }
  }
  
  # Filter by gene selection
  if (gene_selection != "all") {
    if (gene_selection == "structural") {
      structural_genes <- genome_features$genes$gene[genome_features$genes$category == "structural"]
      filtered_mutations <- filtered_mutations %>% filter(GENE_NAME %in% structural_genes)
    } else if (gene_selection == "non-structural") {
      non_structural_genes <- genome_features$genes$gene[genome_features$genes$category == "non-structural"]
      filtered_mutations <- filtered_mutations %>% filter(GENE_NAME %in% non_structural_genes)
    } else {
      selected_genes <- unlist(strsplit(gene_selection, ","))
      filtered_mutations <- filtered_mutations %>% filter(GENE_NAME %in% selected_genes)
    }
  }
  
  cat("Found", nrow(filtered_mutations), "mutations after filtering\n")
  
  # Define colors based on actual genes in the genome
  unique_genes <- unique(genome_features$genes$gene)
  
  # Check if this is an alphavirus or flavivirus
  if (any(c("nsp1", "nsp2", "nsp3", "nsp4") %in% unique_genes)) {
    # Alphavirus gene order
    gene_names <- c("nsp1", "nsp2", "nsp3", "nsp4", "C", "E3", "E2", "6K", "E1")
  } else {
    # Flavivirus gene order
    gene_names <- c("C", "prM", "Env", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5")
  }
  
  # Keep only genes that exist in this genome
  gene_names <- gene_names[gene_names %in% unique_genes]
  
  if (!is.null(custom_colors)) {
    color_values <- unlist(strsplit(custom_colors, ","))
    if (length(color_values) < length(gene_names)) {
      color_values <- brewer.pal(max(3, length(gene_names)), "Paired")
    }
  } else {
    # Use a color palette that works for any number of genes
    if (length(gene_names) <= 12) {
      color_values <- brewer.pal(max(3, length(gene_names)), "Paired")
    } else {
      color_values <- rainbow(length(gene_names))
    }
  }
  
  gene_colors <- setNames(color_values[1:length(gene_names)], gene_names)
  
  # Create the main plot area
  grid.newpage()
  
  # Define layout with reduced space for title
  pushViewport(viewport(layout = grid.layout(
    nrow = 3, 
    ncol = 1,
    heights = unit(c(0.3, 1.5, 2), "null")  # Reduced title space, more space for tables
  )))
  
  # Create genome diagram in center
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  
  genome_line_y <- 0.5
  
  genes_for_plot <- genome_features$genes %>%
    mutate(
      y = genome_line_y,
      x_center = (start + end) / 2
    )
  
  # Calculate structural and non-structural ranges
  str_genes <- genome_features$genes %>% filter(category == "structural")
  non_str_genes <- genome_features$genes %>% filter(category == "non-structural")
  
  str_range <- if(nrow(str_genes) > 0) c(min(str_genes$start), max(str_genes$end)) else NULL
  non_str_range <- if(nrow(non_str_genes) > 0) c(min(non_str_genes$start), max(non_str_genes$end)) else NULL
  
  # Create the genome plot
  p_genome <- ggplot() +
    # Main genome line
    geom_segment(aes(x = 1, xend = genome_features$genome_length, 
                     y = genome_line_y, yend = genome_line_y),
                 linewidth = 3, color = "black") +
    
    # Gene rectangles
    geom_rect(data = genes_for_plot,
              aes(xmin = start, xmax = end, 
                  ymin = y - 0.06, ymax = y + 0.06, fill = gene),
              color = "white", linewidth = 0.5) +
    
    # Gene labels directly on genes
    geom_text(data = genes_for_plot,
              aes(x = x_center, y = y, label = gene),
              size = 4, fontface = "bold", color = "white") +
    
    scale_x_continuous(name = "Genome position (nt)", 
                       breaks = seq(0, ceiling(genome_features$genome_length/1000)*1000, by = 1000),
                       labels = function(x) format(x, big.mark = ",")) +
    scale_y_continuous(name = NULL, limits = c(0.2, 0.8)) +
    scale_fill_manual(values = gene_colors) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12)
    )
  
  # Add UTR regions if they exist
  if (!is.null(genome_features$utrs) && nrow(genome_features$utrs) > 0) {
    p_genome <- p_genome +
      geom_segment(data = genome_features$utrs,
                   aes(x = start, xend = end, y = genome_line_y, yend = genome_line_y),
                   linewidth = 5, color = "grey60") +
      geom_text(data = genome_features$utrs, 
                aes(x = (start + end)/2, y = genome_line_y + 0.1, label = region), 
                size = 3)
  }
  
  # Add structural/non-structural labels if applicable
  if (!is.null(str_range)) {
    p_genome <- p_genome +
      geom_segment(aes(x = str_range[1], xend = str_range[2], 
                       y = genome_line_y - 0.15, yend = genome_line_y - 0.15),
                   color = "#2166ac", linewidth = 3, alpha = 0.5) +
      geom_text(aes(x = mean(str_range), y = genome_line_y - 0.22), 
                label = "Structural proteins", size = 4, fontface = "bold", color = "#2166ac")
  }
  
  if (!is.null(non_str_range)) {
    p_genome <- p_genome +
      geom_segment(aes(x = non_str_range[1], xend = non_str_range[2], 
                       y = genome_line_y - 0.15, yend = genome_line_y - 0.15),
                   color = "#762a83", linewidth = 3, alpha = 0.5) +
      geom_text(aes(x = mean(non_str_range), y = genome_line_y - 0.22), 
                label = "Non-structural proteins", size = 4, fontface = "bold", color = "#762a83")
  }
  
  # Add mutations to the plot
  if (nrow(filtered_mutations) > 0) {
    mutations_for_plot <- filtered_mutations %>%
      mutate(gene_color = sapply(GENE_NAME, function(g) {
        if (g %in% names(gene_colors)) gene_colors[g] else "#999999"
      }))
    
    # Add mutation markers with thinner lines
    p_genome <- p_genome +
      geom_segment(data = mutations_for_plot,
                   aes(x = POS, xend = POS, 
                       y = genome_line_y - 0.08, 
                       yend = genome_line_y - 0.12),
                   color = mutations_for_plot$gene_color, 
                   linewidth = 0.5, alpha = 0.8)
  }
  
  print(p_genome, newpage = FALSE)
  popViewport()
  
  # Add title at top
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text(paste0("Mutations in ", unique(mutations$CHROM)[1], " - ", 
                   genome_features$definition, " (cutoff: ", cutoff*100, "%)"),
            gp = gpar(fontsize = 14, fontface = "bold"))
  popViewport()
  
  # Create gene-specific mutation tables (NON-SYNONYMOUS ONLY)
  if (nrow(filtered_mutations) > 0) {
    gene_tables <- create_gene_mutation_tables(filtered_mutations, gene_colors)
    
    # Position all tables below the genome
    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
    
    # Get all gene tables in the correct order
    all_gene_names <- c("C", "prM", "Env", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5")
    ordered_tables <- gene_tables[names(gene_tables) %in% all_gene_names]
    ordered_tables <- ordered_tables[order(match(names(ordered_tables), all_gene_names))]
    
    if (length(ordered_tables) > 0) {
      # Create layout for all gene tables
      n_cols <- min(5, length(ordered_tables))  # Max 5 columns
      n_rows <- ceiling(length(ordered_tables) / n_cols)
      
      pushViewport(viewport(layout = grid.layout(
        nrow = n_rows,
        ncol = n_cols,
        widths = unit(rep(1, n_cols), "null"),
        heights = unit(rep(1, n_rows), "null")
      )))
      
      for (i in 1:length(ordered_tables)) {
        row <- ceiling(i / n_cols)
        col <- ((i - 1) %% n_cols) + 1
        
        pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
        
        table_data <- ordered_tables[[i]]
        
        # Create table theme
        # Check for TER in entries to make them bold
        fontfaces <- rep("plain", nrow(table_data$data))
        for (j in 1:nrow(table_data$data)) {
          if (grepl("TER", table_data$data$Mutations[j])) {
            fontfaces[j] <- "bold"
          }
        }
        
        table_theme <- ttheme_default(
          base_size = 8,
          core = list(
            fg_params = list(parse = FALSE, fontsize = 7, fontface = fontfaces),
            bg_params = list(fill = "white")
          ),
          colhead = list(
            fg_params = list(fontface = "bold", fontsize = 8),
            bg_params = list(fill = table_data$color, alpha = 0.7)
          )
        )
        
        # Add gene name as header
        table_grob <- tableGrob(
          table_data$data,
          rows = NULL,
          cols = table_data$title,
          theme = table_theme
        )
        
        grid.draw(table_grob)
        popViewport()
      }
      popViewport()
    }
    popViewport()
  }
  
  popViewport()  # Main layout
}

# Main execution
main <- function() {
  # Read mutations file
  cat("Reading mutations file:", opt$input, "\n")
  mutations <- read_mutations(opt$input)
  
  # Get accession number
  if (!is.null(opt$accession)) {
    accession <- opt$accession
  } else {
    accession <- unique(mutations$CHROM)[1]
    if (accession == "CHROM" || accession == "" || is.na(accession)) {
      accession <- "AY532665.1"  # Default WNV reference
    }
  }
  
  # Get genome features
  genome_features <- get_genome_features(accession)
  
  # Debug: Check what was returned
  if (!is.null(genome_features$genes)) {
    cat("Found", nrow(genome_features$genes), "genes in genome\n")
  } else {
    cat("WARNING: No genes found in genome features\n")
  }
  
  # Create visualization
  cat("Creating visualization...\n")
  
  # Save plot
  cat("Saving plot to:", opt$output, "\n")
  output_format <- tolower(tools::file_ext(opt$output))
  
  if (output_format == "pdf") {
    pdf(opt$output, width = opt$width, height = opt$height)
    create_genome_visualization(
      mutations, 
      genome_features, 
      opt$cutoff,
      opt$`mutation-genes`,
      opt$colors
    )
    dev.off()
  } else if (output_format %in% c("png", "jpg", "jpeg", "tiff")) {
    png(opt$output, width = opt$width, height = opt$height, units = "in", res = 300)
    create_genome_visualization(
      mutations, 
      genome_features, 
      opt$cutoff,
      opt$`mutation-genes`,
      opt$colors
    )
    dev.off()
  } else {
    cat("Unrecognized output format. Saving as PDF...\n")
    pdf(paste0(tools::file_path_sans_ext(opt$output), ".pdf"), 
        width = opt$width, height = opt$height)
    create_genome_visualization(
      mutations, 
      genome_features, 
      opt$cutoff,
      opt$`mutation-genes`,
      opt$colors
    )
    dev.off()
  }
  
  # Write mutation table
  filtered_mutations <- mutations %>% filter(Allele_Frequency >= opt$cutoff)
  
  # Map gene names before writing table
  for (i in 1:nrow(filtered_mutations)) {
    gene_name <- tolower(filtered_mutations$GENE_NAME[i])
    
    if (is.na(gene_name)) {
      filtered_mutations$GENE_NAME[i] <- "Intergenic"
      next
    }
    
    # Check if it's polyprotein and map by position
    if (grepl("polyprotein", gene_name)) {
      pos <- filtered_mutations$POS[i]
      mapped <- FALSE
      
      # Find which gene this position falls in
      if (!is.null(genome_features$genes) && nrow(genome_features$genes) > 0) {
        for (j in 1:nrow(genome_features$genes)) {
          if (pos >= genome_features$genes$start[j] && pos <= genome_features$genes$end[j]) {
            filtered_mutations$GENE_NAME[i] <- genome_features$genes$gene[j]
            mapped <- TRUE
            break
          }
        }
      }
      
      if (!mapped) {
        # Check if in UTR regions
        if (!is.null(genome_features$genes) && nrow(genome_features$genes) > 0) {
          if (pos < min(genome_features$genes$start)) {
            filtered_mutations$GENE_NAME[i] <- "5'UTR"
          } else if (pos > max(genome_features$genes$end)) {
            filtered_mutations$GENE_NAME[i] <- "3'UTR"
          } else {
            filtered_mutations$GENE_NAME[i] <- "Intergenic"
          }
        } else {
          filtered_mutations$GENE_NAME[i] <- "Intergenic"
        }
      }
    }
  }
  
  # Filter by gene selection for table output
  if (opt$`mutation-genes` != "all") {
    if (opt$`mutation-genes` == "structural") {
      structural_genes <- genome_features$genes$gene[genome_features$genes$category == "structural"]
      filtered_mutations <- filtered_mutations %>% filter(GENE_NAME %in% structural_genes)
    } else if (opt$`mutation-genes` == "non-structural") {
      non_structural_genes <- genome_features$genes$gene[genome_features$genes$category == "non-structural"]
      filtered_mutations <- filtered_mutations %>% filter(GENE_NAME %in% non_structural_genes)
    } else {
      selected_genes <- unlist(strsplit(opt$`mutation-genes`, ","))
      filtered_mutations <- filtered_mutations %>% filter(GENE_NAME %in% selected_genes)
    }
  }
  
  write_mutation_table(filtered_mutations, opt$output, opt$`mutation-genes`)
  
  cat("Done!\n")
}

# Run the main function
main()