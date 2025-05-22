#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# visualize_mutations.R - Visualize mutations from a SnpEff-annotated LoFreq TSV file
# 
# This script creates a genome visualization showing mapped mutations on the viral genome
# with support for filtering by gene regions and allele frequency.
# -------------------------------------------------------------------------

# Function to install and load required packages
install_and_load_packages <- function() {
  # Define required packages
  cran_packages <- c("optparse", "ggplot2", "dplyr", "readr", "stringr", 
                     "patchwork", "ggrepel", "RColorBrewer")
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
  make_option(c("--width"), type="numeric", default=14, 
              help="Plot width in inches [default= %default]"),
  make_option(c("--height"), type="numeric", default=8, 
              help="Plot height in inches [default= %default]"),
  make_option(c("--highlight-freq"), type="numeric", default=0.5, 
              help="Highlight mutations with frequency above this threshold [default= %default]"),
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
  
  # Handle Allele_Frequency (might be named differently)
  if ("Allele_Frequency" %in% colnames(mutations)) {
    mutations$Allele_Frequency <- as.numeric(mutations$Allele_Frequency)
  } else if ("AF" %in% colnames(mutations)) {
    mutations$Allele_Frequency <- as.numeric(mutations$AF)
  } else {
    # Try to extract from INFO field if it exists
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
      stop("Could not find allele frequency information in the input file")
    }
  }
  
  # Convert other numeric columns
  mutations$POS <- as.numeric(mutations$POS)
  if ("Total_Depth" %in% colnames(mutations)) {
    mutations$Total_Depth <- as.numeric(mutations$Total_Depth)
  } else if ("DP" %in% colnames(mutations)) {
    mutations$Total_Depth <- as.numeric(mutations$DP)
  }
  
  # Extract protein positions if available
  if ("PROTEIN_POSITION_AND_LENGTH" %in% colnames(mutations)) {
    mutations$protein_pos <- as.numeric(str_extract(mutations$PROTEIN_POSITION_AND_LENGTH, "^\\d+"))
  } else if ("HGVSp" %in% colnames(mutations)) {
    # Try to extract from HGVSp (format like p.Ala123Val)
    mutations$protein_pos <- as.numeric(str_extract(mutations$HGVSp, "(?<=p\\.[A-Za-z]{1,3})\\d+"))
  }
  
  return(mutations)
}

# Function to get GenBank file and extract gene positions
get_genome_features <- function(accession) {
  # Create a temp file for the GenBank file
  gb_file <- tempfile(fileext = ".gb")
  
  # Download GenBank file with retry mechanism
  cat("Fetching GenBank record for", accession, "...\n")
  
  # Try to fetch the record with retry mechanism
  max_attempts <- 3
  for (attempt in 1:max_attempts) {
    tryCatch({
      entrez_fetch(db="nucleotide", id=accession, rettype="gb", retmode="text", file=gb_file)
      
      # Verify file was created and has content
      if (file.exists(gb_file) && file.size(gb_file) > 0) {
        break  # Success
      } else {
        stop("Empty GenBank file downloaded")
      }
    }, error = function(e) {
      if (attempt < max_attempts) {
        cat("Attempt", attempt, "failed. Retrying...\n")
        Sys.sleep(2)  # Wait 2 seconds before retry
      } else {
        cat("Error fetching GenBank record:", e$message, "\n")
        cat("Using default WNV gene structure instead.\n")
        
        # Create default WNV gene structure
        write("LOCUS       DEFAULT_WNV          11000 bp    RNA     linear   VRL 01-JAN-2025
FEATURES             Location/Qualifiers
     source          1..11000
                     /organism=\"West Nile virus\"
     CDS             97..465
                     /gene=\"C\"
                     /product=\"capsid protein\"
     CDS             466..966
                     /gene=\"prM\"
                     /product=\"premembrane protein\"
     CDS             967..2469
                     /gene=\"Env\"
                     /product=\"envelope protein E\"
     CDS             2470..3525
                     /gene=\"NS1\"
                     /product=\"nonstructural protein 1\"
     CDS             3526..4218
                     /gene=\"NS2a\"
                     /product=\"nonstructural protein 2A\"
     CDS             4219..4611
                     /gene=\"NS2b\"
                     /product=\"nonstructural protein 2B\"
     CDS             4612..6468
                     /gene=\"NS3\"
                     /product=\"nonstructural protein 3\"
     CDS             6469..6915
                     /gene=\"NS4a\"
                     /product=\"nonstructural protein 4A\"
     CDS             6916..7671
                     /gene=\"NS4b\"
                     /product=\"nonstructural protein 4B\"
     CDS             7672..10395
                     /gene=\"NS5\"
                     /product=\"nonstructural protein 5\"
ORIGIN
        1 agtagttcgc ctgtgtgagc tgacaaactt agttagtgtt tgtgagctgc aaacttgcga
       61 gagtctcccg gaaatagtgg gtgtttatct agaaacaaca attaatgtgg aaggcctcgt
      121 cccctcggac cgcctgatcg ctgcgccacc tgctgtctct gcgaatcact gtgggaatcc
      181 agcccggtgg ctgctctaat ggtcaccacg aatgcaacat acgacgctct tgtaaaactg
      //", gb_file)
      }
    })
  }
  
  # Verify the file exists and has content
  if (!file.exists(gb_file) || file.size(gb_file) == 0) {
    cat("GenBank file not created properly. Using default WNV gene structure.\n")
    # Create default WNV gene structure
    write("LOCUS       DEFAULT_WNV          11000 bp    RNA     linear   VRL 01-JAN-2025
FEATURES             Location/Qualifiers
     source          1..11000
                     /organism=\"West Nile virus\"
     CDS             97..465
                     /gene=\"C\"
                     /product=\"capsid protein\"
     CDS             466..966
                     /gene=\"prM\"
                     /product=\"premembrane protein\"
     CDS             967..2469
                     /gene=\"Env\"
                     /product=\"envelope protein E\"
     CDS             2470..3525
                     /gene=\"NS1\"
                     /product=\"nonstructural protein 1\"
     CDS             3526..4218
                     /gene=\"NS2a\"
                     /product=\"nonstructural protein 2A\"
     CDS             4219..4611
                     /gene=\"NS2b\"
                     /product=\"nonstructural protein 2B\"
     CDS             4612..6468
                     /gene=\"NS3\"
                     /product=\"nonstructural protein 3\"
     CDS             6469..6915
                     /gene=\"NS4a\"
                     /product=\"nonstructural protein 4A\"
     CDS             6916..7671
                     /gene=\"NS4b\"
                     /product=\"nonstructural protein 4B\"
     CDS             7672..10395
                     /gene=\"NS5\"
                     /product=\"nonstructural protein 5\"
ORIGIN
        1 agtagttcgc ctgtgtgagc tgacaaactt agttagtgtt tgtgagctgc aaacttgcga
       61 gagtctcccg gaaatagtgg gtgtttatct agaaacaaca attaatgtgg aaggcctcgt
      121 cccctcggac cgcctgatcg ctgcgccacc tgctgtctct gcgaatcact gtgggaatcc
      181 agcccggtgg ctgctctaat ggtcaccacg aatgcaacat acgacgctct tgtaaaactg
      //", gb_file)
  }
  
  # Read the GenBank file directly
  gb_content <- readLines(gb_file)
  
  # Parse the GenBank file to get features
  features <- list()
  in_features <- FALSE
  current_feature <- NULL
  feature_type <- NULL
  feature_location <- NULL
  gene_name <- NULL
  product_name <- NULL
  genome_length <- NULL
  definition <- NULL
  
  for (line in gb_content) {
    if (grepl("^LOCUS", line)) {
      # Extract genome length from LOCUS line
      parts <- strsplit(line, "\\s+")[[1]]
      for (i in 1:length(parts)) {
        if (grepl("^[0-9]+$", parts[i])) {
          genome_length <- as.numeric(parts[i])
          break
        }
      }
    }
    
    if (grepl("^DEFINITION", line)) {
      # Extract definition line (organism description)
      definition <- gsub("^DEFINITION\\s+", "", line)
      # Handle continuation lines for definition
      next_line_idx <- which(gb_content == line) + 1
      if (next_line_idx <= length(gb_content)) {
        while (next_line_idx <= length(gb_content) && 
               grepl("^\\s{12}", gb_content[next_line_idx]) && 
               !grepl("^ACCESSION|^VERSION|^KEYWORDS", gb_content[next_line_idx])) {
          definition <- paste(definition, trimws(gb_content[next_line_idx]))
          next_line_idx <- next_line_idx + 1
        }
      }
    }
    
    if (grepl("^FEATURES", line)) {
      in_features <- TRUE
      next
    }
    
    if (in_features) {
      if (grepl("^\\s{5}[a-zA-Z]", line)) {
        # Start of a new feature
        if (!is.null(current_feature) && !is.null(feature_type) && feature_type == "CDS") {
          # Save previous CDS feature
          features[[length(features) + 1]] <- list(
            type = feature_type,
            location = feature_location,
            gene = ifelse(is.null(gene_name), paste0("gene_", length(features) + 1), gene_name),
            product = ifelse(is.null(product_name), ifelse(is.null(gene_name), paste0("gene_", length(features) + 1), gene_name), product_name)
          )
        }
        
        # Reset for new feature
        current_feature <- line
        feature_type <- trimws(strsplit(line, "\\s+")[[1]][1])
        feature_location <- gsub("^\\s*[a-zA-Z]+\\s+(.*)", "\\1", line)
        gene_name <- NULL
        product_name <- NULL
      } else if (grepl("^\\s{21}/gene=", line)) {
        # Gene name
        gene_name <- gsub("^\\s*\\/gene=\"(.*)\".*", "\\1", line)
      } else if (grepl("^\\s{21}/product=", line)) {
        # Product name
        product_name <- gsub("^\\s*\\/product=\"(.*)\".*", "\\1", line)
      } else if (grepl("^\\s{21}", line) && !is.null(current_feature)) {
        # Continuation of a feature or qualification
        if (grepl("^\\s{21}[0-9]", line) || grepl("^\\s{21}complement", line) || grepl("^\\s{21}join", line)) {
          # Continuation of location
          feature_location <- paste0(feature_location, gsub("^\\s*", "", line))
        } else if (grepl("\\/gene=", line)) {
          gene_name <- gsub("^\\s*\\/gene=\"(.*)\".*", "\\1", line)
        } else if (grepl("\\/product=", line)) {
          product_name <- gsub("^\\s*\\/product=\"(.*)\".*", "\\1", line)
        }
      }
    }
    
    if (grepl("^ORIGIN", line)) {
      in_features <- FALSE
      # Save the last feature if it's a CDS
      if (!is.null(current_feature) && !is.null(feature_type) && feature_type == "CDS") {
        features[[length(features) + 1]] <- list(
          type = feature_type,
          location = feature_location,
          gene = ifelse(is.null(gene_name), paste0("gene_", length(features) + 1), gene_name),
          product = ifelse(is.null(product_name), ifelse(is.null(gene_name), paste0("gene_", length(features) + 1), gene_name), product_name)
        )
      }
      break
    }
  }
  
  # Process features to extract gene coordinates
  genes_data <- data.frame(
    gene = character(),
    start = numeric(),
    end = numeric(),
    product = character(),
    strand = character(),
    stringsAsFactors = FALSE
  )
  
  if (length(features) > 0) {
    for (i in 1:length(features)) {
      if (features[[i]]$type == "CDS") {
        loc <- features[[i]]$location
        
        # Determine strand
        strand <- "+"
        if (grepl("complement", loc)) {
          strand <- "-"
          loc <- gsub("complement\\((.*)\\)", "\\1", loc)
        }
        
        # Extract coordinates
        if (grepl("join", loc)) {
          # Handle joined features
          loc <- gsub("join\\((.*)\\)", "\\1", loc)
          parts <- strsplit(loc, ",")[[1]]
          starts <- numeric()
          ends <- numeric()
          
          for (part in parts) {
            coords <- as.numeric(strsplit(gsub("[^0-9]", " ", part), "\\s+")[[1]])
            coords <- coords[coords > 0]
            if (length(coords) >= 2) {
              starts <- c(starts, coords[1])
              ends <- c(ends, coords[2])
            }
          }
          
          start_pos <- min(starts)
          end_pos <- max(ends)
        } else {
          # Handle simple features
          coords <- as.numeric(strsplit(gsub("[^0-9]", " ", loc), "\\s+")[[1]])
          coords <- coords[coords > 0]
          if (length(coords) >= 2) {
            start_pos <- coords[1]
            end_pos <- coords[2]
          } else {
            next
          }
        }
        
        genes_data <- rbind(genes_data, data.frame(
          gene = features[[i]]$gene,
          start = start_pos,
          end = end_pos,
          product = features[[i]]$product,
          strand = strand,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Convert to numeric
  genes_data$start <- as.numeric(genes_data$start)
  genes_data$end <- as.numeric(genes_data$end)
  
  # If no genes were found, create default WNV gene structure
  if (nrow(genes_data) == 0) {
    cat("No gene features found in GenBank file. Using default WNV gene structure.\n")
    genes_data <- data.frame(
      gene = c("C", "prM", "Env", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5"),
      start = c(97, 466, 967, 2470, 3526, 4219, 4612, 6469, 6916, 7672),
      end = c(465, 966, 2469, 3525, 4218, 4611, 6468, 6915, 7671, 10395),
      product = c("capsid protein", "premembrane protein", "envelope protein E", 
                  "nonstructural protein 1", "nonstructural protein 2A", "nonstructural protein 2B",
                  "nonstructural protein 3", "nonstructural protein 4A", "nonstructural protein 4B", 
                  "nonstructural protein 5"),
      strand = rep("+", 10),
      stringsAsFactors = FALSE
    )
    genome_length <- 11000
    definition <- "West Nile virus, complete genome"
  }
  
  # Create a lookup table for WNV genes to categorize them as structural or non-structural
  gene_categories <- data.frame(
    gene = c("C", "prM", "Env", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5"),
    category = c(rep("structural", 3), rep("non-structural", 7)),
    stringsAsFactors = FALSE
  )
  
  # Try to match genes with known WNV genes
  for (i in 1:nrow(genes_data)) {
    product <- tolower(genes_data$product[i])
    gene <- tolower(genes_data$gene[i])
    
    # Check for known genes in product or gene name
    if (grepl("capsid|core", product)) {
      genes_data$gene[i] <- "C"
    } else if (grepl("premembrane|pr[^a-z]?m", product) || grepl("^pr[^a-z]?m$", gene)) {
      genes_data$gene[i] <- "prM"
    } else if (grepl("envelope|env", product) || grepl("^env$", gene)) {
      genes_data$gene[i] <- "Env"
    } else if (grepl("ns1", product) || grepl("^ns1$", gene)) {
      genes_data$gene[i] <- "NS1"
    } else if (grepl("ns2a", product) || grepl("^ns2a$", gene)) {
      genes_data$gene[i] <- "NS2a"
    } else if (grepl("ns2b", product) || grepl("^ns2b$", gene)) {
      genes_data$gene[i] <- "NS2b"
    } else if (grepl("ns3", product) || grepl("^ns3$", gene)) {
      genes_data$gene[i] <- "NS3"
    } else if (grepl("ns4a", product) || grepl("^ns4a$", gene)) {
      genes_data$gene[i] <- "NS4a"
    } else if (grepl("ns4b", product) || grepl("^ns4b$", gene)) {
      genes_data$gene[i] <- "NS4b"
    } else if (grepl("ns5", product) || grepl("^ns5$", gene)) {
      genes_data$gene[i] <- "NS5"
    }
  }
  
  # Add category for each gene
  genes_data$category <- sapply(genes_data$gene, function(g) {
    idx <- which(gene_categories$gene == g)
    if (length(idx) > 0) {
      return(gene_categories$category[idx])
    } else {
      return("unknown")
    }
  })
  
  # Get UTR regions
  first_cds_start <- min(genes_data$start)
  last_cds_end <- max(genes_data$end)
  
  # If we don't have a genome length, estimate it as 120% of the end of the last CDS
  if (is.null(genome_length) || is.na(genome_length)) {
    genome_length <- ceiling(last_cds_end * 1.2)
    cat("Genome length not found, estimated as:", genome_length, "\n")
  }
  
  utr_data <- data.frame(
    region = c("5'UTR", "3'UTR"),
    start = c(1, last_cds_end + 1),
    end = c(first_cds_start - 1, genome_length),
    stringsAsFactors = FALSE
  )
  
  # Return the combined data
  return(list(
    genes = genes_data,
    utrs = utr_data,
    genome_length = genome_length,
    definition = ifelse(is.null(definition), "West Nile virus, complete genome", definition)
  ))
}

# Function to filter mutations by gene selection
filter_mutations_by_genes <- function(mutations, genome_features, gene_selection) {
  if (gene_selection == "all") {
    return(mutations)
  } else if (gene_selection == "structural") {
    structural_genes <- genome_features$genes$gene[genome_features$genes$category == "structural"]
    return(mutations %>% filter(GENE_NAME %in% structural_genes))
  } else if (gene_selection == "non-structural") {
    non_structural_genes <- genome_features$genes$gene[genome_features$genes$category == "non-structural"]
    return(mutations %>% filter(GENE_NAME %in% non_structural_genes))
  } else {
    # Specific gene selection (comma-separated list)
    selected_genes <- unlist(strsplit(gene_selection, ","))
    return(mutations %>% filter(GENE_NAME %in% selected_genes))
  }
}

# Function to create genome visualization
create_genome_visualization <- function(mutations, genome_features, cutoff, gene_selection, custom_colors = NULL, highlight_freq = 0.5) {
  # Filter mutations by allele frequency
  filtered_mutations <- mutations %>% filter(Allele_Frequency >= cutoff)
  
  # Set default gene_selection if it's empty
  if (is.null(gene_selection) || length(gene_selection) == 0) {
    gene_selection <- "all"
  }
  
  # Set default highlight_freq if it's NULL or invalid
  if (is.null(highlight_freq) || is.na(highlight_freq) || !is.numeric(highlight_freq)) {
    highlight_freq <- 0.5
  }
  
  # Map gene names to standard names if needed
  for (i in 1:nrow(filtered_mutations)) {
    # Skip NAs
    if (is.na(filtered_mutations$GENE_NAME[i])) next
    
    # Try to match gene names based on content
    gene_name <- filtered_mutations$GENE_NAME[i]
    
    if (grepl("capsid|^c$", gene_name, ignore.case = TRUE)) {
      filtered_mutations$GENE_NAME[i] <- "C"
    } else if (grepl("prm|prec?m|premembrane", gene_name, ignore.case = TRUE)) {
      filtered_mutations$GENE_NAME[i] <- "prM"
    } else if (grepl("^env$|envelope|^e$", gene_name, ignore.case = TRUE)) {
      filtered_mutations$GENE_NAME[i] <- "Env"
    } else if (grepl("ns1", gene_name, ignore.case = TRUE)) {
      filtered_mutations$GENE_NAME[i] <- "NS1"
    } else if (grepl("ns2a", gene_name, ignore.case = TRUE)) {
      filtered_mutations$GENE_NAME[i] <- "NS2a"
    } else if (grepl("ns2b", gene_name, ignore.case = TRUE)) {
      filtered_mutations$GENE_NAME[i] <- "NS2b"
    } else if (grepl("ns3", gene_name, ignore.case = TRUE)) {
      filtered_mutations$GENE_NAME[i] <- "NS3"
    } else if (grepl("ns4a", gene_name, ignore.case = TRUE)) {
      filtered_mutations$GENE_NAME[i] <- "NS4a"
    } else if (grepl("ns4b", gene_name, ignore.case = TRUE)) {
      filtered_mutations$GENE_NAME[i] <- "NS4b"
    } else if (grepl("ns5", gene_name, ignore.case = TRUE)) {
      filtered_mutations$GENE_NAME[i] <- "NS5"
    } else if (grepl("polyprotein", gene_name, ignore.case = TRUE)) {
      # For polyprotein, try to determine the gene based on position
      if ("POS" %in% names(filtered_mutations)) {
        pos <- filtered_mutations$POS[i]
        for (j in 1:nrow(genome_features$genes)) {
          if (pos >= genome_features$genes$start[j] && 
              pos <= genome_features$genes$end[j]) {
            filtered_mutations$GENE_NAME[i] <- genome_features$genes$gene[j]
            break
          }
        }
      }
    }
  }
  
  # Filter mutations by gene selection
  orig_count <- nrow(filtered_mutations)
  
  if (gene_selection != "all") {
    if (gene_selection == "structural") {
      structural_genes <- genome_features$genes$gene[genome_features$genes$category == "structural"]
      filtered_mutations <- filtered_mutations %>% filter(GENE_NAME %in% structural_genes)
    } else if (gene_selection == "non-structural") {
      non_structural_genes <- genome_features$genes$gene[genome_features$genes$category == "non-structural"]
      filtered_mutations <- filtered_mutations %>% filter(GENE_NAME %in% non_structural_genes)
    } else {
      # Specific gene selection (comma-separated list)
      selected_genes <- unlist(strsplit(gene_selection, ","))
      filtered_mutations <- filtered_mutations %>% filter(GENE_NAME %in% selected_genes)
    }
  }
  
  # Check if we have mutations after filtering
  if (nrow(filtered_mutations) == 0) {
    cat("No mutations met the criteria after filtering by gene and frequency.\n")
    cat("Original mutations count: ", orig_count, "\n")
    cat("Gene selection: ", gene_selection, "\n")
    # Print unique gene names to help debug
    if (orig_count > 0) {
      cat("Available gene names in data: ", paste(unique(mutations$GENE_NAME), collapse=", "), "\n")
    }
  } else {
    cat("After filtering, kept", nrow(filtered_mutations), "mutations out of", orig_count, "\n")
  }
  
  # Calculate structural and non-structural gene ranges
  str_genes <- genome_features$genes %>% filter(category == "structural")
  non_str_genes <- genome_features$genes %>% filter(category == "non-structural")
  
  str_range <- c(min(str_genes$start), max(str_genes$end))
  non_str_range <- c(min(non_str_genes$start), max(non_str_genes$end))
  
  # Define colors for genes
  gene_names <- c("C", "prM", "Env", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5")
  
  if (!is.null(custom_colors)) {
    color_values <- unlist(strsplit(custom_colors, ","))
    if (length(color_values) < length(gene_names)) {
      warning("Not enough custom colors provided. Using default colors.")
      color_values <- brewer.pal(length(gene_names), "Paired")
    }
  } else {
    color_values <- c("#4575b4", "#74add1", "#abd9e9", "#fdae61", "#f46d43", "#d73027", "#a50026", "#762a83", "#9970ab", "#c2a5cf")
  }
  
  gene_colors <- setNames(color_values, gene_names)
  
  # Prepare data for plotting
  # 1. Genome layout with gene regions
  genes_for_plot <- genome_features$genes %>%
    mutate(y = ifelse(category == "structural", 0.8, 0.2))
  
  # 2. High-level structural vs non-structural regions
  category_regions <- data.frame(
    category = c("Structural proteins", "Non-structural proteins"),
    start = c(str_range[1], non_str_range[1]),
    end = c(str_range[2], non_str_range[2]),
    y = c(0.8, 0.2),
    stringsAsFactors = FALSE
  )
  
  # 3. UTR regions
  utrs_for_plot <- genome_features$utrs %>%
    mutate(y = 0.5)
  
  # 4. Mutations
  if (nrow(filtered_mutations) > 0) {
    # First create a simplified data frame with essential columns
    mutations_for_plot <- data.frame(
      POS = filtered_mutations$POS,
      GENE_NAME = ifelse(is.na(filtered_mutations$GENE_NAME) | filtered_mutations$GENE_NAME == "", 
                         "Intergenic", filtered_mutations$GENE_NAME),
      Allele_Frequency = filtered_mutations$Allele_Frequency,
      gene_color = rep("gray50", nrow(filtered_mutations)),
      highlight = rep(FALSE, nrow(filtered_mutations)),
      label = rep(NA, nrow(filtered_mutations)),
      stringsAsFactors = FALSE
    )
    
    # Now add colors
    for (i in 1:nrow(mutations_for_plot)) {
      gene <- mutations_for_plot$GENE_NAME[i]
      if (gene %in% names(gene_colors)) {
        mutations_for_plot$gene_color[i] <- gene_colors[gene]
      } else {
        mutations_for_plot$gene_color[i] <- "#999999"  # Default gray for unrecognized genes
      }
    }
    
    # Add highlight flags
    mutations_for_plot$highlight <- mutations_for_plot$Allele_Frequency >= highlight_freq
    
    # Add labels
    for (i in 1:nrow(mutations_for_plot)) {
      # Only label highlighted mutations
      if (!mutations_for_plot$highlight[i]) next
      
      if (i <= nrow(filtered_mutations) && 
          !is.na(filtered_mutations$protein_pos[i]) && 
          !is.na(filtered_mutations$HGVSp[i]) && 
          filtered_mutations$HGVSp[i] != "") {
        mutations_for_plot$label[i] <- paste0(
          mutations_for_plot$GENE_NAME[i], ":", 
          filtered_mutations$HGVSp[i], " (", 
          round(mutations_for_plot$Allele_Frequency[i]*100), "%)"
        )
      } else {
        mutations_for_plot$label[i] <- paste0(
          mutations_for_plot$GENE_NAME[i], ":", 
          mutations_for_plot$POS[i], " (", 
          round(mutations_for_plot$Allele_Frequency[i]*100), "%)"
        )
      }
    }
  } else {
    # Create an empty data frame for plotting
    mutations_for_plot <- data.frame(
      POS = numeric(0),
      GENE_NAME = character(0),
      Allele_Frequency = numeric(0),
      gene_color = character(0),
      highlight = logical(0),
      label = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  # Create the genome layout plot
  # 1. Base layout with UTRs
  p_genome <- ggplot() +
    # Add UTR regions
    geom_rect(data = utrs_for_plot, 
              aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55, fill = region),
              color = "black", fill = "grey80") +
    geom_text(data = utrs_for_plot, 
              aes(x = (start + end)/2, y = 0.6, label = region), 
              size = 4) +
    
    # Add high-level category regions (structural vs non-structural)
    geom_rect(data = category_regions,
              aes(xmin = start, xmax = end, ymin = y - 0.12, ymax = y + 0.12, fill = category),
              color = "black", fill = NA) +
    geom_text(data = category_regions,
              aes(x = (start + end)/2, y = y + 0.22, label = category),
              size = 5, fontface = "bold") +
              
    # Create a background region showing the full genome length
    geom_rect(aes(xmin = 1, xmax = genome_features$genome_length, ymin = 0.45, ymax = 0.55),
              fill = NA, color = "black") +
              
    # Add gene regions
    geom_rect(data = genes_for_plot,
              aes(xmin = start, xmax = end, ymin = y - 0.08, ymax = y + 0.08, fill = gene),
              color = "black") +
    geom_text(data = genes_for_plot,
              aes(x = (start + end)/2, y = y, label = gene),
              size = 4, fontface = "bold") +
              
    # Add structural/non-structural protein lengths
    geom_text(data = data.frame(
      x = c((str_range[1] + str_range[2])/2, (non_str_range[1] + non_str_range[2])/2),
      y = c(0.8 - 0.18, 0.2 - 0.18),
      label = c(paste0("(", max(str_genes$end) - min(str_genes$start) + 1, "nt/", 
                      ceiling((max(str_genes$end) - min(str_genes$start) + 1)/3), "aa)"),
                paste0("(", max(non_str_genes$end) - min(non_str_genes$start) + 1, "nt/", 
                      ceiling((max(non_str_genes$end) - min(non_str_genes$start) + 1)/3), "aa)")
      )
    ), aes(x = x, y = y, label = label), size = 4) +
              
    # Add a scale bar at the bottom
    scale_x_continuous(name = "Genome position (nt)", 
                       breaks = seq(0, ceiling(genome_features$genome_length/1000)*1000, by = 1000),
                       labels = function(x) format(x, big.mark = ",")) +
    scale_y_continuous(name = NULL, limits = c(-2, 1.2)) +
    scale_fill_manual(values = c(gene_colors, "5'UTR" = "grey80", "3'UTR" = "grey80",
                                "Structural proteins" = NA, "Non-structural proteins" = NA)) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  # 2. Create mutation plot
  p_mutations <- ggplot() +
    # Add gene regions as background
    geom_rect(data = genes_for_plot,
              aes(xmin = start, xmax = end, ymin = -1.9, ymax = -0.1, fill = gene),
              alpha = 0.1)
  
  # Only add mutation lines and labels if there are mutations
  if (nrow(mutations_for_plot) > 0) {
    p_mutations <- p_mutations +
      # Add mutations as vertical lines
      geom_segment(data = mutations_for_plot,
                  aes(x = POS, xend = POS, y = -1.9, yend = -0.1, color = GENE_NAME, 
                      alpha = Allele_Frequency, linewidth = ifelse(highlight, 1.2, 0.8)),
                  lineend = "round")
    
    # Only add labels if there are any non-NA labels
    if (nrow(mutations_for_plot %>% filter(!is.na(label))) > 0) {
      p_mutations <- p_mutations +
        # Add mutation labels
        geom_text_repel(data = mutations_for_plot %>% filter(!is.na(label)),
                       aes(x = POS, y = -0.5, label = label, color = GENE_NAME),
                       size = 3, box.padding = 0.5, segment.color = "grey50",
                       min.segment.length = 0, max.overlaps = 30, direction = "y")
    }
  } else {
    # Add a note when no mutations are found
    p_mutations <- p_mutations +
      annotate("text", x = genome_features$genome_length/2, y = -1, 
               label = paste0("No mutations found with frequency ≥ ", cutoff*100, "%"), 
               size = 4, fontface = "italic")
  }
  
  p_mutations <- p_mutations +
                   
    # Set scales
    scale_x_continuous(name = NULL, 
                      breaks = seq(0, ceiling(genome_features$genome_length/1000)*1000, by = 1000),
                      labels = function(x) format(x, big.mark = ",")) +
    scale_y_continuous(name = NULL, limits = c(-2, 0)) +
    scale_color_manual(values = gene_colors, name = "Gene") +
    scale_alpha_continuous(range = c(0.5, 1), name = "Allele\nFrequency") +
    scale_linewidth_identity() +
    
    # Customize theme
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  # Combine plots
  p_combined <- p_genome / p_mutations + plot_layout(heights = c(1, 1.5))
  
  # Add title with definition
  accession <- unique(mutations$CHROM)[1]
  definition <- genome_features$definition
  # Truncate definition if too long (keep first 80 characters)
  if (nchar(definition) > 80) {
    definition <- paste0(substr(definition, 1, 77), "...")
  }
  
  final_plot <- p_combined + 
    plot_annotation(
      title = paste0("Mutations in ", accession, " - ", definition, " (cutoff: ", cutoff*100, "%)"),
      theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
    )
  
  return(final_plot)
}

# Main execution
main <- function() {
  # Read mutations file
  cat("Reading mutations file:", opt$input, "\n")
  mutations <- read_mutations(opt$input)
  
  # Get accession number
  if (!is.null(opt$accession)) {
    accession <- opt$accession
    cat("Using provided accession:", accession, "\n")
  } else {
    accession <- unique(mutations$CHROM)[1]
    cat("Using accession from CHROM column:", accession, "\n")
    
    # Check if it's a valid accession or needs to be replaced with a default
    if (accession == "CHROM" || accession == "") {
      # Not a valid accession, use a default WNV reference
      cat("CHROM doesn't appear to be a valid accession. Using default WNV reference: AY532665.1\n")
      accession <- "AY532665.1"  # West Nile virus reference genome for NY99
    }
  }
  
  # Get genome features
  genome_features <- get_genome_features(accession)
  
  # Map gene names in mutations to known WNV genes if needed
  # This handles cases where SnpEff uses different gene naming
  known_genes <- genome_features$genes$gene
  mutations$GENE_NAME_ORIG <- mutations$GENE_NAME
  
  for (i in 1:nrow(mutations)) {
    gene <- mutations$GENE_NAME[i]
    
    # Try to match with known genes
    if (!is.na(gene) && !(gene %in% known_genes)) {
      # Check if it's a gene ID that needs mapping
      for (j in 1:length(known_genes)) {
        if (grepl(gene, genome_features$genes$product[j], ignore.case = TRUE) ||
            grepl(gene, genome_features$genes$gene[j], ignore.case = TRUE)) {
          mutations$GENE_NAME[i] <- known_genes[j]
          break
        }
      }
    }
  }
  
  # Create visualization
  cat("Creating visualization...\n")
  plot <- create_genome_visualization(
    mutations, 
    genome_features, 
    opt$cutoff,
    opt$mutation_genes,
    opt$colors,
    opt$highlight_freq
  )
  
  # Save plot
  cat("Saving plot to:", opt$output, "\n")
  output_format <- tolower(tools::file_ext(opt$output))
  
  if (output_format == "pdf") {
    ggsave(opt$output, plot, width = opt$width, height = opt$height, device = cairo_pdf)
  } else if (output_format %in% c("png", "jpg", "jpeg", "tiff")) {
    ggsave(opt$output, plot, width = opt$width, height = opt$height, dpi = 300)
  } else {
    cat("Unrecognized output format. Saving as PDF...\n")
    ggsave(paste0(tools::file_path_sans_ext(opt$output), ".pdf"), 
           plot, width = opt$width, height = opt$height, device = cairo_pdf)
  }
  
  cat("Done!\n")
}

# Run the main function
main()