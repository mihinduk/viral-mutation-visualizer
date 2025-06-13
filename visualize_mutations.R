#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# visualize_mutations_v5.R - Non-synonymous mutations only with separate output table
# -------------------------------------------------------------------------

# Function to install and load required packages
install_and_load_packages <- function() {
  # Define required packages
  cran_packages <- c("optparse", "ggplot2", "dplyr", "tidyr", "readr", "stringr", 
                     "patchwork", "ggrepel", "RColorBrewer", "gridExtra", "grid")
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
  definition <- "West Nile virus, complete genome"
  
  # Calculate UTR regions
  first_cds_start <- min(genes_data$start)
  last_cds_end <- max(genes_data$end)
  
  utr_data <- data.frame(
    region = c("5'UTR", "3'UTR"),
    start = c(1, last_cds_end + 1),
    end = c(first_cds_start - 1, genome_length),
    stringsAsFactors = FALSE
  )
  
  return(list(
    genes = genes_data,
    utrs = utr_data,
    genome_length = genome_length,
    definition = definition
  ))
}

# Function to format amino acid change
format_aa_change <- function(hgvsp) {
  # Vectorized function to handle multiple values
  sapply(hgvsp, function(x) {
    if (is.na(x) || x == "") return("")
    
    # Remove "p." prefix
    aa_change <- gsub("^p\\.", "", x)
  
  # Keep three-letter amino acid codes but add spaces
  # This regex adds a space before any digit that follows letters
  # and a space after any digit that precedes letters
  aa_change <- gsub("([A-Za-z])(\\d)", "\\1 \\2", aa_change)
  aa_change <- gsub("(\\d)([A-Za-z])", "\\1 \\2", aa_change)
  
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
          # Add space between amino acid letter and position number
          aa_display <- gsub("([A-Z])([0-9])", "\\1 \\2", aa_display)
          entry <- paste0(entry, "\t", aa_display)
        } else {
          # Add space between amino acid letter and position number
          aa_display <- gsub("([A-Z])([0-9])", "\\1 \\2", gene_muts$aa_change[i])
          entry <- paste0(entry, "\t", aa_display)
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

# Function to create genome-only visualization (no mutations or no non-synonymous mutations)
create_genome_only_visualization <- function(genome_features, mutations, cutoff, message) {
  # Define colors
  gene_names <- c("C", "prM", "Env", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5")
  color_values <- c("#4575b4", "#74add1", "#abd9e9", "#fdae61", "#f46d43", 
                    "#d73027", "#a50026", "#762a83", "#9970ab", "#c2a5cf")
  gene_colors <- setNames(color_values, gene_names)
  
  # Create the main plot area
  grid.newpage()
  
  # Define layout with space for message
  pushViewport(viewport(layout = grid.layout(
    nrow = 3, 
    ncol = 1,
    heights = unit(c(0.3, 1.5, 0.5), "null")
  )))
  
  # Add title at top
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  accession <- unique(mutations$CHROM)[1]
  grid.text(paste0("Mutations in ", accession, " - ", 
                   genome_features$definition, " (cutoff: ", cutoff*100, "%)"),
            gp = gpar(fontsize = 14, fontface = "bold"))
  popViewport()
  
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
  
  str_range <- c(min(str_genes$start), max(str_genes$end))
  non_str_range <- c(min(non_str_genes$start), max(non_str_genes$end))
  
  # Create the genome plot
  p_genome <- ggplot() +
    # Main genome line
    geom_segment(aes(x = 1, xend = genome_features$genome_length, 
                     y = genome_line_y, yend = genome_line_y),
                 linewidth = 3, color = "black") +
    
    # UTR regions
    geom_segment(data = genome_features$utrs,
                 aes(x = start, xend = end, y = genome_line_y, yend = genome_line_y),
                 linewidth = 5, color = "grey60") +
    geom_text(data = genome_features$utrs, 
              aes(x = (start + end)/2, y = genome_line_y + 0.1, label = region), 
              size = 3) +
    
    # Gene rectangles
    geom_rect(data = genes_for_plot,
              aes(xmin = start, xmax = end, 
                  ymin = y - 0.06, ymax = y + 0.06, fill = gene),
              color = "white", linewidth = 0.5) +
    
    # Gene labels directly on genes
    geom_text(data = genes_for_plot,
              aes(x = x_center, y = y, label = gene),
              size = 4, fontface = "bold", color = "white") +
    
    # Structural/Non-structural labels
    geom_segment(aes(x = str_range[1], xend = str_range[2], 
                     y = genome_line_y - 0.15, yend = genome_line_y - 0.15),
                 color = "#2166ac", linewidth = 3, alpha = 0.5) +
    geom_text(aes(x = mean(str_range), y = genome_line_y - 0.22), 
              label = "Structural proteins", size = 4, fontface = "bold", color = "#2166ac") +
    
    geom_segment(aes(x = non_str_range[1], xend = non_str_range[2], 
                     y = genome_line_y - 0.15, yend = genome_line_y - 0.15),
                 color = "#762a83", linewidth = 3, alpha = 0.5) +
    geom_text(aes(x = mean(non_str_range), y = genome_line_y - 0.22), 
              label = "Non-structural proteins", size = 4, fontface = "bold", color = "#762a83") +
    
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
  
  print(p_genome, newpage = FALSE)
  popViewport()
  
  # Add message at bottom
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
  grid.text(message, gp = gpar(fontsize = 16, fontface = "italic", col = "darkred"))
  popViewport()
  
  popViewport()  # Main layout
}

# Function to create genome visualization with multiple tables
create_genome_visualization <- function(mutations, genome_features, cutoff, gene_selection, custom_colors = NULL) {
  # Filter mutations by allele frequency
  filtered_mutations <- mutations %>% filter(Allele_Frequency >= cutoff)
  
  # Check if we have any mutations after filtering
  if (nrow(filtered_mutations) == 0) {
    cat("No mutations found with frequency >=", cutoff*100, "%\n")
    # Create genome-only visualization
    create_genome_only_visualization(genome_features, mutations, cutoff, "No mutations at this Allele Frequency")
    return(invisible(NULL))
  }
  
  # Check for non-synonymous mutations to determine what message to show
  nonsyn_check <- filtered_mutations %>%
    filter(
      grepl("missense|nonsense|stop_gained|stop_lost", EFFECT, ignore.case = TRUE) |
      (!is.na(HGVSp) & HGVSp != "" & !grepl("=$", HGVSp)) |
      grepl("HIGH", PUTATIVE_IMPACT, ignore.case = TRUE)
    ) %>%
    filter(!grepl("^synonymous_variant$", EFFECT, ignore.case = TRUE))
  
  # If no non-synonymous mutations but we have mutations, show special message
  if (nrow(nonsyn_check) == 0) {
    cat("No non-synonymous mutations found with frequency >=", cutoff*100, "%\n")
    create_genome_only_visualization(genome_features, mutations, cutoff, "No non-synonymous mutations at this Allele Frequency")
    return(invisible(NULL))
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
      for (j in 1:nrow(genome_features$genes)) {
        if (pos >= genome_features$genes$start[j] && pos <= genome_features$genes$end[j]) {
          filtered_mutations$GENE_NAME[i] <- genome_features$genes$gene[j]
          mapped <- TRUE
          break
        }
      }
      
      if (!mapped) {
        # Check if in UTR regions
        if (pos < min(genome_features$genes$start)) {
          filtered_mutations$GENE_NAME[i] <- "5'UTR"
        } else if (pos > max(genome_features$genes$end)) {
          filtered_mutations$GENE_NAME[i] <- "3'UTR"
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
  
  # Define colors
  gene_names <- c("C", "prM", "Env", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5")
  
  if (!is.null(custom_colors)) {
    color_values <- unlist(strsplit(custom_colors, ","))
    if (length(color_values) < length(gene_names)) {
      color_values <- brewer.pal(length(gene_names), "Paired")
    }
  } else {
    color_values <- c("#4575b4", "#74add1", "#abd9e9", "#fdae61", "#f46d43", 
                      "#d73027", "#a50026", "#762a83", "#9970ab", "#c2a5cf")
  }
  
  gene_colors <- setNames(color_values, gene_names)
  
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
  
  str_range <- c(min(str_genes$start), max(str_genes$end))
  non_str_range <- c(min(non_str_genes$start), max(non_str_genes$end))
  
  # Create the genome plot
  p_genome <- ggplot() +
    # Main genome line
    geom_segment(aes(x = 1, xend = genome_features$genome_length, 
                     y = genome_line_y, yend = genome_line_y),
                 linewidth = 3, color = "black") +
    
    # UTR regions
    geom_segment(data = genome_features$utrs,
                 aes(x = start, xend = end, y = genome_line_y, yend = genome_line_y),
                 linewidth = 5, color = "grey60") +
    geom_text(data = genome_features$utrs, 
              aes(x = (start + end)/2, y = genome_line_y + 0.1, label = region), 
              size = 3) +
    
    # Gene rectangles
    geom_rect(data = genes_for_plot,
              aes(xmin = start, xmax = end, 
                  ymin = y - 0.06, ymax = y + 0.06, fill = gene),
              color = "white", linewidth = 0.5) +
    
    # Gene labels directly on genes
    geom_text(data = genes_for_plot,
              aes(x = x_center, y = y, label = gene),
              size = 4, fontface = "bold", color = "white") +
    
    # Structural/Non-structural labels
    geom_segment(aes(x = str_range[1], xend = str_range[2], 
                     y = genome_line_y - 0.15, yend = genome_line_y - 0.15),
                 color = "#2166ac", linewidth = 3, alpha = 0.5) +
    geom_text(aes(x = mean(str_range), y = genome_line_y - 0.22), 
              label = "Structural proteins", size = 4, fontface = "bold", color = "#2166ac") +
    
    geom_segment(aes(x = non_str_range[1], xend = non_str_range[2], 
                     y = genome_line_y - 0.15, yend = genome_line_y - 0.15),
                 color = "#762a83", linewidth = 3, alpha = 0.5) +
    geom_text(aes(x = mean(non_str_range), y = genome_line_y - 0.22), 
              label = "Non-structural proteins", size = 4, fontface = "bold", color = "#762a83") +
    
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
  
  # Only process gene mapping and write table if there are mutations
  if (nrow(filtered_mutations) > 0) {
    # Map gene names before writing table
    for (i in 1:nrow(filtered_mutations)) {
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
      for (j in 1:nrow(genome_features$genes)) {
        if (pos >= genome_features$genes$start[j] && pos <= genome_features$genes$end[j]) {
          filtered_mutations$GENE_NAME[i] <- genome_features$genes$gene[j]
          mapped <- TRUE
          break
        }
      }
      
      if (!mapped) {
        # Check if in UTR regions
        if (pos < min(genome_features$genes$start)) {
          filtered_mutations$GENE_NAME[i] <- "5'UTR"
        } else if (pos > max(genome_features$genes$end)) {
          filtered_mutations$GENE_NAME[i] <- "3'UTR"
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
  } else {
    cat("No mutations to write to table at cutoff", opt$cutoff*100, "%\n")
  }
  
  cat("Done!\n")
}

# Run the main function
main()