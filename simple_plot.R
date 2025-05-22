#!/usr/bin/env Rscript

# A simplified script to visualize WNV mutations

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)

# Read data
input_file <- "/Users/handley_lab/Handley Lab Dropbox/virome/Diamond_lab_isolate_seq/2025_04_02_Diamond_multi_viral/WNV/WNV_Ug_NY/NovaSeq_N917_I13847_Cell_Culture_RNA_Diamond_Scheaffer_WNV_NY99_6_36_West_Nile_Virus.snpEFF.ann_200.tsv"
output_file <- "WNV_NY99_6_36_mutations.pdf"
cutoff <- 0.90

# Read the TSV file
cat("Reading input file:", input_file, "\n")
mutations <- read_tsv(input_file)

# Filter by frequency
filtered_mutations <- mutations %>% 
  filter(Allele_Frequency >= cutoff)

cat("Found", nrow(filtered_mutations), "mutations with frequency >=", cutoff, "\n")

# Define WNV genome structure
genes <- data.frame(
  gene = c("C", "prM", "Env", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5"),
  start = c(97, 466, 967, 2470, 3526, 4219, 4612, 6469, 6916, 7672),
  end = c(465, 966, 2469, 3525, 4218, 4611, 6468, 6915, 7671, 10395),
  category = c(rep("structural", 3), rep("non-structural", 7)),
  stringsAsFactors = FALSE
)

# Define UTRs
utrs <- data.frame(
  region = c("5'UTR", "3'UTR"),
  start = c(1, 10396),
  end = c(96, 11000),
  stringsAsFactors = FALSE
)

# Define colors for genes
gene_colors <- c(
  "C" = "#4575b4", 
  "prM" = "#74add1", 
  "Env" = "#abd9e9", 
  "NS1" = "#fdae61", 
  "NS2a" = "#f46d43", 
  "NS2b" = "#d73027", 
  "NS3" = "#a50026", 
  "NS4a" = "#762a83", 
  "NS4b" = "#9970ab", 
  "NS5" = "#c2a5cf"
)

# Create the genome layout plot
p_genome <- ggplot() +
  # Add genome background
  geom_rect(aes(xmin = 1, xmax = 11000, ymin = 0.45, ymax = 0.55),
            fill = NA, color = "black") +
  
  # Add UTR regions
  geom_rect(data = utrs, 
            aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55, fill = region),
            color = "black", fill = "grey80") +
  
  # Add gene regions
  geom_rect(data = genes,
            aes(xmin = start, xmax = end, 
                ymin = ifelse(category == "structural", 0.72, 0.12), 
                ymax = ifelse(category == "structural", 0.88, 0.28), 
                fill = gene),
            color = "black") +
  
  # Add category labels
  annotate("text", x = mean(c(genes$start[1], genes$end[3])), y = 1.0, 
           label = "Structural proteins", size = 4) +
  annotate("text", x = mean(c(genes$start[4], genes$end[10])), y = 0.0, 
           label = "Non-structural proteins", size = 4) +
  
  # Add gene labels
  geom_text(data = genes,
            aes(x = (start + end)/2, 
                y = ifelse(category == "structural", 0.8, 0.2), 
                label = gene),
            size = 3) +
  
  # Add scales
  scale_x_continuous(name = "Genome position (nt)", 
                     breaks = seq(0, 11000, by = 1000),
                     labels = function(x) format(x, big.mark = ",")) +
  scale_y_continuous(name = NULL, limits = c(-0.8, 1.2)) +
  scale_fill_manual(values = c(gene_colors, "5'UTR" = "grey80", "3'UTR" = "grey80")) +
  
  # Theme
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5)
  ) +
  ggtitle("WNV NY99 Genome Structure")

# Process mutations for plotting
mutations_plot_data <- filtered_mutations %>%
  mutate(
    position = POS,
    height = Allele_Frequency,
    gene_name = GENE_NAME
  )

# Create the mutations plot 
p_mutations <- ggplot() +
  # Add gene regions as background
  geom_rect(data = genes,
            aes(xmin = start, xmax = end, ymin = -0.7, ymax = -0.1, fill = gene),
            alpha = 0.1) +
            
  # Set scales
  scale_x_continuous(name = "Genome position (nt)", 
                    breaks = seq(0, 11000, by = 1000),
                    labels = function(x) format(x, big.mark = ",")) +
  scale_y_continuous(name = "Allele Frequency", limits = c(-0.7, 1)) +
  scale_fill_manual(values = gene_colors, name = "Gene") +
  
  # Theme
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  ggtitle(paste0("Mutations with frequency ≥ ", cutoff*100, "%"))

# Add mutations if available
if (nrow(mutations_plot_data) > 0) {
  # Add mutation lines
  p_mutations <- p_mutations +
    geom_segment(data = mutations_plot_data,
                aes(x = position, xend = position, y = -0.7, yend = height,
                    color = gene_name),
                alpha = 0.7, size = 0.5) +
    geom_point(data = mutations_plot_data,
               aes(x = position, y = height, color = gene_name),
               alpha = 0.8, size = 2) +
    scale_color_manual(values = gene_colors)
} else {
  # Add a note when no mutations are found
  p_mutations <- p_mutations +
    annotate("text", x = 5500, y = 0, 
             label = paste0("No mutations found with frequency ≥ ", cutoff*100, "%"), 
             size = 4, fontface = "italic")
}

# Combine plots
library(gridExtra)
p_combined <- grid.arrange(p_genome, p_mutations, nrow = 2, heights = c(1, 1.5))

# Save the plot
cat("Saving plot to:", output_file, "\n")
ggsave(output_file, p_combined, width = 10, height = 8)
cat("Done!\n")