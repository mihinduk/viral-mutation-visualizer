#!/bin/bash

# Example usage for visualize_mutations.R

# Make the script executable
chmod +x visualize_mutations.R

# Example 1: Display all mutations above 5% frequency
./visualize_mutations.R \
  --input mutations.tsv \
  --output all_mutations.pdf \
  --cutoff 0.05 \
  --mutation-genes all

# Example 2: Display only mutations in structural proteins above 10% frequency
./visualize_mutations.R \
  --input mutations.tsv \
  --output structural_mutations.pdf \
  --cutoff 0.1 \
  --mutation-genes structural

# Example 3: Display only mutations in non-structural proteins above 10% frequency
./visualize_mutations.R \
  --input mutations.tsv \
  --output non_structural_mutations.pdf \
  --cutoff 0.1 \
  --mutation-genes non-structural

# Example 4: Display mutations in specific genes
./visualize_mutations.R \
  --input mutations.tsv \
  --output specific_genes_mutations.pdf \
  --cutoff 0.1 \
  --mutation-genes NS3,NS5

# Example 5: Customize plot appearance
./visualize_mutations.R \
  --input mutations.tsv \
  --output custom_style_mutations.pdf \
  --cutoff 0.05 \
  --width 16 \
  --height 10 \
  --highlight-freq 0.4 \
  --colors "#4575b4,#74add1,#abd9e9,#fdae61,#f46d43,#d73027,#a50026,#762a83,#9970ab,#c2a5cf"

# Example 6: Generate PNG output instead of PDF
./visualize_mutations.R \
  --input mutations.tsv \
  --output mutations.png \
  --cutoff 0.05 \
  --mutation-genes all

# Example 7: Manually specify an accession number
./visualize_mutations.R \
  --input mutations.tsv \
  --output mutations_custom_accession.pdf \
  --accession NC_009942.1 \
  --cutoff 0.05 \
  --mutation-genes all