# Viral Mutation Visualization

This tool creates visualizations of viral genome mutations from SnpEff-annotated LoFreq output files. It displays mutations in the context of the viral genome structure, highlighting gene boundaries, coding regions, and mutation locations.

## Features

- Automatic retrieval of genome information from GenBank using the accession number
- Visualization of complete viral genome organization (structural and non-structural proteins)
- Display of mutations as vertical lines with customizable frequency cutoffs
- Filtering options for specific genes or gene categories (structural vs. non-structural)
- Customizable visualization parameters (colors, dimensions, etc.)
- Support for multiple output formats (PDF, PNG, JPG, TIFF)

## Installation

### Requirements

The script requires R (version 4.0 or higher) with the following packages:

```
optparse
ggplot2
dplyr
readr
rentrez
genbankr
seqinr
stringr
patchwork
ggrepel
RColorBrewer
```

### Installing Dependencies

This project uses `renv` for package management. To set up the required packages:

1. Install R (version 4.0 or higher)
2. Run the setup script:

```R
# From R console
source("setup_renv.R")
```

This will install all necessary packages in a project-specific library using `renv`.

#### Manual Installation (Alternative)

If you prefer not to use `renv`, you can install packages manually:

```R
install.packages(c("optparse", "ggplot2", "dplyr", "readr", "stringr", "patchwork", "ggrepel", "RColorBrewer"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("rentrez", "seqinr"))
```

## Usage

The script is run from the command line with several options:

```
Rscript visualize_mutations.R --input <input-file> --output <output-file> [options]
```

### Required Arguments

- `--input`: Path to the TSV file with SnpEff-annotated LoFreq output
- `--output`: Path for the output visualization file (PDF, PNG, JPG, or TIFF)

### Optional Arguments

- `--cutoff`: Allele frequency cutoff for filtering mutations (default: 0.05)
- `--mutation-genes`: Genes to display mutations for (default: "all")
  - Options: "all", "structural", "non-structural", or specific gene names (comma-separated)
- `--colors`: Custom color scheme (comma-separated hex values)
- `--width`: Plot width in inches (default: 14)
- `--height`: Plot height in inches (default: 8)
- `--highlight-freq`: Highlight mutations with frequency above this threshold (default: 0.5)
- `--accession`: Manually specify GenBank accession (overrides accession in CHROM column)

## Example Usage

### Display All Mutations Above 5% Frequency

```
Rscript visualize_mutations.R --input mutations.tsv --output mutation_plot.pdf --cutoff 0.05 --mutation-genes all
```

### Display Only Mutations in Structural Proteins Above 10% Frequency

```
Rscript visualize_mutations.R --input mutations.tsv --output mutation_plot.pdf --cutoff 0.1 --mutation-genes structural
```

### Display Mutations in Specific Genes

```
Rscript visualize_mutations.R --input mutations.tsv --output mutation_plot.pdf --cutoff 0.1 --mutation-genes NS3,NS5
```

### Customize Plot Appearance

```
Rscript visualize_mutations.R --input mutations.tsv --output mutation_plot.pdf --cutoff 0.05 --width 16 --height 10 --highlight-freq 0.4 --colors "#FF0000,#00FF00,#0000FF,#FFFF00,#FF00FF,#00FFFF,#000000,#FFFFFF,#888888,#444444"
```

## Input File Format

The script expects a TSV file with the following columns (from SnpEff-annotated LoFreq output):

```
CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, Total_Depth, Allele_Frequency, strand_bias, DP4,
EFFECT, PUTATIVE_IMPACT, GENE_NAME, GENE_ID, FEATURE_TYPE, FEATURE_ID, TRANSCRIPT_TYPE,
HGVSc, HGVSp, cDNA_POSITION_AND_LENGTH, CDS_POSITION_AND_LENGTH, PROTEIN_POSITION_AND_LENGTH, ERROR
```

## Output Visualization

The output visualization includes:

1. Top panel: Genome organization
   - Complete viral genome layout with 5' and 3' UTRs
   - Structural and non-structural protein regions
   - Individual gene locations with proper boundaries

2. Bottom panel: Mutation mapping
   - Vertical lines representing mutation positions
   - Color-coded by gene
   - Opacity reflecting allele frequency
   - Labels for significant mutations (above highlight threshold)
   - Different line weights for highlighting important mutations

## Troubleshooting

- If the GenBank fetching fails, you can manually specify the accession with `--accession`
- If gene names in your TSV don't match the standard viral gene names, the tool attempts to map them automatically
- For large genomes, consider adjusting the plot dimensions with `--width` and `--height`
- If mutation labels overlap too much, adjust the `max.overlaps` parameter in the script

## Known Limitations

- The tool is optimized for viral genomes, particularly flaviviruses like WNV
- Very large numbers of mutations may cause label overcrowding
- Automatically matches genes based on name patterns, which may require manual correction for some viruses