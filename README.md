# Viral Mutation Visualization

This tool creates visualizations of viral genome mutations from SnpEff-annotated LoFreq output files. It displays mutations in the context of the viral genome structure, highlighting gene boundaries, coding regions, and mutation locations.

## Available Versions

- **visualize_mutations.R**: Original version optimized for West Nile virus (WNV) with hardcoded genome structure
- **visualize_mutations_dynamic.R**: Dynamic version that fetches genome information from NCBI for any virus

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
# For West Nile virus (fast, no internet required)
Rscript visualize_mutations.R --input <input-file> --output <output-file> [options]

# For any virus (requires internet connection for NCBI lookup)
Rscript visualize_mutations_dynamic.R --input <input-file> --output <output-file> [options]
```

### Note on Rscript Path

On some systems (especially macOS), you may need to use the full path to Rscript:

```
/usr/local/bin/Rscript visualize_mutations.R --input <input-file> --output <output-file> [options]
```

To find where Rscript is installed on your system:
```
which Rscript
```

Alternatively, you can make the scripts directly executable:
```
chmod +x visualize_mutations.R visualize_mutations_dynamic.R
./visualize_mutations.R --input <input-file> --output <output-file> [options]
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
# Option 1: Using Rscript directly (may need full path on some systems)
Rscript visualize_mutations.R --input mutations.tsv --output mutation_plot.pdf --cutoff 0.05 --mutation-genes all

# Option 2: Using the executable script
./visualize_mutations.R --input mutations.tsv --output mutation_plot.pdf --cutoff 0.05 --mutation-genes all
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

## Additional Tools

### parse_snpeff_vcf.py

This script filters SnpEff-annotated VCF files by depth and allele frequency, creating TSV files suitable for visualization.

#### Usage

```bash
python3 parse_snpeff_vcf.py -i <input.vcf> -d <depth> [-f <frequency>] [-q <qual>] [-o <output.tsv>] [-O <output-dir>]
```

#### Arguments

- `-i/--input`: Input VCF file from SnpEff annotation (required)
- `-d/--depth`: Minimum depth requirement (required)
- `-f/--min-freq`: Minimum allele frequency for filtered VCF output (optional)
- `-q/--min-qual`: Minimum QUAL score for filtering (optional)
- `-o/--output`: Output TSV filename (default: auto-generated)
- `-O/--output-dir`: Output directory (default: same as input file)

#### Example

```bash
# Filter by depth 200 and create TSV
python3 parse_snpeff_vcf.py -i sample.snpEFF.ann.vcf -d 200

# Filter by depth 200 and QUAL score >= 49314
python3 parse_snpeff_vcf.py -i sample.snpEFF.ann.vcf -d 200 -q 49314

# Also create filtered VCF with AF >= 0.9
python3 parse_snpeff_vcf.py -i sample.snpEFF.ann.vcf -d 200 -f 0.9 -O results/

# Combine all filters: depth, QUAL, and frequency
python3 parse_snpeff_vcf.py -i sample.snpEFF.ann.vcf -d 200 -q 49314 -f 0.9
```

### consensus_to_proteins.py

This script extracts individual protein sequences from a consensus genome sequence using GenBank annotations.

#### Installation

Requires Biopython:
```bash
pip3 install biopython
```

#### Usage

```bash
python3 consensus_to_proteins.py -c <consensus.fa> [-g <genbank.gb> | -a <accession>] -p <prefix> [-o <output-dir>]
```

#### Arguments

- `-c/--consensus`: Consensus genome FASTA file (required)
- `-g/--genbank`: GenBank file with annotations (either this or -a required)
- `-a/--accession`: GenBank accession number to fetch (either this or -g required)
- `-p/--prefix`: Prefix for output protein files (required)
- `-o/--output-dir`: Output directory for protein FASTA files (default: current directory)

#### Example

```bash
# Using local GenBank file
python3 consensus_to_proteins.py -c consensus.fa -g genes.gbk -p WNV_sample -o proteins/

# Using GenBank accession
python3 consensus_to_proteins.py -c consensus.fa -a AY532665 -p WNV_sample -o proteins/
```

For flaviviruses with polyproteins, the script automatically cleaves them into individual proteins (C, prM, E, NS1, NS2A, NS2B, NS3, NS4A, NS4B, NS5).

### vcf_to_consensus.py

This script generates a consensus FASTA sequence from a VCF file by applying variants to a reference genome.

#### Requirements

Requires bcftools (includes bgzip and tabix). **Note: bcftools must be installed separately as it's not a Python package.**

**Recommended: Use Conda for easy installation**
```bash
# Create environment from provided file
conda env create -f environment.yml
conda activate viral-mutation-viz

# Or install just bcftools in existing environment
conda install -c bioconda bcftools
```

**Alternative installation methods:**
```bash
# macOS with Homebrew
brew install bcftools

# Linux (Debian/Ubuntu)
sudo apt-get install bcftools

# Linux (RHEL/CentOS)
sudo yum install bcftools
```

The script will check for bcftools and provide installation instructions if not found.

#### Usage

```bash
python3 vcf_to_consensus.py -i <input.vcf> [-r <reference.fa> | -a <accession>] [-o <output.fa>]
```

#### Arguments

- `-i/--input`: Input VCF file (required)
- `-r/--reference`: Local reference FASTA file (either this or -a required)
- `-a/--accession`: GenBank accession to fetch reference (either this or -r required)
- `-o/--output`: Output consensus FASTA file (default: <input>_consensus.fa)
- `--keep-temp`: Keep temporary compressed files

#### Example

```bash
# Using local reference
python3 vcf_to_consensus.py -i filtered.vcf -r reference.fa -o consensus.fa

# Using GenBank accession
python3 vcf_to_consensus.py -i filtered.vcf -a AY532665 -o consensus.fa

# Auto-generate output filename
python3 vcf_to_consensus.py -i sample_200_AF_0.9.vcf -a HM152775
```

## Complete Workflow Example

1. Filter VCF by depth and frequency:
```bash
python3 parse_snpeff_vcf.py -i sample.snpEFF.ann.vcf -d 200 -f 0.9
```

2. Create consensus sequence from filtered VCF:
```bash
python3 vcf_to_consensus.py -i sample_200_AF_0.9.vcf -a AY532665 -o consensus.fa
```

3. Extract protein sequences:
```bash
python3 consensus_to_proteins.py -c consensus.fa -a AY532665 -p sample -o proteins/
```

4. Visualize mutations:
```bash
Rscript visualize_mutations.R --input sample_200.tsv --output mutations.pdf --cutoff 0.05
```

## Known Limitations

- The tool is optimized for viral genomes, particularly flaviviruses like WNV
- Very large numbers of mutations may cause label overcrowding
- Automatically matches genes based on name patterns, which may require manual correction for some viruses