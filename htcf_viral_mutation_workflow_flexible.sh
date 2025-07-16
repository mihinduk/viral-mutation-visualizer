#\!/bin/bash
# HTCF-adapted viral mutation visualizer workflow
# This script demonstrates the complete workflow with QUAL score filtering

# Default values
DEFAULT_DEPTH=200
DEFAULT_VIS_CUTOFF=0.05
DEFAULT_ENV_NAME="viral_genomics"

# Example usage:
# ./htcf_viral_mutation_workflow.sh input.vcf reference.fasta accession min_qual_score output_prefix [vis_cutoff] [min_depth]

set -e

# Store arguments first to avoid conflicts
INPUT_VCF="$1"
REFERENCE="$2"
ACCESSION="$3"
MIN_QUAL="$4"
OUTPUT_PREFIX="$5"
VIS_CUTOFF="${6:-$DEFAULT_VIS_CUTOFF}"
MIN_DEPTH="${7:-$DEFAULT_DEPTH}"

# Environment variables (can be overridden)
CONDA_ENV="${VIRAL_ENV_NAME:-$DEFAULT_ENV_NAME}"
MAMBA_PATH="${MAMBA_PATH:-/home/mihindu/miniforge3/bin/mamba}"
CONDA_PATH="${CONDA_PATH:-/ref/sahlab/software/anaconda3/bin/conda}"

# Extract pipeline directory from the original script path
if [ -n "$SLURM_JOB_ID" ]; then
    # Extract original script path from SLURM job details
    ORIGINAL_SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID  < /dev/null |  grep -oP 'Command=\K[^ ]+')
    PIPELINE_DIR="$(dirname "$ORIGINAL_SCRIPT_PATH")"
else
    # Running outside SLURM
    PIPELINE_DIR="$(dirname "$0")"
fi

# Check arguments
if [ $# -lt 5 ]; then
    echo "Usage: htcf_viral_mutation_workflow.sh <input_vcf> <reference_fasta> <accession> <min_qual_score> <output_prefix> [vis_cutoff] [min_depth]"
    echo ""
    echo "Required arguments:"
    echo "  input_vcf       - Input VCF file from snpEff annotation"
    echo "  reference_fasta - Reference genome FASTA file"
    echo "  accession       - GenBank accession number (e.g., HM440560.1)"
    echo "  min_qual_score  - Minimum QUAL score for filtering (e.g., 39234)"
    echo "  output_prefix   - Prefix for output files"
    echo ""
    echo "Optional arguments:"
    echo "  vis_cutoff      - Visualization allele frequency cutoff (default: $DEFAULT_VIS_CUTOFF = 5%)"
    echo "  min_depth       - Minimum read depth for filtering (default: $DEFAULT_DEPTH)"
    echo ""
    echo "Environment variables (optional):"
    echo "  VIRAL_ENV_NAME  - Conda environment name (default: $DEFAULT_ENV_NAME)"
    echo "  MAMBA_PATH      - Path to mamba executable"
    echo "  CONDA_PATH      - Path to conda executable"
    echo ""
    echo "Example: htcf_viral_mutation_workflow.sh sample.vcf reference.fasta HM440560.1 49314 sample_filtered 0.01 100"
    exit 1
fi

echo "Starting viral mutation visualization pipeline..."
echo "Input VCF: $INPUT_VCF"
echo "Reference: $REFERENCE"
echo "Accession: $ACCESSION"
echo "Minimum QUAL score: $MIN_QUAL"
echo "Minimum depth: $MIN_DEPTH"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Visualization AF cutoff: $VIS_CUTOFF ($(echo "scale=1; $VIS_CUTOFF * 100" | bc)%)"
echo "Using conda environment: $CONDA_ENV"

# Create output directory
OUTPUT_DIR="${OUTPUT_PREFIX}_results"
mkdir -p "$OUTPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# Setup conda environment
echo ""
echo "Setting up conda environment..."

# Check if conda environment exists
if $CONDA_PATH env list | grep -q "$CONDA_ENV"; then
    echo "Environment $CONDA_ENV already exists"
else
    echo "Error: Environment $CONDA_ENV not found. Please create it first."
    exit 1
fi

# Use mamba run approach
MAMBA_CMD="$MAMBA_PATH run -n $CONDA_ENV"

# Step 1: Parse snpEff VCF and filter by QUAL score and depth
echo ""
echo "Step 1: Parsing snpEff VCF and filtering by QUAL >= $MIN_QUAL and depth >= $MIN_DEPTH..."
$MAMBA_CMD python3 "${PIPELINE_DIR}/parse_snpeff_vcf.py" -i "$INPUT_VCF" -d "$MIN_DEPTH" -q "$MIN_QUAL" -o "${OUTPUT_PREFIX}_mutations.tsv" -O "$OUTPUT_DIR"

# Find the actual filtered VCF file created by parse_snpeff_vcf.py
FILTERED_VCF=$(find "$OUTPUT_DIR" -name "*_filtered.vcf" 2>/dev/null | tail -1)
if [ -z "$FILTERED_VCF" ]; then
    echo "Error: Could not find filtered VCF file in $OUTPUT_DIR"
    exit 1
fi
echo "Using filtered VCF: $FILTERED_VCF"

# Step 2: Generate consensus sequence from filtered VCF
echo ""
echo "Step 2: Generating consensus sequence from filtered VCF..."
$MAMBA_CMD python3 "${PIPELINE_DIR}/vcf_to_consensus_M2.py" -i "$FILTERED_VCF" -r "$REFERENCE" -o "${OUTPUT_DIR}/${OUTPUT_PREFIX}_consensus.fa"

# Step 3: Generate protein sequences
echo ""
echo "Step 3: Generating protein sequences..."
$MAMBA_CMD python3 "${PIPELINE_DIR}/consensus_to_proteins.py" -c "${OUTPUT_DIR}/${OUTPUT_PREFIX}_consensus.fa" -p "${OUTPUT_PREFIX}" -a "$ACCESSION" -o "${OUTPUT_DIR}"

# Step 4: Visualize mutations
echo ""
echo "Step 4: Creating mutation visualization..."
# Find the actual TSV file created by parse_snpeff_vcf.py
TSV_FILE=$(find "$OUTPUT_DIR" -name "*_mutations*.tsv" 2>/dev/null | grep -v "_table.tsv" | tail -1)
if [ -z "$TSV_FILE" ]; then
    echo "Error: Could not find mutations TSV file in $OUTPUT_DIR"
    exit 1
fi
echo "Using mutations TSV: $TSV_FILE"

echo "NOTE: ALL mutations are saved in the TSV file. Only mutations with AF >= $VIS_CUTOFF are displayed in the plot."
$MAMBA_CMD python3 "${PIPELINE_DIR}/visualize_mutations_python.py" --input "$TSV_FILE" --output "${OUTPUT_DIR}/${OUTPUT_PREFIX}_mutations" --accession "$ACCESSION" --cutoff "$VIS_CUTOFF"

echo ""
echo "Pipeline complete\! Output files in $OUTPUT_DIR:"
echo "  - Filtered mutations TSV: $TSV_FILE (contains ALL mutations passing QUAL >= $MIN_QUAL and depth >= $MIN_DEPTH)"
echo "  - Filtered VCF: $FILTERED_VCF"
echo "  - Consensus sequence: ${OUTPUT_DIR}/${OUTPUT_PREFIX}_consensus.fa"
echo "  - Protein sequences: ${OUTPUT_DIR}/${OUTPUT_PREFIX}_*.fasta"
echo "  - Visualization: ${OUTPUT_DIR}/${OUTPUT_PREFIX}_mutations.png"
echo "    (showing only mutations with AF >= $VIS_CUTOFF)"
echo "  - Mutations table: ${OUTPUT_DIR}/${OUTPUT_PREFIX}_mutations_mutations_table.tsv"
echo "    (showing only mutations with AF >= $VIS_CUTOFF for the plot)"
