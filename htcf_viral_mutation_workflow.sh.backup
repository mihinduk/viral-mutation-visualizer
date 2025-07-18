#\!/bin/bash
# HTCF-adapted viral mutation visualizer workflow
# This script demonstrates the complete workflow with QUAL score filtering

# Example usage:
# ./htcf_viral_mutation_workflow.sh input.vcf reference.fasta accession min_qual_score output_prefix [vis_cutoff]

set -e

# Store arguments first to avoid conflicts
INPUT_VCF="$1"
REFERENCE="$2"
ACCESSION="$3"
MIN_QUAL="$4"
OUTPUT_PREFIX="$5"
VIS_CUTOFF="${6:-0.05}"  # Default to 0.05 (5%) if not specified

# Extract pipeline directory from the original script path
# Use the working solution from CLAUDE.md for SLURM path resolution
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
    echo "Usage: htcf_viral_mutation_workflow.sh <input_vcf> <reference_fasta> <accession> <min_qual_score> <output_prefix> [vis_cutoff]"
    echo ""
    echo "Required arguments:"
    echo "  input_vcf       - Input VCF file from snpEff annotation"
    echo "  reference_fasta - Reference genome FASTA file"
    echo "  accession       - GenBank accession number (e.g., HM440560.1)"
    echo "  min_qual_score  - Minimum QUAL score for filtering (e.g., 39234)"
    echo "  output_prefix   - Prefix for output files"
    echo ""
    echo "Optional arguments:"
    echo "  vis_cutoff      - Visualization allele frequency cutoff (default: 0.05 = 5%)"
    echo "                    Note: ALL mutations are saved in the TSV file regardless of this cutoff"
    echo ""
    echo "Example: htcf_viral_mutation_workflow.sh sample.vcf reference.fasta HM440560.1 49314 sample_filtered 0.01"
    exit 1
fi

echo "Starting viral mutation visualization pipeline..."
echo "Input VCF: $INPUT_VCF"
echo "Reference: $REFERENCE"
echo "Accession: $ACCESSION"
echo "Minimum QUAL score: $MIN_QUAL"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Visualization AF cutoff: $VIS_CUTOFF ($(echo "scale=1; $VIS_CUTOFF * 100" | bc)%)"

# Create output directory
OUTPUT_DIR="${OUTPUT_PREFIX}_results"
mkdir -p "$OUTPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# Setup conda environment using mamba run approach like shotgun_viral_genomics
echo ""
echo "Setting up conda environment..."

# Create environment if it doesn't exist using full path to conda
if /ref/sahlab/software/anaconda3/bin/conda env list | grep -q "viral_genomics"; then
    echo "Environment viral_genomics already exists"
else
    echo "Creating conda environment..."
    /ref/sahlab/software/anaconda3/bin/conda env create -f "${PIPELINE_DIR}/environment_M2.yml"
fi

# Use mamba run approach instead of conda activate
MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics"

# Step 1: Parse snpEff VCF and filter by QUAL score
echo ""
echo "Step 1: Parsing snpEff VCF and filtering by QUAL >= $MIN_QUAL..."
$MAMBA_CMD python3 "${PIPELINE_DIR}/parse_snpeff_vcf.py" -i "$INPUT_VCF" -d 200 -q "$MIN_QUAL" -o "${OUTPUT_PREFIX}_mutations.tsv" -O "$OUTPUT_DIR"

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
TSV_FILE=$(find "$OUTPUT_DIR" -name "*_mutations*.tsv" 2>/dev/null | tail -1)
if [ -z "$TSV_FILE" ]; then
    echo "Error: Could not find mutations TSV file in $OUTPUT_DIR"
    exit 1
fi
echo "Using mutations TSV: $TSV_FILE"

# Check if we should use R or Python version
if command -v Rscript &> /dev/null && [ -f "${PIPELINE_DIR}/visualize_mutations_dynamic.R" ]; then
    echo "Using R visualization..."
    $MAMBA_CMD Rscript "${PIPELINE_DIR}/visualize_mutations_dynamic.R" "$TSV_FILE" "${OUTPUT_DIR}/${OUTPUT_PREFIX}_mutations.pdf"
else
    echo "Using Python visualization..."
    echo "NOTE: ALL mutations are saved in the TSV file. Only mutations with AF >= $VIS_CUTOFF are displayed in the plot."
    $MAMBA_CMD python3 "${PIPELINE_DIR}/visualize_mutations_python.py" --input "$TSV_FILE" --output "${OUTPUT_DIR}/${OUTPUT_PREFIX}_mutations" --accession "$ACCESSION" --cutoff "$VIS_CUTOFF"
fi

echo ""
echo "Pipeline complete\! Output files in $OUTPUT_DIR:"
echo "  - Filtered mutations TSV: $TSV_FILE (contains ALL mutations passing QUAL/depth filters)"
echo "  - Filtered VCF: $FILTERED_VCF"
echo "  - Consensus sequence: ${OUTPUT_DIR}/${OUTPUT_PREFIX}_consensus.fa"
echo "  - Protein sequences: ${OUTPUT_DIR}/${OUTPUT_PREFIX}_*.fasta"
echo "  - Visualization: ${OUTPUT_DIR}/${OUTPUT_PREFIX}_mutations.pdf (or .png/.svg)"
echo "    (showing only mutations with AF >= $VIS_CUTOFF)"
echo "  - Mutations table: ${OUTPUT_DIR}/${OUTPUT_PREFIX}_mutations_mutations_table.tsv"
echo "    (showing only mutations with AF >= $VIS_CUTOFF for the plot)"
