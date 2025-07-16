#\!/bin/bash
#SBATCH --job-name=viral_mut_viz
#SBATCH --output=viral_mut_viz_%j.out
#SBATCH --error=viral_mut_viz_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=sahlab

# Extract pipeline directory from the original sbatch command
ORIGINAL_SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID  < /dev/null |  grep -oP "Command=\K[^ ]+")
PIPELINE_DIR="$(dirname "$ORIGINAL_SCRIPT_PATH")"

# Activate conda using HTCF-specific path
source /ref/sahlab/software/anaconda3/bin/activate

# Create and activate the viral-mutation-viz environment
conda env create -f "${PIPELINE_DIR}/environment.yml" --force
conda activate viral-mutation-viz

# Install R packages
Rscript "${PIPELINE_DIR}/fix_R_packages.R"

# Example usage - customize these parameters
INPUT_FILE="${1:-sample_mutations.tsv}"
OUTPUT_FILE="${2:-mutations_output.pdf}"
CUTOFF="${3:-0.05}"

# Run the visualization
"${PIPELINE_DIR}/visualize_mutations.R" \
    --input "${INPUT_FILE}" \
    --output "${OUTPUT_FILE}" \
    --cutoff "${CUTOFF}" \
    --mutation-genes all

echo "Viral mutation visualization completed successfully\!"
