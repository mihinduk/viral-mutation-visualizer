#\!/bin/bash
# Test the updated visualization scripts with Powassan virus data

# Set paths
MUTATION_SUMMARY="/scratch/sahlab/kathie/2024_10_28_diamond_rerun/analysis/ivar_variant_analysis/powassan_samples_mutation_summary.csv"
DEPTH_DIR="/scratch/sahlab/kathie/2024_10_28_diamond_rerun/output/coverage_depth"
OUTPUT_DIR="test_output_20250718_083242"

# Create output directory
mkdir -p $OUTPUT_DIR

# Run mutation visualization
echo "Running mutation visualization..."
python visualize_mutations_python_outdir.py \
    -i $MUTATION_SUMMARY \
    -o $OUTPUT_DIR/mutations \
    --reference HM440560.1

# Run depth visualization for a sample
echo "Running depth visualization..."
SAMPLE="N1027-10_S10"  # Example sample
python visualize_depth.py \
    -d $DEPTH_DIR/${SAMPLE}_depth.txt \
    -o $OUTPUT_DIR/depth_${SAMPLE}.png \
    --reference HM440560.1

echo "Test complete\! Check output in: $OUTPUT_DIR"
ls -la $OUTPUT_DIR/

