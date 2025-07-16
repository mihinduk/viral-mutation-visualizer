#\!/bin/bash
# Quick test script to verify HTCF setup

echo "Testing HTCF viral mutation visualizer setup..."

# Test conda activation
echo "1. Testing conda activation..."
source /ref/sahlab/software/anaconda3/bin/activate
if [ $? -eq 0 ]; then
    echo "   ✓ Conda activation successful"
else
    echo "   ✗ Conda activation failed"
    exit 1
fi

# Test environment creation
echo "2. Testing environment creation..."
conda env create -f environment.yml --force
if [ $? -eq 0 ]; then
    echo "   ✓ Environment creation successful"
else
    echo "   ✗ Environment creation failed"
    exit 1
fi

# Test R packages installation
echo "3. Testing R packages..."
conda activate viral-mutation-viz
Rscript fix_R_packages.R
if [ $? -eq 0 ]; then
    echo "   ✓ R packages installation successful"
else
    echo "   ✗ R packages installation failed"
    exit 1
fi

# Test with sample data
echo "4. Testing visualization with sample data..."
./visualize_mutations.R --input sample_mutations.tsv --output test_output.pdf --cutoff 0.05 --mutation-genes all
if [ $? -eq 0 ]; then
    echo "   ✓ Visualization test successful"
    echo "   Output file: test_output.pdf"
else
    echo "   ✗ Visualization test failed"
    exit 1
fi

echo ""
echo "🎉 All tests passed\! The pipeline is ready for HTCF."
echo ""
echo "Usage examples:"
echo "1. Interactive test: ./test_setup.sh"
echo "2. SLURM submission: sbatch /full/path/to/submit_viral_mutation_viz.sh input.tsv output.pdf 0.05"
