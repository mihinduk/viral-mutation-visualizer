# Python Version - No More R Package Headaches! 🐍

## Quick Start (For Bench Scientists)

**The Python version solves all the R package installation issues!** 

### Step 1: Install Python packages (ONE TIME ONLY)
```bash
pip install pandas numpy matplotlib
```

Or if you use conda:
```bash
conda install pandas numpy matplotlib
```

### Step 2: Use the Python script instead of R
```bash
# Instead of this (R version with package issues):
# Rscript visualize_mutations.R --input data.tsv --output plot.pdf

# Use this (Python version that just works):
python visualize_mutations_python.py --input data.tsv --output plot.png --cutoff 0.05
```

## Why Python Version is Better

✅ **No R package installation nightmares**  
✅ **Works on M2 Macs without issues**  
✅ **Identical visual output to R version**  
✅ **Only 3 simple Python packages needed**  
✅ **Much easier for bench scientists**  

## Command Line Usage

```bash
# Basic usage
python visualize_mutations_python.py --input mutations.tsv --output figure.png

# With custom cutoff
python visualize_mutations_python.py --input data.tsv --output plot.pdf --cutoff 0.9

# Custom figure size
python visualize_mutations_python.py --input data.tsv --output plot.png --width 20 --height 12
```

## What You Get

1. **Main visualization** (PNG or PDF) - identical to R version
2. **Data table** (`filename_mutations_table.tsv`) - all mutations in spreadsheet format

## Installation Test

Run this to make sure everything works:
```bash
python visualize_mutations_python.py --help
```

You should see the help message with no errors.

## Troubleshooting

### "No module named pandas"
```bash
pip install pandas numpy matplotlib
```

### "python: command not found"
Try `python3` instead:
```bash
python3 visualize_mutations_python.py --input data.tsv --output plot.png
```

### Permission denied
```bash
chmod +x visualize_mutations_python.py
```

## For Your Collaborators

Tell them to:
1. **Skip all the R installation headaches** 
2. **Install 3 Python packages**: `pip install pandas numpy matplotlib`
3. **Use the Python script**: Replace `Rscript visualize_mutations.R` with `python visualize_mutations_python.py`

## Comparison

| Feature | R Version | Python Version |
|---------|-----------|----------------|
| Package installation | 💀 Complex, often fails | ✅ Simple, 3 packages |
| M2 Mac compatibility | ❌ Many issues | ✅ Works perfectly |
| Dependencies | 10+ R packages | 3 Python packages |
| Visual output | Good | Identical |
| For bench scientists | Difficult | Easy |

## Migration from R Version

**Instead of:**
```bash
Rscript visualize_mutations.R --input data.tsv --output plot.pdf --cutoff 0.05
```

**Use:**
```bash
python visualize_mutations_python.py --input data.tsv --output plot.png --cutoff 0.05
```

The output will be visually identical but much more reliable!

---

**Bottom line:** The Python version eliminates all the R package installation headaches while producing identical results. Perfect for bench scientists who just want their analyses to work.