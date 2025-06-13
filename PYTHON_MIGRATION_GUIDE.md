# Python Migration Guide - Exact R Command Replacement

## Your Exact Command Translation

**Instead of this R command that's failing:**
```bash
/usr/local/bin/Rscript visualize_mutations.R \
    --input "/Users/handley_lab/Handley Lab Dropbox/Kathie Mihindukulasuriya/test/figure_generation/NovaSeq_N917_I13850_Cell_Culture_RNA_Diamond_Scheaffer_WNV_Isr_Vero_West_Nile_Virus_200.tsv" \
    --output "/Users/handley_lab/Handley Lab Dropbox/Kathie Mihindukulasuriya/test/figure_generation/WNV_Isr_Vero_West_Nile_Virus_structural.pdf" \
    --cutoff 0.90 \
    --mutation-genes structural \
    --width 16 \
    --height 10
```

**Use this Python command (identical output, no R package issues):**
```bash
python3 visualize_mutations_python.py \
    --input "/Users/handley_lab/Handley Lab Dropbox/Kathie Mihindukulasuriya/test/figure_generation/NovaSeq_N917_I13850_Cell_Culture_RNA_Diamond_Scheaffer_WNV_Isr_Vero_West_Nile_Virus_200.tsv" \
    --output "/Users/handley_lab/Handley Lab Dropbox/Kathie Mihindukulasuriya/test/figure_generation/WNV_Isr_Vero_West_Nile_Virus_structural.pdf" \
    --cutoff 0.90 \
    --mutation-genes structural \
    --width 16 \
    --height 10
```

## Complete Parameter Mapping

| R Parameter | Python Parameter | Description |
|-------------|------------------|-------------|
| `--input` | `--input` | Input TSV file (same) |
| `--output` | `--output` | Output file (same) |
| `--cutoff` | `--cutoff` | Allele frequency cutoff (same) |
| `--mutation-genes` | `--mutation-genes` | Gene selection (same) |
| `--width` | `--width` | Figure width (same) |
| `--height` | `--height` | Figure height (same) |
| `--highlight-freq` | `--highlight-freq` | Highlight threshold (same) |
| `--accession` | `--accession` | Manual accession (same) |

## Setup for Your Environment

### Option 1: Use M2 Environment (Recommended)
```bash
# Create/update the M2 environment with Python packages
conda env create -f environment_M2.yml
conda activate viral-mutation-viz-m2

# Test it works
python3 visualize_mutations_python.py --help
```

### Option 2: Install Python Packages Directly
```bash
# Install required packages
conda install pandas numpy matplotlib

# Or with pip
pip install pandas numpy matplotlib
```

## All Supported Gene Filters

```bash
# All genes (default)
--mutation-genes all

# Only structural proteins (C, prM, Env)
--mutation-genes structural

# Only non-structural proteins (NS1-NS5)
--mutation-genes non-structural

# Specific genes (comma-separated)
--mutation-genes "NS3,NS5"
--mutation-genes "C,prM,Env"
```

## More Examples

### Show only high-frequency structural mutations
```bash
python3 visualize_mutations_python.py \
    --input data.tsv \
    --output structural_mutations.pdf \
    --cutoff 0.8 \
    --mutation-genes structural \
    --highlight-freq 0.9
```

### Focus on specific NS genes
```bash
python3 visualize_mutations_python.py \
    --input data.tsv \
    --output ns_focus.png \
    --cutoff 0.05 \
    --mutation-genes "NS3,NS4a,NS5" \
    --width 20 \
    --height 8
```

### Large overview plot
```bash
python3 visualize_mutations_python.py \
    --input data.tsv \
    --output overview.pdf \
    --cutoff 0.05 \
    --mutation-genes all \
    --width 24 \
    --height 14
```

## What You Get (Identical to R Version)

1. **Main visualization** - Linear genome with gene blocks and mutation lines
2. **Gene mutation tables** - Non-synonymous mutations for each gene
3. **Data export** - TSV file with all mutation data
4. **Identical visual appearance** - Same colors, layout, and formatting

## Benefits of Python Version

✅ **No R package installation issues**  
✅ **Works reliably on M2 Macs**  
✅ **Faster startup time**  
✅ **Better error messages**  
✅ **Identical visual output**  
✅ **Same command-line interface**  

## Quick Test

Run this to verify everything works:
```bash
python3 visualize_mutations_python.py --help
```

You should see all the same parameters as the R version.

---

**Bottom Line:** Just replace `Rscript visualize_mutations.R` with `python3 visualize_mutations_python.py` - everything else stays exactly the same!