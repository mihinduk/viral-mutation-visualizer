# Easy Setup Instructions for Bench Scientists

## Quick Fix for R Package Errors

If you're getting errors like "There is no package named ggplot2" or "cowplot had a non-zero status", follow these simple steps:

### Step 1: Run the Automatic Fix (EASIEST)

**On Windows:**
1. Double-click the file `fix_R_packages.R` 
2. Wait 5-10 minutes for installation to complete
3. You're done!

**On Mac/Linux:**
1. Open Terminal
2. Navigate to the viral-mutation-visualizer folder
3. Run: `Rscript fix_R_packages.R`
4. Wait 5-10 minutes for installation to complete

### Step 2: Test That It Worked

Try running your visualization command again:
```bash
Rscript visualize_mutations.R --input your_data.tsv --output plot.pdf
```

## What This Fix Does

- ✅ Removes the problematic `cowplot` package (it wasn't being used anyway)
- ✅ Automatically installs all required R packages
- ✅ Tests that everything works correctly
- ✅ Gives you a clear success/failure message

## If You Still Have Problems

1. **Make sure R is installed**: Type `R --version` in terminal/command prompt
2. **Check internet connection**: Package installation requires internet
3. **Try running the fix again**: Sometimes packages need multiple attempts
4. **Contact your bioinformatics support**: Show them the error messages

## Alternative: Use RStudio (If You Have It)

1. Open RStudio
2. Copy and paste this into the console:

```r
# Install all packages at once
install.packages(c("optparse", "ggplot2", "dplyr", "tidyr", "readr", "stringr", 
                   "patchwork", "ggrepel", "RColorBrewer", "gridExtra"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("rentrez", "seqinr"))
```

3. Press Enter and wait for installation to complete

## Success Message

When everything works, you should see:
```
🎉 SUCCESS! All packages installed correctly.
You can now run the viral mutation visualization scripts.
```

---

**Need help?** This fix is designed to work automatically. If it doesn't work on your system, contact your IT or bioinformatics support and show them this file.