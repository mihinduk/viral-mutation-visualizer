#\!/usr/bin/env python3
import re

# Read the file
with open('visualize_mutations_python_outdir.py', 'r') as f:
    lines = f.readlines()

# Find and replace lines that use GENE_COORDS with gene_coords
for i in range(len(lines)):
    # Replace GENE_COORDS with gene_coords (from dynamic config)
    if 'GENE_COORDS[' in lines[i]:
        lines[i] = lines[i].replace('GENE_COORDS', 'gene_coords')
    
    # Fix the gene drawing loop to use gene_coords instead of GENE_COORDS
    if 'for gene, (start, end) in GENE_COORDS.items():' in lines[i]:
        lines[i] = '    for gene, (start, end) in gene_coords.items():\n'

# Write the updated file
with open('visualize_mutations_python_outdir.py', 'w') as f:
    f.writelines(lines)

print('Fixed gene coordinate references')
