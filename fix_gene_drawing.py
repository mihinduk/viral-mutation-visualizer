#\!/usr/bin/env python3
import re

# Read the file
with open('visualize_mutations_python_outdir.py', 'r') as f:
    content = f.read()

# Fix the UTR drawing to be more robust
# Replace the 5' UTR section
old_utr5 = # Draw 5 UTR
    utr5_rect = Rectangle((0, gene_y), gene_coords[C][0]-1, gene_height, 
                         facecolor=lightgray, edgecolor=black, linewidth=0.5)
    ax.add_patch(utr5_rect)"""

new_utr5 = """# Draw 5 UTR
