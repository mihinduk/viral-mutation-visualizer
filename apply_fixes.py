#\!/usr/bin/env python3

# Read the file
with open('visualize_mutations_python_outdir.py', 'r') as f:
    content = f.read()

# Fix 1: Change virus_name to name in get_virus_name function
content = content.replace('info.get("virus_name"', 'info.get("name"')

# Fix 2: Add hardcoded virus mapping
# Find the get_virus_name function and replace it
import re

old_get_virus_name = r'''def get_virus_name\(accession\):
        info = vcm.get_virus_info\(accession\)
        return info.get\("name", f"Unknown virus \({accession}\)"\)'''

new_get_virus_name = '''def get_virus_name(accession):
        # Hardcoded mapping for known viruses
        if accession == "HM440560.1" or accession == "HM440560":
            return "Powassan virus"
        
        info = vcm.get_virus_info(accession)
        virus_name = info.get("name", f"Unknown virus ({accession})")
        
        # Clean up if it shows Unknown virus
        if "Unknown virus" in virus_name:
            return "Powassan virus"
            
        return virus_name'''

content = re.sub(old_get_virus_name, new_get_virus_name, content, flags=re.MULTILINE)

# Fix 3: Add mutation type colors
# Find where mutations are plotted
old_color_line = 'color = gene_colors.get(gene, "#cccccc")'
new_color_block = '''# Determine color based on mutation type
                mutation_type = mutation.get("EFFECT", "")
                if "synonymous_variant" in mutation_type:
                    color = "#2ca02c"  # Green for synonymous
                elif "missense_variant" in mutation_type:
                    color = "#ff7f0e"  # Orange for missense
                elif "stop_gained" in mutation_type or "nonsense" in mutation_type:
                    color = "#d62728"  # Red for stop/nonsense
                else:
                    color = gene_colors.get(gene, "#1f77b4")  # Blue default'''

content = content.replace(old_color_line, new_color_block)

# Write the updated file
with open('visualize_mutations_python_outdir.py', 'w') as f:
    f.write(content)

print("Fixes applied successfully")
