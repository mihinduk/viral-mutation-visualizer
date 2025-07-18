#!/usr/bin/env python3
"""
visualize_mutations_python.py
Python version of the viral mutation visualization tool.
Creates identical output to the R version but uses only standard Python packages.
Designed for M2 Mac compatibility and ease of use for bench scientists.
"""

import argparse
# Use dynamic virus configuration if available
try:
    from virus_config_framework import VirusConfigManager
    vcm = VirusConfigManager()
    # Test if it can actually load data properly
    test_info = vcm.get_virus_info('NC_075022.1')
    if not test_info.get('colors') or len(test_info.get('gene_coords', {})) < 5:
        # Framework is not working properly, fall back to simple manager
        raise ImportError('Framework not loading proper data')
    
    def get_virus_name(accession):
        info = vcm.get_virus_info(accession)
        return info.get("name", f"Unknown virus ({accession})")
    
    def get_gene_info(accession):
        info = vcm.get_virus_info(accession)
        return (info.get("gene_coords", {}), 
                info.get("colors", {}),
                info.get("structural_genes", []),
                info.get("nonstructural_genes", []))
    
    print("Using dynamic virus configuration framework")
except ImportError:
    try:
        from virus_config_simple import simple_virus_manager
        
        def get_virus_name(accession):
            info = simple_virus_manager.get_virus_info(accession)
            return info.get("name", f"Unknown virus ({accession})")
        
        def get_gene_info(accession):
            info = simple_virus_manager.get_virus_info(accession)
            return (info.get("gene_coords", {}), 
                    info.get("colors", {}),
                    info.get("structural_genes", []),
                    info.get("nonstructural_genes", []))
        
        print("Using simple virus configuration manager")
    except ImportError:
        from quick_virus_fix import get_virus_name, get_gene_info
        print("Using quick fix virus configuration")

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
import warnings

# Suppress matplotlib warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)

def check_dependencies():
    """Check if all required Python packages are available"""
    required_packages = {
        'pandas': 'pandas',
        'numpy': 'numpy', 
        'matplotlib': 'matplotlib'
    }
    
    missing = []
    versions = {}
    
    for package_name, import_name in required_packages.items():
        try:
            module = __import__(import_name)
            if hasattr(module, '__version__'):
                versions[package_name] = module.__version__
            else:
                versions[package_name] = 'unknown'
        except ImportError:
            missing.append(package_name)
    
    if missing:
        print("❌ Missing required packages:")
        for pkg in missing:
            print(f"   - {pkg}")
        print("\nTo install missing packages:")
        print("   pip install " + " ".join(missing))
        print("\nOr if using conda:")
        print("   conda install " + " ".join(missing))
        return False
    
    print("✅ All required packages found:")
    for pkg, version in versions.items():
        print(f"   - {pkg}: {version}")
    
    return True

# Gene coordinates for West Nile virus (1-based)
GENE_COORDS = {
    'C': (97, 465),
    'prM': (466, 966), 
    'Env': (967, 2469),
    'NS1': (2470, 3525),
    'NS2a': (3526, 4218),
    'NS2b': (4219, 4611),
    'NS3': (4612, 6468),
    'NS4a': (6469, 6915),
    'NS4b': (6916, 7671),
    'NS5': (7672, 10395)
}

# Colors matching the R version exactly
GENE_COLORS = {
    'C': '#4575b4',
    'prM': '#74add1', 
    'Env': '#abd9e9',
    'NS1': '#fdae61',
    'NS2a': '#f46d43',
    'NS2b': '#d73027',
    'NS3': '#a50026',
    'NS4a': '#762a83',
    'NS4b': '#9970ab',
    'NS5': '#c2a5cf'
}

# Structural vs non-structural genes
structural_genes = ['C', 'prM', 'Env']
NONstructural_genes = ['NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5']

# Genome length for West Nile virus
GENOME_LENGTH = 11029

def map_position_to_gene(position, accession):
    """Map a genomic position to the corresponding gene"""
    gene_coords, _, _, _ = get_gene_info(accession)
    for gene, (start, end) in gene_coords.items():
        if start <= position <= end:
            return gene
    return None

def parse_amino_acid_change(hgvs_p):
    """Parse amino acid change from HGVSp notation"""
    if pd.isna(hgvs_p) or hgvs_p == '':
        return None
    
    # Handle various HGVSp formats
    if 'p.' in str(hgvs_p):
        aa_change = str(hgvs_p).split('p.')[1]
        # Convert * to TER for stop codons
        aa_change = aa_change.replace('*', 'TER')
        return aa_change
    
    return None

def filter_mutations(df, cutoff, effect_filter=None):
    """Filter mutations by allele frequency and optionally by effect"""
    # Convert allele frequency to numeric if it's not already
    if 'Allele_Frequency' in df.columns:
        df['Allele_Frequency'] = pd.to_numeric(df['Allele_Frequency'], errors='coerce')
        df = df[df['Allele_Frequency'] >= cutoff]
    
    # Filter by effect type if specified
    if effect_filter:
        if 'EFFECT' in df.columns:
            df = df[df['EFFECT'].isin(effect_filter)]
    
    return df

def create_genome_diagram(ax, mutations_df, title, gene_filter="all", highlight_freq=0.5, accession=None):
    """Create the linear genome diagram with gene blocks"""
    # Get virus-specific gene information
    gene_coords, gene_colors, _, _ = get_gene_info(accession)

    """Create the linear genome diagram with gene blocks"""
    # Get dynamic gene info for this virus
    gene_coords, gene_colors, structural_genes, nonstructural_genes = get_gene_info(accession)    
    # If we only have Polyprotein, use hardcoded coordinates for known viruses
    if len(gene_coords) == 1 and 'Polyprotein' in gene_coords:
        if accession == 'HM440560.1' or accession == 'HM440560':
            # Hardcoded Powassan virus gene coordinates
            gene_coords = {
                'C': [110, 357],
                'prM': [358, 858],
                'Env': [859, 2361],
                'NS1': [2362, 3417],
                'NS2a': [3418, 4110],
                'NS2b': [4111, 4503],
                'NS3': [4504, 6360],
                'NS4a': [6361, 6807],
                'NS4b': [6808, 7563],
                'NS5': [7564, 10287]
            }
            # Use standard flavivirus colors
            gene_colors = {
                'C': '#4575b4',
                'prM': '#74add1',
                'Env': '#abd9e9',
                'NS1': '#fee090',
                'NS2a': '#fdae61',
                'NS2b': '#f46d43',
                'NS3': '#d73027',
                'NS4a': '#a50026',
                'NS4b': '#762a83',
                'NS5': '#5aae61'
            }
            structural_genes = ['C', 'prM', 'Env']
            nonstructural_genes = ['NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5']
            print('Using hardcoded Powassan virus gene coordinates')

        elif accession == 'NC_075022.1':
            # Hardcoded VEEV (alphavirus) gene coordinates
            gene_coords = {
                'nsP1': [84, 1685],
                'nsP2': [1686, 4226],
                'nsP3': [4227, 5945],
                'nsP4': [5946, 7664],
                'C': [7667, 7900],
                'E3': [7901, 8059],
                'E2': [8060, 9394],
                'E1': [9395, 11023]
            }
            # Alphavirus-specific colors
            gene_colors = {
                'nsP1': '#d62728',    # Red
                'nsP2': '#ff7f0e',    # Orange  
                'nsP3': '#2ca02c',    # Green
                'nsP4': '#1f77b4',    # Blue
                'C': '#9467bd',       # Purple
                'E3': '#8c564b',      # Brown
                'E2': '#e377c2',      # Pink
                'E1': '#7f7f7f'       # Gray
            }
            structural_genes = ['C', 'E3', 'E2', 'E1']
            nonstructural_genes = ['nsP1', 'nsP2', 'nsP3', 'nsP4']
            print('Using hardcoded VEEV gene coordinates')
    
    
    # Set up the plot
    ax.set_xlim(0, GENOME_LENGTH)
    ax.set_ylim(-0.8, 1.2)
    
    # Draw gene blocks
    gene_height = 0.3
    gene_y = 0.35
    
    # Draw 5' UTR
    # Find first gene start
    first_gene_start = min(coord[0] for coord in gene_coords.values()) if gene_coords else 1
    if first_gene_start > 1:
        utr5_rect = Rectangle((0, gene_y), first_gene_start-1, gene_height, 
                         facecolor='lightgray', edgecolor='black', linewidth=0.5)
        ax.add_patch(utr5_rect)
    
    # Draw genes
    for gene, (start, end) in gene_coords.items():
        width = end - start + 1
        rect = Rectangle((start-1, gene_y), width, gene_height,
                        facecolor=gene_colors.get(gene, "#808080"), edgecolor='black', linewidth=0.5)
        ax.add_patch(rect)
        
        # Add gene label
        gene_center = start + width/2 - 1
        ax.text(gene_center, gene_y + gene_height/2, gene, 
               ha='center', va='center', fontsize=9, fontweight='bold')
    
    # Draw 3' UTR
    # Find last gene end
    last_gene_end = max(coord[1] for coord in gene_coords.values()) if gene_coords else GENOME_LENGTH
    utr3_start = last_gene_end
    if utr3_start < GENOME_LENGTH:
        utr3_rect = Rectangle((utr3_start, gene_y), GENOME_LENGTH - utr3_start, gene_height,
                             facecolor='lightgray', edgecolor='black', linewidth=0.5)
        ax.add_patch(utr3_rect)
    
    # Add structural/non-structural labels
    # Find structural gene boundaries
    struct_coords = [(start, end) for gene, (start, end) in gene_coords.items() if gene in structural_genes]
    if struct_coords:
        struct_start = min(coord[0] for coord in struct_coords) - 1
        struct_end = max(coord[1] for coord in struct_coords) - 1
    else:
        struct_start = 0
        struct_end = 0
    struct_center = (struct_start + struct_end) / 2
    
    # Find non-structural gene boundaries
    nonstruct_coords = [(start, end) for gene, (start, end) in gene_coords.items() if gene in nonstructural_genes]
    if nonstruct_coords:
        nonstruct_start = min(coord[0] for coord in nonstruct_coords) - 1
        nonstruct_end = max(coord[1] for coord in nonstruct_coords) - 1
    else:
        nonstruct_start = 0
        nonstruct_end = 0  
    nonstruct_center = (nonstruct_start + nonstruct_end) / 2
    
    # Structural proteins label and underline
    ax.text(struct_center, 0.8, 'Structural proteins', ha='center', va='center', 
           fontsize=12, fontweight='bold', color='#4575b4')
    ax.plot([struct_start, struct_end], [0.75, 0.75], color='#4575b4', linewidth=3)
    
    # Non-structural proteins label and underline  
    ax.text(nonstruct_center, 0.8, 'Non-structural proteins', ha='center', va='center',
           fontsize=12, fontweight='bold', color='#d73027')
    ax.plot([nonstruct_start, nonstruct_end], [0.75, 0.75], color='#d73027', linewidth=3)
    
    # Draw mutation lines
    if not mutations_df.empty:
        # Apply gene filtering to mutations for display
        display_mutations = filter_genes_for_display(mutations_df.copy(), gene_filter)
        
        for _, mutation in display_mutations.iterrows():
            pos = mutation['POS']
            gene = map_position_to_gene(pos, accession)
            if gene:
                # Determine color based on mutation type
                mutation_type = mutation.get("EFFECT", "")
                if "synonymous_variant" in mutation_type:
                    color = "#2ca02c"  # Green for synonymous
                elif "missense_variant" in mutation_type:
                    color = "#ff7f0e"  # Orange for missense
                elif "stop_gained" in mutation_type or "nonsense" in mutation_type:
                    color = "#d62728"  # Red for stop/nonsense
                else:
                    color = gene_colors.get(gene, "#1f77b4")  # Blue default
                # Highlight high-frequency mutations with thicker lines
                freq = mutation.get('Allele_Frequency', 0)
                if freq >= highlight_freq:
                    linewidth = 3.0
                    alpha = 1.0
                else:
                    linewidth = 1.5
                    alpha = 0.8
                ax.plot([pos, pos], [0.25, 0.35], color=color, linewidth=linewidth, alpha=alpha)
    
    # Format x-axis
    ax.set_xlabel('Genome Position', fontsize=12, fontweight='bold')
    ax.set_xticks(np.arange(0, GENOME_LENGTH+1, 1000))
    ax.set_xticklabels([f'{x:,}' for x in np.arange(0, GENOME_LENGTH+1, 1000)])
    
    # Remove y-axis
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Add title
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

def filter_genes_for_display(mutations_df, gene_filter):
    """Filter mutations based on gene selection"""
    if gene_filter == 'all':
        return mutations_df
    elif gene_filter == 'structural':
        # Filter for structural genes
        mutations_df['Gene'] = mutations_df['POS'].apply(lambda pos: map_position_to_gene(pos, accession))
        return mutations_df[mutations_df['Gene'].isin(structural_genes)]
    elif gene_filter == 'non-structural':
        # Filter for non-structural genes  
        mutations_df['Gene'] = mutations_df['POS'].apply(lambda pos: map_position_to_gene(pos, accession))
        return mutations_df[mutations_df['Gene'].isin(NONstructural_genes)]
    else:
        # Comma-separated gene list
        selected_genes = [g.strip() for g in gene_filter.split(',')]
        mutations_df['Gene'] = mutations_df['POS'].apply(lambda pos: map_position_to_gene(pos, accession))
        return mutations_df[mutations_df['Gene'].isin(selected_genes)]

def create_mutation_tables(fig, mutations_df, start_row=0.4, gene_filter="all", accession=None):
    """Create tables showing non-synonymous mutations for each gene"""
    
    # Filter for non-synonymous mutations only
    nonsyn_effects = ['missense_variant', 'nonsense_variant', 'stop_gained', 'stop_lost']
    nonsyn_mutations = mutations_df[mutations_df['EFFECT'].isin(nonsyn_effects)].copy()
    
    # Apply gene filtering
    nonsyn_mutations = filter_genes_for_display(nonsyn_mutations, gene_filter)
    
    if nonsyn_mutations.empty:
        # Add text saying no non-synonymous mutations found
        fig.text(0.5, 0.3, 'No non-synonymous mutations found above the specified cutoff',
                ha='center', va='center', fontsize=14, style='italic')
        return
    
    # Add gene column based on position
    nonsyn_mutations['Gene'] = nonsyn_mutations['POS'].apply(lambda pos: map_position_to_gene(pos, accession))
    nonsyn_mutations = nonsyn_mutations.dropna(subset=['Gene'])
    
    # Parse amino acid changes
    nonsyn_mutations['AA_Change'] = nonsyn_mutations['HGVSp'].apply(parse_amino_acid_change)
    
    # Group by gene
    # Get dynamic gene list for this virus
    gene_coords, gene_colors, _, _ = get_gene_info(accession)
    genes_with_mutations = []
    for gene in gene_coords.keys():
        gene_mutations = nonsyn_mutations[nonsyn_mutations['Gene'] == gene]
        if not gene_mutations.empty:
            genes_with_mutations.append((gene, gene_mutations))
    
    if not genes_with_mutations:
        fig.text(0.5, 0.3, 'No non-synonymous mutations found above the specified cutoff',
                ha='center', va='center', fontsize=14, style='italic')
        return
    
    # Create tables in a grid layout (max 5 columns)
    max_cols = 5
    n_genes = len(genes_with_mutations)
    n_rows = (n_genes + max_cols - 1) // max_cols
    
    # Calculate table positions
    table_width = 0.18
    table_spacing = 0.02
    total_width = min(n_genes, max_cols) * table_width + (min(n_genes, max_cols) - 1) * table_spacing
    start_x = (1 - total_width) / 2
    
    for i, (gene, gene_mutations) in enumerate(genes_with_mutations):
        row = i // max_cols
        col = i % max_cols
        
        # Calculate position
        x = start_x + col * (table_width + table_spacing)
        y = start_row - row * 0.25
        
        # Create table data
        table_data = []
        for _, mut in gene_mutations.iterrows():
            nt_change = f"{mut['POS']}{mut['REF']}>{mut['ALT']}"
            aa_change = mut['AA_Change'] if mut['AA_Change'] else 'Unknown'
            table_data.append([nt_change, aa_change])
        
        # Create table
        if table_data:
            # Create subplot for this table
            table_ax = fig.add_axes([x, y-0.15, table_width, 0.15])
            table_ax.axis('off')
            
            # Create table
            table = table_ax.table(cellText=table_data,
                                 colLabels=['Position', 'AA Change'],
                                 cellLoc='center',
                                 loc='center',
                                 colWidths=[0.5, 0.5])
            
            # Style the table
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 1.5)
            
            # Color the header
            for i in range(2):
                table[(0, i)].set_facecolor(gene_colors.get(gene, "#808080"))
                table[(0, i)].set_text_props(weight='bold', color='white')
            
            # Add gene name above table
            fig.text(x + table_width/2, y + 0.02, gene, 
                    ha='center', va='center', fontsize=12, 
                    fontweight='bold', color=gene_colors.get(gene, "#808080"))

def save_mutations_table(mutations_df, output_path, cutoff, accession=None):
    """Save all mutations to a TSV file"""
    if mutations_df.empty:
        return
    
    # Prepare output data
    output_data = []
    for _, mut in mutations_df.iterrows():
        gene = map_position_to_gene(mut['POS'], accession)
        aa_change = parse_amino_acid_change(mut.get('HGVSp', ''))
        
        freq_percent = f"{mut['Allele_Frequency']*100:.1f}%" if 'Allele_Frequency' in mut else "Unknown"
        
        output_data.append({
            'Position': mut['POS'],
            'Gene': gene if gene else 'Unknown',
            'Nucleotide_Change': f"{mut['REF']}>{mut['ALT']}",
            'Amino_Acid_Change': aa_change if aa_change else 'N/A',
            'Effect': mut.get('EFFECT', 'Unknown'),
            'Impact': mut.get('PUTATIVE_IMPACT', 'Unknown'),
            'Frequency_Percent': freq_percent
        })
    
    # Create DataFrame and save
    output_df = pd.DataFrame(output_data)
    output_df.to_csv(output_path, sep='\t', index=False)
    print(f"✅ Mutations table saved to: {output_path}")


def generate_html_report(mutations_df, image_path, html_path, accession, cutoff):
    """Generate HTML report with embedded image and mutation data"""
    import os
    
    # Get virus info for the report
    virus_name = get_virus_name(accession)
    
    # Count mutations
    total_mutations = len(mutations_df)
    
    # Create HTML content
    html_content = f"""<\!DOCTYPE html>
<html>
<head>
    <title>Viral Mutation Report - {virus_name}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background-color: #f8f9fa; padding: 20px; border-radius: 5px; margin-bottom: 20px; }}
        .image-container {{ text-align: center; margin: 20px 0; }}
        .image-container img {{ max-width: 100%; height: auto; border: 1px solid #ddd; }}
        .stats {{ display: flex; justify-content: space-around; margin: 20px 0; }}
        .stat-box {{ background-color: #e9ecef; padding: 15px; border-radius: 5px; text-align: center; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Viral Mutation Analysis Report</h1>
        <h2>{virus_name} ({accession})</h2>
        <p>Generated on: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    
    <div class="stats">
        <div class="stat-box">
            <h3>Total Mutations</h3>
            <p style="font-size: 24px; font-weight: bold;">{total_mutations}</p>
        </div>
        <div class="stat-box">
            <h3>Cutoff Threshold</h3>
            <p style="font-size: 24px; font-weight: bold;">{cutoff*100:.1f}%</p>
        </div>
    </div>
    
    <div class="image-container">
        <h3>Genome Visualization</h3>
        <img src="{os.path.basename(image_path)}" alt="Viral mutation visualization">
    </div>
    
    <div style="margin-top: 40px; font-size: 12px; color: #666;">
        <p>Generated with Viral Mutation Visualizer  < /dev/null |  Family-based configuration system</p>
    </div>
</body>
</html>"""
    
    # Write HTML file
    with open(html_path, 'w') as f:
        f.write(html_content)


def main():
    parser = argparse.ArgumentParser(
        description='Visualize viral genome mutations (Python version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --input mutations.tsv --output plot.png --cutoff 0.05
  %(prog)s --input data.tsv --output figure.pdf --cutoff 0.9 --mutation-genes structural
  %(prog)s --input data.tsv --output plot.png --cutoff 0.05 --mutation-genes "NS3,NS5" --width 20 --height 12
  %(prog)s --input data.tsv --output plot.pdf --cutoff 0.9 --mutation-genes non-structural --highlight-freq 0.8
        """
    )
    
    parser.add_argument('--input', required=True,
                       help='Input TSV file with SnpEff-annotated LoFreq output')
    parser.add_argument('--output', required=True, 
                       help='Output file path (PNG or PDF)')
    parser.add_argument('--cutoff', type=float, default=0.05,
                       help='Allele frequency cutoff (default: 0.05)')
    parser.add_argument('--mutation-genes', default='all',
                       help='Genes to display mutations for: "all", "structural", "non-structural", or comma-separated gene names (default: all)')
    parser.add_argument('--width', type=float, default=16,
                       help='Figure width in inches (default: 16)')
    parser.add_argument('--height', type=float, default=10,
                       help='Figure height in inches (default: 10)')
    parser.add_argument('--highlight-freq', type=float, default=0.5,
                       help='Highlight mutations with frequency above this threshold (default: 0.5)')
    parser.add_argument('--accession', 
                       help='Manually specify GenBank accession (overrides accession in CHROM column)')
    parser.add_argument("--outdir", help="Output directory (default: current directory)")
    
    args = parser.parse_args()
    
    print("Viral Mutation Visualizer (Python Version)")
    print("=" * 50)
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    # Check input file
    if not Path(args.input).exists():
        print(f"❌ Input file not found: {args.input}")
        sys.exit(1)
    
    print(f"📁 Input file: {args.input}")
    print(f"🎯 Output file: {args.output}")
    print(f"📊 Allele frequency cutoff: {args.cutoff}")
    
    try:
        # Read input data
        print("\n📖 Reading input data...")
        df = pd.read_csv(args.input, sep='\t')
        print(f"   Found {len(df)} total mutations")
        
        # Required columns
        required_cols = ['POS', 'REF', 'ALT']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"❌ Missing required columns: {missing_cols}")
            sys.exit(1)
        
        # Filter mutations by cutoff
        filtered_df = filter_mutations(df, args.cutoff)
        print(f"   {len(filtered_df)} mutations above cutoff {args.cutoff}")
        
        # Create the visualization
        print("\n🎨 Creating visualization...")
        fig = plt.figure(figsize=(args.width, args.height))
        
        # Main genome diagram (top 60% of figure)
        genome_ax = fig.add_axes([0.1, 0.5, 0.8, 0.4])
        
        # Get accession from command line or CHROM column
        if args.accession:
            accession = args.accession
        else:
            accession = df['CHROM'].iloc[0] if 'CHROM' in df.columns and not df.empty else "Unknown"
        
        # Create title based on gene filter
        gene_desc = ""
        if args.mutation_genes == 'structural':
            gene_desc = " (structural proteins)"
        elif args.mutation_genes == 'non-structural':
            gene_desc = " (non-structural proteins)"
        elif args.mutation_genes != 'all':
            gene_desc = f" ({args.mutation_genes})"
        
        virus_name = get_virus_name(accession)
        title = f"Mutations in {accession} - {virus_name}, complete genome{gene_desc} (cutoff: {args.cutoff*100:.1f}%)"
        
        create_genome_diagram(genome_ax, filtered_df, title, args.mutation_genes, args.highlight_freq, accession)
        
        # Mutation tables (bottom 40% of figure)
        create_mutation_tables(fig, filtered_df, start_row=0.4, gene_filter=args.mutation_genes, accession=accession)
        
        # Save the figure
        print(f"\n💾 Saving figure...")
        fig.savefig(args.output, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        
        # Handle output directory  
        if args.outdir:
            # Create output directory if specified
            outdir = Path(args.outdir)
            outdir.mkdir(parents=True, exist_ok=True)
            
            # Use output filename but place in specified directory
            output_filename = Path(args.output).name
            final_output_path = outdir / output_filename
            
            # Table path in the same directory
            output_base = Path(args.output).stem
            table_path = outdir / f"{output_base}_mutations_table.tsv"
            
            print(f"\\n📁 Using output directory: {outdir}")
            print(f"   Full path: {outdir.absolute()}")
        else:
            # Use original output path
            final_output_path = Path(args.output)
            output_base = final_output_path.stem
            output_dir = final_output_path.parent
            table_path = output_dir / f"{output_base}_mutations_table.tsv"
        
        # Save the figure
        print(f"\\n💾 Saving figure...")
        fig.savefig(final_output_path, dpi=300, bbox_inches="tight", 
                   facecolor="white", edgecolor="none")
        print(f"✅ Figure saved to: {final_output_path}")
        
        # Save mutations table
        save_mutations_table(filtered_df, table_path, args.cutoff, accession)
        
        plt.close()
        
        print("\\n🎉 Visualization complete\!")
        print(f"   Main figure: {final_output_path}")
        print(f"   Data table: {table_path}")
        
    except Exception as e:
        print(f"\\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
