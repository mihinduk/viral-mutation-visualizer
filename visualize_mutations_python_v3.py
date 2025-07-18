#\!/usr/bin/env python3
"""
Viral Mutation Visualizer with Extensible Virus Database
Supports automatic discovery of new virus families
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import sys
import os
from virus_config_framework import get_virus_info

def check_dependencies():
    """Check if all required packages are available"""
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
        return False
    
    print("✅ All required packages found:")
    for pkg, version in versions.items():
        print(f"   - {pkg}: {version}")
    
    return True

def map_position_to_gene(position, virus_info):
    """Map a genomic position to the corresponding gene"""
    gene_coords = virus_info.get('gene_coords', {})
    
    for gene, coords in gene_coords.items():
        if isinstance(coords, list) and len(coords) == 2:
            start, end = coords
            if start <= position <= end:
                return gene
    return None

def setup_gene_visualization(virus_info):
    """Setup gene coordinates and colors for visualization"""
    gene_coords = virus_info.get('gene_coords', {})
    family = virus_info.get('family', 'unknown')
    
    # Load family-specific styling
    from virus_config_framework import virus_manager
    styling = virus_manager.get_family_styling(family)
    
    # Get colors, with defaults for missing genes
    colors = styling.get('colors', {})
    default_colors = ['#4575b4', '#74add1', '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026', '#762a83']
    
    gene_colors = {}
    for i, gene in enumerate(gene_coords.keys()):
        if gene in colors:
            gene_colors[gene] = colors[gene]
        else:
            gene_colors[gene] = default_colors[i % len(default_colors)]
    
    # Get structural/non-structural classification
    structural_genes = styling.get('structural_genes', [])
    nonstructural_genes = styling.get('nonstructural_genes', [])
    
    # If not defined, make educated guess
    if not structural_genes and not nonstructural_genes:
        for gene in gene_coords.keys():
            if gene.lower() in ['c', 'prm', 'env', 'e', 'e1', 'e2', 'e3', 's', 'm', 'n']:
                structural_genes.append(gene)
            else:
                nonstructural_genes.append(gene)
    
    return gene_coords, gene_colors, structural_genes, nonstructural_genes

def create_visualization(df, accession, args):
    """Create the mutation visualization"""
    # Get virus information
    virus_info = get_virus_info(accession)
    virus_name = virus_info['name']
    
    # Setup gene visualization
    gene_coords, gene_colors, structural_genes, nonstructural_genes = setup_gene_visualization(virus_info)
    
    if not gene_coords:
        print(f"⚠️  Warning: No gene coordinates available for {virus_name}")
        print(f"   Consider running: python3 add_virus.py {accession}")
    
    # Map mutations to genes
    df['Gene'] = df['POS'].apply(lambda pos: map_position_to_gene(pos, virus_info))
    
    # Filter by allele frequency cutoff
    filtered_df = df[df['AF'] >= args.cutoff].copy()
    
    print(f"📖 Reading input data...")
    print(f"   Found {len(df)} total mutations")
    print(f"   {len(filtered_df)} mutations above cutoff {args.cutoff}")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(args.width, args.height))
    
    # Determine genome length
    genome_length = virus_info.get('genome_length')
    if not genome_length and gene_coords:
        genome_length = max([coords[1] for coords in gene_coords.values() if isinstance(coords, list)])
    if not genome_length:
        genome_length = df['POS'].max() + 1000 if not df.empty else 11000
    
    # Draw gene structure
    y_gene = 0.5
    gene_height = 0.3
    
    if gene_coords:
        # Draw structural vs non-structural regions
        if structural_genes:
            struct_start = min([gene_coords[g][0] for g in structural_genes if g in gene_coords])
            struct_end = max([gene_coords[g][1] for g in structural_genes if g in gene_coords])
            ax.add_patch(plt.Rectangle((0, y_gene-gene_height/2), struct_end, gene_height, 
                                     facecolor='lightblue', alpha=0.3, zorder=0))
            ax.text(struct_end/2, y_gene+gene_height, 'Structural proteins', 
                   ha='center', va='bottom', fontweight='bold', color='blue')
        
        if nonstructural_genes:
            nonstruct_start = min([gene_coords[g][0] for g in nonstructural_genes if g in gene_coords])
            nonstruct_end = max([gene_coords[g][1] for g in nonstructural_genes if g in gene_coords])
            ax.add_patch(plt.Rectangle((nonstruct_start, y_gene-gene_height/2), 
                                     nonstruct_end-nonstruct_start, gene_height,
                                     facecolor='lightcoral', alpha=0.3, zorder=0))
            ax.text((nonstruct_start+nonstruct_end)/2, y_gene+gene_height, 'Non-structural proteins',
                   ha='center', va='bottom', fontweight='bold', color='red')
        
        # Draw individual genes
        for gene, coords in gene_coords.items():
            if isinstance(coords, list) and len(coords) == 2:
                start, end = coords
                color = gene_colors.get(gene, '#cccccc')
                
                # Draw gene rectangle
                ax.add_patch(plt.Rectangle((start, y_gene-gene_height/2), end-start, gene_height,
                                         facecolor=color, edgecolor='black', linewidth=0.5))
                
                # Add gene label
                ax.text((start+end)/2, y_gene, gene, ha='center', va='center', 
                       fontweight='bold', fontsize=10)
    
    # Plot mutations
    if not filtered_df.empty:
        for _, mut in filtered_df.iterrows():
            pos = mut['POS']
            af = mut['AF']
            
            # Plot mutation line
            ax.plot([pos, pos], [y_gene-gene_height/2-0.1, y_gene-gene_height/2-0.3], 
                   'k-', linewidth=2, alpha=0.7)
            
            # Plot mutation point
            ax.plot(pos, y_gene-gene_height/2-0.3, 'ro', markersize=8*af+2, alpha=0.8)
    
    # Set up the plot
    ax.set_xlim(0, genome_length)
    ax.set_ylim(-0.5, 1.2)
    ax.set_xlabel('Genome Position')
    ax.set_ylabel('')
    
    # Create title with gene info
    if args.mutation_genes and args.mutation_genes != "all":
        gene_desc = f" ({args.mutation_genes} genes)"
    else:
        gene_desc = ""
    
    title = f"Mutations in {accession} - {virus_name}, complete genome{gene_desc} (cutoff: {args.cutoff*100:.1f}%)")
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Remove y-axis
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Add mutation table
    if not filtered_df.empty:
        # Group by gene for table
        gene_mutations = {}
        for _, mut in filtered_df.iterrows():
            gene = mut.get('Gene', 'Unknown')
            if gene not in gene_mutations:
                gene_mutations[gene] = []
            
            gene_mutations[gene].append({
                'Position': f"{mut['POS']}{mut['REF']}>{mut['ALT']}",
                'AA_Change': mut.get('HGVSP', '').replace('p.', '') if 'HGVSP' in mut else 'N/A'
            })
        
        # Create tables for each gene with mutations
        table_y = -0.4
        col_width = 0.15
        
        for gene, mutations in gene_mutations.items():
            if gene and gene \!= 'Unknown':
                # Table header
                ax.text(0.5, table_y, gene, transform=ax.transAxes, fontweight='bold', 
                       ha='center', va='top', fontsize=12, 
                       bbox=dict(boxstyle='round,pad=0.3', facecolor=gene_colors.get(gene, '#cccccc')))
                
                # Table content
                table_data = []
                headers = ['Position', 'AA Change']
                
                for i, mut in enumerate(mutations[:5]):  # Limit to 5 mutations per gene
                    table_data.append([mut['Position'], mut['AA_Change']])
                
                if table_data:
                    table = ax.table(cellText=table_data, colLabels=headers,
                                   cellLoc='center', loc='center',
                                   bbox=[0.1, table_y-0.15, 0.35, 0.1])
                    table.auto_set_font_size(False)
                    table.set_fontsize(8)
                    table.scale(1, 1.5)
    
    plt.tight_layout()
    return fig

def main():
    parser = argparse.ArgumentParser(
        description='Visualize viral genome mutations with extensible virus database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --input mutations.tsv --output plot.png --cutoff 0.05
  %(prog)s --input data.tsv --output figure.pdf --accession NC_075022.1 --cutoff 0.001
"""
    )
    
    parser.add_argument('--input', required=True, help='Input TSV file with SnpEff-annotated output')
    parser.add_argument('--output', required=True, help='Output file path (PNG or PDF)')
    parser.add_argument('--cutoff', type=float, default=0.05, help='Allele frequency cutoff (default: 0.05)')
    parser.add_argument('--mutation-genes', default='all', help='Genes to display: "all", "structural", "non-structural", or comma-separated names')
    parser.add_argument('--width', type=float, default=16, help='Figure width in inches (default: 16)')
    parser.add_argument('--height', type=float, default=10, help='Figure height in inches (default: 10)')
    parser.add_argument('--highlight-freq', type=float, default=0.5, help='Highlight mutations above this frequency (default: 0.5)')
    parser.add_argument('--accession', help='Manually specify GenBank accession (overrides accession in data)')
    
    args = parser.parse_args()
    
    print("Viral Mutation Visualizer (Extensible Database Version)")
    print("=" * 55)
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    print(f"📁 Input file: {args.input}")
    print(f"🎯 Output file: {args.output}")
    print(f"📊 Allele frequency cutoff: {args.cutoff}")
    
    # Read input data
    try:
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        print(f"❌ Error reading input file: {e}")
        sys.exit(1)
    
    # Get accession
    if args.accession:
        accession = args.accession
    elif 'CHROM' in df.columns:
        accession = df['CHROM'].iloc[0]
    else:
        print("❌ No accession found. Please specify with --accession")
        sys.exit(1)
    
    print(f"🧬 Processing accession: {accession}")
    
    # Create visualization
    print("🎨 Creating visualization...")
    fig = create_visualization(df, accession, args)
    
    # Save figure
    print("💾 Saving figure...")
    fig.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save filtered data table
    filtered_df = df[df['AF'] >= args.cutoff] if 'AF' in df.columns else df
    table_output = args.output.replace('.png', '_mutations_table.tsv').replace('.pdf', '_mutations_table.tsv')
    filtered_df.to_csv(table_output, sep='\t', index=False)
    
    print(f"✅ Figure saved to: {args.output}")
    print(f"✅ Mutations table saved to: {table_output}")
    print("\n🎉 Visualization complete\!")
    print(f"   Main figure: {args.output}")
    print(f"   Data table: {table_output}")

if __name__ == '__main__':
    main()
