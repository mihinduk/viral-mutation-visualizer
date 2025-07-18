#\!/usr/bin/env python3
"""
visualize_depth.py
Visualizes read depth across viral genomes with gene annotations.
Integrates with the viral mutation visualizer framework for consistent styling.
"""

import argparse
import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import subprocess
import warnings
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

# Suppress matplotlib warnings
warnings.filterwarnings('ignore', category=UserWarning)
# Load known virus configurations
import json
try:
    with open('virus_configs/known_viruses.json', 'r') as f:
        KNOWN_VIRUSES = json.load(f)
except:
    KNOWN_VIRUSES = {}


def check_dependencies():
    """Check if all required tools and packages are available"""
    # Check Python packages
    required_packages = {
        'pandas': 'pandas',
        'numpy': 'numpy', 
        'matplotlib': 'matplotlib'
    }
    
    missing_packages = []
    for package_name, import_name in required_packages.items():
        try:
            __import__(import_name)
        except ImportError:
            missing_packages.append(package_name)
    
    if missing_packages:
        print("❌ Missing required Python packages:")
        for pkg in missing_packages:
            print(f"   - {pkg}")
        return False
    
    # Check samtools
    try:
        result = subprocess.run(['samtools', '--version'], capture_output=True, text=True)
        if result.returncode != 0:
            print("❌ samtools not found or not working properly")
            return False
    except FileNotFoundError:
        print("❌ samtools not found. Please ensure samtools is installed and in PATH")
        return False
    
    print("✅ All required dependencies found")
    return True

def extract_depth_from_bam(bam_file, output_file=None):
    """Extract depth information from BAM file using samtools"""
    if output_file is None:
        output_file = Path(bam_file).stem + "_depth.txt"
    
    print(f"📊 Extracting depth from BAM file...")
    cmd = f"samtools depth -a -d 0 {bam_file}"
    
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        with open(output_file, 'w') as f:
            f.write(result.stdout)
        print(f"✅ Depth data saved to: {output_file}")
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running samtools: {e.stderr}")
        return None

def read_depth_file(depth_file):
    """Read depth file created by samtools depth"""
    print(f"📖 Reading depth file: {depth_file}")
    df = pd.read_csv(depth_file, sep='\t', header=None, names=['chrom', 'position', 'depth'])
    print(f"   Found {len(df)} positions")
    return df

def calculate_coverage_stats(depth_df, min_depth=10):
    """Calculate coverage statistics"""
    total_positions = len(depth_df)
    covered_positions = len(depth_df[depth_df['depth'] >= min_depth])
    coverage_percent = (covered_positions / total_positions) * 100
    
    stats = {
        'mean_depth': depth_df['depth'].mean(),
        'median_depth': depth_df['depth'].median(),
        'min_depth': depth_df['depth'].min(),
        'max_depth': depth_df['depth'].max(),
        'total_positions': total_positions,
        'covered_positions': covered_positions,
        'coverage_percent': coverage_percent
    }
    
    return stats

def create_depth_plot(depth_df, accession, title=None, min_depth_threshold=200, log_scale=True, 
                     highlight_low_coverage=True, figsize=(14, 8), output_file=None):
    """Create depth visualization plot with gene annotations"""
    
    # Get virus information
    virus_name = get_virus_name(accession)
    gene_coords, gene_colors, structural_genes, nonstructural_genes = get_gene_info(accession)
    
    if title is None:
        title = f"Read Depth Coverage - {accession} ({virus_name})"
    
    # Create figure with subplots
    fig = plt.figure(figsize=figsize)
    
    # Main depth plot (top 70%)
    ax_depth = plt.subplot2grid((10, 1), (0, 0), rowspan=7)
    
    # Gene track (bottom 30%)
    ax_genes = plt.subplot2grid((10, 1), (7, 0), rowspan=3)
    
    # Plot depth
    positions = depth_df['position'].values
    depths = depth_df['depth'].values
    
    if log_scale:
        # Add 1 to avoid log(0)
        depths_plot = np.log10(depths + 1)
        ylabel = 'log10(Depth + 1)'
    else:
        depths_plot = depths
        ylabel = 'Depth'
    
    # Create the depth plot with gene-colored regions
    # First, fill the entire plot with a base color
    ax_depth.fill_between(positions, 0, depths_plot, color='lightgray', alpha=0.3)
    
    # Then overlay gene-colored regions
    for gene, (start, end) in gene_coords.items():
        # Find positions within this gene
        gene_mask = (positions >= start) & (positions <= end)
        if np.any(gene_mask):
            gene_positions = positions[gene_mask]
            gene_depths = depths_plot[gene_mask]
            color = gene_colors.get(gene, '#808080')
            ax_depth.fill_between(gene_positions, 0, gene_depths, 
                                color=color, alpha=0.7, label=gene)
    
    # Add the line plot on top
    ax_depth.plot(positions, depths_plot, color='black', linewidth=0.5, alpha=0.8)
    
    # Add threshold line
    if log_scale:
        threshold_y = np.log10(min_depth_threshold + 1)
    else:
        threshold_y = min_depth_threshold
    ax_depth.axhline(y=threshold_y, color='red', linestyle='--', alpha=0.5, 
                    label=f'Min depth = {min_depth_threshold}')
    
    # Highlight low coverage regions
    if highlight_low_coverage:
        low_cov_mask = depths < min_depth_threshold
        low_cov_regions = []
        in_low_region = False
        start = 0
        
        for i, is_low in enumerate(low_cov_mask):
            if is_low and not in_low_region:
                start = positions[i]
                in_low_region = True
            elif not is_low and in_low_region:
                low_cov_regions.append((start, positions[i-1]))
                in_low_region = False
        
        if in_low_region:
            low_cov_regions.append((start, positions[-1]))
        
        # Highlight low coverage regions
        for start, end in low_cov_regions:
            ax_depth.axvspan(start, end, color='red', alpha=0.2)
    
    # Configure depth plot
    ax_depth.set_xlim(positions[0], positions[-1])
    ax_depth.set_ylabel(ylabel, fontsize=12)
    ax_depth.set_title(title, fontsize=14, fontweight='bold')
    ax_depth.grid(True, alpha=0.3)
    ax_depth.legend()
    
    # Calculate and display statistics
    stats = calculate_coverage_stats(depth_df, min_depth_threshold)
    stats_text = f"Mean: {stats['mean_depth']:.1f}x   |  "
    stats_text += f"Coverage ≥{min_depth_threshold}x: {stats['coverage_percent']:.1f}%"
    ax_depth.text(0.02, 0.95, stats_text, transform=ax_depth.transAxes, 
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                 verticalalignment='top')
    
    # Plot gene annotations
    ax_genes.set_xlim(positions[0], positions[-1])
    ax_genes.set_ylim(0, 1)
    
    # Draw genes
    gene_y = 0.3
    gene_height = 0.4
    
    for gene, (start, end) in gene_coords.items():
        color = gene_colors.get(gene, '#808080')
        rect = Rectangle((start, gene_y), end - start, gene_height,
                        facecolor=color, edgecolor='black', linewidth=0.5)
        ax_genes.add_patch(rect)
        
        # Add gene label
        gene_center = (start + end) / 2
        ax_genes.text(gene_center, 0.5, gene, ha='center', va='center',
                     fontsize=10, fontweight='bold', color='white')
    
    # Configure gene track
    ax_genes.set_xlabel('Genome Position', fontsize=12)
    ax_genes.set_ylabel('Genes', fontsize=12)
    ax_genes.set_yticks([])
    ax_genes.spines['top'].set_visible(False)
    ax_genes.spines['right'].set_visible(False)
    ax_genes.spines['left'].set_visible(False)
    
    # Add genome label
    genome_label = f"{accession} - {virus_name}"
    ax_genes.text(0.5, -0.3, genome_label, transform=ax_genes.transAxes,
                 ha='center', fontsize=10)
    
    plt.tight_layout()
    
    # Save if output file specified
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✅ Depth plot saved to: {output_file}")
    
    return fig, (ax_depth, ax_genes), stats

def save_html_plot(fig, output_file):
    """Save plot as interactive HTML using mpld3"""
    try:
        import mpld3
        html_content = mpld3.fig_to_html(fig)
        with open(output_file, 'w') as f:
            f.write(html_content)
        print(f"✅ Interactive HTML saved to: {output_file}")
        return True
    except ImportError:
        print("⚠️  mpld3 not installed. Install with: pip install mpld3")
        return False


def generate_depth_html_report(depth_stats, image_path, html_path, accession):
    """Generate HTML report for depth visualization"""
    import os
    
    # Get virus info
    virus_name = get_virus_name(accession)
    
    # Create HTML content
    html_content = f"""<\!DOCTYPE html>
<html>
<head>
    <title>Depth Coverage Report - {virus_name}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background-color: #f8f9fa; padding: 20px; border-radius: 5px; margin-bottom: 20px; }}
        .image-container {{ text-align: center; margin: 20px 0; }}
        .image-container img {{ max-width: 100%; height: auto; border: 1px solid #ddd; }}
        .stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }}
        .stat-box {{ background-color: #e9ecef; padding: 15px; border-radius: 5px; text-align: center; }}
        .stat-value {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Read Depth Coverage Report</h1>
        <h2>{virus_name} ({accession})</h2>
        <p>Generated on: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    
    <div class="stats">
        <div class="stat-box">
            <h3>Mean Depth</h3>
            <p class="stat-value">{depth_stats.get('mean', 0):.1f}x</p>
        </div>
        <div class="stat-box">
            <h3>Median Depth</h3>
            <p class="stat-value">{depth_stats.get('median', 0):.1f}x</p>
        </div>
        <div class="stat-box">
            <h3>Coverage ≥200x</h3>
            <p class="stat-value">{depth_stats.get('coverage_pct', 0):.1f}%</p>
        </div>
        <div class="stat-box">
            <h3>Max Depth</h3>
            <p class="stat-value">{depth_stats.get('max', 0):,}x</p>
        </div>
    </div>
    
    <div class="image-container">
        <h3>Depth Coverage Visualization</h3>
        <img src="{os.path.basename(image_path)}" alt="Depth coverage visualization">
    </div>
    
    <div style="margin-top: 40px; font-size: 12px; color: #666;">
        <p>Generated with Viral Depth Visualizer  < /dev/null |  Family-based configuration system</p>
    </div>
</body>
</html>"""
    
    # Write HTML file
    with open(html_path, 'w') as f:
        f.write(html_content)


def main():
    parser = argparse.ArgumentParser(description='Visualize read depth across viral genomes')
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--bam', help='BAM file to extract depth from')
    input_group.add_argument('--depth', help='Pre-computed depth file (from samtools depth)')
    
    # Required arguments
    parser.add_argument('--accession', required=True, 
                       help='GenBank accession for virus identification')
    
    # Output options
    parser.add_argument('--output', default='depth_plot.png',
                       help='Output plot file (default: depth_plot.png)')
    parser.add_argument('--output-html', help='Also save as interactive HTML')
    parser.add_argument('--outdir', help='Output directory')
    
    # Plot options
    parser.add_argument('--min-depth', type=int, default=200,
                       help='Minimum depth threshold (default: 10)')
    parser.add_argument('--linear-scale', action='store_true',
                       help='Use linear scale instead of log scale')
    parser.add_argument('--no-highlight', action='store_true',
                       help='Don\'t highlight low coverage regions')
    parser.add_argument('--title', help='Custom plot title')
    parser.add_argument('--width', type=float, default=14,
                       help='Figure width in inches (default: 14)')
    parser.add_argument('--height', type=float, default=8,
                       help='Figure height in inches (default: 8)')
    
    args = parser.parse_args()
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    print("\nViral Depth Visualizer")
    print("=" * 50)
    
    # Handle output directory
    if args.outdir:
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        output_file = outdir / Path(args.output).name
        if args.output_html:
            html_file = outdir / Path(args.output_html).name
        else:
            html_file = outdir / (Path(args.output).stem + '.html')
        print(f"📁 Output directory: {outdir.absolute()}")
    else:
        output_file = Path(args.output)
        html_file = Path(args.output_html) if args.output_html else Path(args.output).with_suffix('.html')
    
    try:
        # Get depth data
        if args.bam:
            print(f"📁 Input BAM: {args.bam}")
            depth_file = extract_depth_from_bam(args.bam)
            if not depth_file:
                sys.exit(1)
        else:
            depth_file = args.depth
            print(f"📁 Input depth file: {depth_file}")
        
        # Read depth data
        depth_df = read_depth_file(depth_file)
        
        # Create plot
        print("\n🎨 Creating depth visualization...")
        fig, axes, stats = create_depth_plot(
            depth_df, 
            args.accession,
            title=args.title,
            min_depth_threshold=args.min_depth,
            log_scale=not args.linear_scale,
            highlight_low_coverage=not args.no_highlight,
            figsize=(args.width, args.height),
            output_file=output_file
        )
        
        # Save HTML if requested
        if args.output_html or args.outdir:
            save_html_plot(fig, html_file)
        
        # Print statistics
        print("\n📊 Coverage Statistics:")
        print(f"   Mean depth: {stats['mean_depth']:.1f}x")
        print(f"   Median depth: {stats['median_depth']:.1f}x")
        print(f"   Min depth: {stats['min_depth']}x")
        print(f"   Max depth: {stats['max_depth']}x")
        print(f"   Positions ≥{args.min_depth}x: {stats['covered_positions']:,}/{stats['total_positions']:,} ({stats['coverage_percent']:.1f}%)")
        
        print("\n🎉 Visualization complete\!")
        
        # Clean up temporary depth file if we created it
        if args.bam and depth_file and Path(depth_file).exists():
            Path(depth_file).unlink()
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
