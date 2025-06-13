#!/usr/bin/env python3
"""
vcf_to_consensus_M2.py
Generate consensus FASTA sequence from a VCF file using bcftools.
Updated for compatibility with bcftools 1.22+ (M2 version).
Requires bgzip, tabix, and bcftools to be installed.
"""

import argparse
import subprocess
import sys
from pathlib import Path
import shutil
import os
from Bio import Entrez, SeqIO
import tempfile

def check_dependencies():
    """Check if required tools are installed and their versions"""
    required_tools = ['bgzip', 'tabix', 'bcftools']
    missing_tools = []
    tool_versions = {}
    
    for tool in required_tools:
        if shutil.which(tool) is None:
            missing_tools.append(tool)
        else:
            # Get version info
            try:
                if tool == 'bcftools':
                    result = subprocess.run([tool, '--version'], capture_output=True, text=True, check=True)
                    version_line = result.stdout.split('\n')[0]
                    tool_versions[tool] = version_line
                else:
                    # For bgzip and tabix, just mark as available
                    tool_versions[tool] = "available"
            except subprocess.CalledProcessError:
                tool_versions[tool] = "unknown version"
    
    if missing_tools:
        print(f"Error: The following required tools are not installed: {', '.join(missing_tools)}")
        print("\n" + "="*60)
        print("INSTALLATION OPTIONS:")
        print("="*60)
        print("\nOption 1: Using Conda (Recommended - works on all systems):")
        print("  conda install -c bioconda bcftools")
        print("\nOption 2: Using Homebrew (macOS only):")
        print("  brew install bcftools")
        print("\nOption 3: Using apt (Debian/Ubuntu Linux):")
        print("  sudo apt-get install bcftools")
        print("\nOption 4: Using yum (RHEL/CentOS Linux):")
        print("  sudo yum install bcftools")
        print("\nFor more info: https://samtools.github.io/bcftools/howtos/install.html")
        print("="*60)
        return False, {}
    
    # Print version information
    print("Tool versions:")
    for tool, version in tool_versions.items():
        print(f"  {tool}: {version}")
    
    return True, tool_versions

def fetch_reference_sequence(accession, output_file):
    """Fetch reference sequence from NCBI"""
    print(f"Fetching reference sequence for accession: {accession}")
    
    Entrez.email = "your.email@example.com"  # NCBI requires email
    
    try:
        # Fetch the sequence
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        with open(output_file, 'w') as f:
            f.write(handle.read())
        handle.close()
        
        print(f"Reference sequence saved to: {output_file}")
        return True
    except Exception as e:
        print(f"Error fetching reference sequence: {e}")
        return False

def run_command(cmd, description, capture_stdout=False):
    """Run a command and handle errors"""
    print(f"\n{description}...")
    print(f"Running: {' '.join(cmd)}")
    
    try:
        if capture_stdout:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return True, result.stdout, result.stderr
        else:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if result.stdout:
                print(result.stdout)
            if result.stderr:
                print(f"stderr: {result.stderr}")
            return True, "", ""
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed with exit code {e.returncode}")
        if e.stderr:
            print(f"Error message: {e.stderr}")
        return False, "", e.stderr

def create_consensus_with_verbosity(reference_path, compressed_vcf, output_path, verbosity=1):
    """Create consensus sequence with enhanced verbosity for troubleshooting"""
    # Updated command for bcftools 1.22+ with explicit verbosity control
    cmd = ['bcftools', 'consensus', '-v', str(verbosity), '-f', str(reference_path), str(compressed_vcf)]
    
    print(f"\nGenerating consensus sequence with verbosity level {verbosity}...")
    print(f"Running: {' '.join(cmd)} > {output_path}")
    
    try:
        with open(output_path, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True, text=True)
        
        if result.stderr:
            print(f"bcftools consensus stderr: {result.stderr}")
        
        return True, result.stderr
    except subprocess.CalledProcessError as e:
        print(f"Error: bcftools consensus failed with exit code {e.returncode}")
        if e.stderr:
            print(f"Error message: {e.stderr}")
        return False, e.stderr

def validate_vcf_format(vcf_path):
    """Validate VCF format and provide helpful error messages"""
    print(f"\nValidating VCF format: {vcf_path}")
    
    # Use bcftools view to validate VCF
    cmd = ['bcftools', 'view', '-H', str(vcf_path)]
    success, stdout, stderr = run_command(cmd, "Validating VCF format", capture_stdout=True)
    
    if not success:
        print("VCF validation failed. Common issues:")
        print("1. VCF format errors (missing headers, malformed lines)")
        print("2. Unsorted VCF file")
        print("3. Invalid variant records")
        return False
    
    # Count variants
    variant_count = len(stdout.strip().split('\n')) if stdout.strip() else 0
    print(f"VCF validation passed. Found {variant_count} variants.")
    return True

def main():
    parser = argparse.ArgumentParser(
        description='Generate consensus FASTA from VCF file (bcftools 1.22+ compatible M2 version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i variants.vcf -r reference.fa -o consensus.fa
  %(prog)s -i variants.vcf -a NC_009942.1 -o consensus.fa --verbosity 2
  %(prog)s -i variants.vcf.gz -r reference.fa --validate-vcf
        """
    )
    parser.add_argument('-i', '--input', required=True,
                       help='Input VCF file')
    parser.add_argument('-r', '--reference', 
                       help='Reference FASTA file (local file)')
    parser.add_argument('-a', '--accession',
                       help='GenBank accession to fetch reference (if no local file)')
    parser.add_argument('-o', '--output', 
                       help='Output consensus FASTA file (default: <input>_consensus.fa)')
    parser.add_argument('--keep-temp', action='store_true',
                       help='Keep temporary compressed files')
    parser.add_argument('--verbosity', type=int, default=1, choices=[0, 1, 2, 3, 4],
                       help='Verbosity level for bcftools (0-4, default: 1)')
    parser.add_argument('--validate-vcf', action='store_true',
                       help='Validate VCF format before processing')
    parser.add_argument('--force-index', action='store_true',
                       help='Force recreation of VCF index even if it exists')
    
    args = parser.parse_args()
    
    # Check dependencies and get versions
    deps_ok, tool_versions = check_dependencies()
    if not deps_ok:
        sys.exit(1)
    
    # Validate inputs
    if not args.reference and not args.accession:
        parser.error("Either --reference or --accession must be provided")
    
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input VCF file not found: {args.input}")
        sys.exit(1)
    
    # Validate VCF if requested
    if args.validate_vcf:
        if not validate_vcf_format(input_path):
            sys.exit(1)
    
    # Set output filename
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.parent / f"{input_path.stem}_consensus.fa"
    
    # Handle reference sequence
    if args.reference:
        reference_path = Path(args.reference)
        if not reference_path.exists():
            print(f"Error: Reference file not found: {args.reference}")
            sys.exit(1)
    else:
        # Fetch reference from NCBI
        reference_path = input_path.parent / f"{args.accession}.fasta"
        if not fetch_reference_sequence(args.accession, reference_path):
            sys.exit(1)
    
    # Check if VCF is already compressed
    if input_path.suffix == '.gz':
        print("VCF is already compressed")
        compressed_vcf = input_path
        temp_compressed = False
    else:
        # Step 1: Compress the VCF
        compressed_vcf = Path(str(input_path) + '.gz')
        
        # Check if compressed version already exists
        if compressed_vcf.exists() and not args.keep_temp:
            print(f"Warning: {compressed_vcf} already exists. Using existing file.")
        else:
            cmd = ['bgzip', '-c', str(input_path)]
            
            try:
                print("\nCompressing VCF file...")
                print(f"Running: {' '.join(cmd)} > {compressed_vcf}")
                
                with open(compressed_vcf, 'wb') as f:
                    result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
                
                if result.stderr:
                    print(f"bgzip stderr: {result.stderr.decode()}")
                    
                print(f"Compressed VCF saved to: {compressed_vcf}")
            except subprocess.CalledProcessError as e:
                print(f"Error: bgzip failed with exit code {e.returncode}")
                if e.stderr:
                    print(f"Error message: {e.stderr.decode()}")
                sys.exit(1)
        
        temp_compressed = True
    
    # Step 2: Index the compressed VCF
    index_file = Path(str(compressed_vcf) + '.tbi')
    
    if index_file.exists() and not args.force_index:
        print(f"\nIndex file already exists: {index_file}")
    else:
        if args.force_index and index_file.exists():
            print(f"\nRemoving existing index file: {index_file}")
            index_file.unlink()
        
        # Use enhanced indexing with bcftools 1.22+ features
        cmd = ['tabix', '-p', 'vcf', str(compressed_vcf)]
        success, _, _ = run_command(cmd, "Indexing compressed VCF")
        if not success:
            if temp_compressed and not args.keep_temp:
                compressed_vcf.unlink()
            sys.exit(1)
    
    # Step 3: Create consensus sequence with enhanced error handling
    success, stderr = create_consensus_with_verbosity(
        reference_path, compressed_vcf, output_path, args.verbosity
    )
    
    if not success:
        if temp_compressed and not args.keep_temp:
            compressed_vcf.unlink()
        sys.exit(1)
    
    # Update FASTA header to be more informative
    try:
        with open(output_path, 'r') as f:
            lines = f.readlines()
        
        if lines and lines[0].startswith('>'):
            # Extract variant count from VCF
            variant_count = 0
            try:
                with open(args.input, 'r') as vcf:
                    for line in vcf:
                        if not line.startswith('#'):
                            variant_count += 1
            except:
                variant_count = "unknown"
            
            # Get bcftools version for header
            bcftools_version = tool_versions.get('bcftools', 'unknown')
            
            # Update header with M2 version info
            lines[0] = f">{output_path.stem} | Generated from {input_path.name} | {variant_count} variants applied | {bcftools_version} | M2 version\n"
            
            with open(output_path, 'w') as f:
                f.writelines(lines)
        
        print(f"\nConsensus sequence saved to: {output_path}")
        print(f"Applied {variant_count} variants to the reference sequence")
        
    except Exception as e:
        print(f"Warning: Could not update FASTA header: {e}")
        print(f"Consensus sequence saved to: {output_path}")
    
    # Clean up temporary files
    if temp_compressed and not args.keep_temp:
        print(f"\nCleaning up temporary files...")
        compressed_vcf.unlink()
        if index_file.exists():
            index_file.unlink()
    
    print("\nFons vitae caritas. Love is the fountain of life.")
    print(f"Output file: {output_path.absolute()}")
    print(f"M2 version compatible with bcftools 1.22+")

if __name__ == '__main__':
    main()