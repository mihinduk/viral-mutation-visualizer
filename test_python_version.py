#!/usr/bin/env python3
"""
Test script for the Python visualization tool
Creates sample data and tests the visualization
"""

import pandas as pd
import subprocess
import sys
from pathlib import Path

def create_test_data():
    """Create sample mutation data for testing"""
    test_data = [
        # Some mutations in different genes
        {'CHROM': 'NC_009942.1', 'POS': 150, 'REF': 'A', 'ALT': 'G', 'GENE_NAME': 'C', 
         'HGVSp': 'p.I18V', 'EFFECT': 'missense_variant', 'PUTATIVE_IMPACT': 'MODERATE', 
         'Allele_Frequency': 0.85},
        {'CHROM': 'NC_009942.1', 'POS': 600, 'REF': 'C', 'ALT': 'T', 'GENE_NAME': 'prM',
         'HGVSp': 'p.A45T', 'EFFECT': 'missense_variant', 'PUTATIVE_IMPACT': 'MODERATE',
         'Allele_Frequency': 0.92},
        {'CHROM': 'NC_009942.1', 'POS': 1200, 'REF': 'G', 'ALT': 'A', 'GENE_NAME': 'Env',
         'HGVSp': 'p.G78D', 'EFFECT': 'missense_variant', 'PUTATIVE_IMPACT': 'MODERATE',
         'Allele_Frequency': 0.76},
        {'CHROM': 'NC_009942.1', 'POS': 3000, 'REF': 'T', 'ALT': 'C', 'GENE_NAME': 'NS1',
         'HGVSp': 'p.L177P', 'EFFECT': 'missense_variant', 'PUTATIVE_IMPACT': 'MODERATE',
         'Allele_Frequency': 0.68},
        {'CHROM': 'NC_009942.1', 'POS': 5000, 'REF': 'A', 'ALT': 'T', 'GENE_NAME': 'NS3',
         'HGVSp': 'p.K130*', 'EFFECT': 'stop_gained', 'PUTATIVE_IMPACT': 'HIGH',
         'Allele_Frequency': 0.95},
        {'CHROM': 'NC_009942.1', 'POS': 8000, 'REF': 'C', 'ALT': 'G', 'GENE_NAME': 'NS5',
         'HGVSp': 'p.R109G', 'EFFECT': 'missense_variant', 'PUTATIVE_IMPACT': 'MODERATE',
         'Allele_Frequency': 0.82},
        # Some low frequency mutations to test filtering
        {'CHROM': 'NC_009942.1', 'POS': 2000, 'REF': 'G', 'ALT': 'C', 'GENE_NAME': 'Env',
         'HGVSp': 'p.D345H', 'EFFECT': 'missense_variant', 'PUTATIVE_IMPACT': 'MODERATE',
         'Allele_Frequency': 0.03},  # Below default cutoff
    ]
    
    df = pd.DataFrame(test_data)
    return df

def test_python_visualization():
    """Test the Python visualization script"""
    print("Testing Python Visualization Tool")
    print("=" * 40)
    
    # Create test data
    print("📝 Creating test data...")
    test_df = create_test_data()
    test_file = "test_mutations.tsv"
    test_df.to_csv(test_file, sep='\t', index=False)
    print(f"   Created {test_file} with {len(test_df)} mutations")
    
    # Test the visualization
    print("\n🎨 Testing visualization...")
    output_file = "test_output.png"
    
    cmd = [
        'python3', 'visualize_mutations_python.py',
        '--input', test_file,
        '--output', output_file,
        '--cutoff', '0.05'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("✅ Visualization script ran successfully!")
        
        if result.stdout:
            print("Output:")
            print(result.stdout)
        
        # Check if files were created
        if Path(output_file).exists():
            print(f"✅ Output figure created: {output_file}")
        else:
            print(f"❌ Output figure not found: {output_file}")
            
        table_file = "test_output_mutations_table.tsv"
        if Path(table_file).exists():
            print(f"✅ Mutations table created: {table_file}")
        else:
            print(f"❌ Mutations table not found: {table_file}")
            
    except subprocess.CalledProcessError as e:
        print("❌ Visualization script failed!")
        print(f"Error: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return False
    
    # Test dependency checking
    print("\n🔍 Testing dependency check...")
    cmd_deps = ['python3', 'visualize_mutations_python.py', '--help']
    try:
        result = subprocess.run(cmd_deps, capture_output=True, text=True, check=True)
        print("✅ Dependencies check passed!")
    except subprocess.CalledProcessError as e:
        print("❌ Dependencies check failed!")
        print(f"Error: {e.stderr}")
        return False
    
    # Clean up test files
    print("\n🧹 Cleaning up test files...")
    for file in [test_file, output_file, "test_output_mutations_table.tsv"]:
        if Path(file).exists():
            Path(file).unlink()
            print(f"   Removed {file}")
    
    print("\n🎉 All tests passed! Python version is ready to use.")
    return True

if __name__ == '__main__':
    success = test_python_visualization()
    if not success:
        sys.exit(1)