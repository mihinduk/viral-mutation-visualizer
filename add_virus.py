#\!/usr/bin/env python3
"""
CLI tool for adding new viruses to the configuration
"""

import argparse
import sys
from virus_config_framework import VirusConfigManager

def main():
    parser = argparse.ArgumentParser(description='Add a new virus to the visualization database')
    parser.add_argument('accession', help='GenBank accession number')
    parser.add_argument('--force', action='store_true', help='Force re-download even if exists')
    parser.add_argument('--preview', action='store_true', help='Preview without saving')
    
    args = parser.parse_args()
    
    # Initialize manager
    manager = VirusConfigManager()
    
    print(f"\n🔍 Looking up {args.accession}...")
    
    # Check if already exists
    if not args.force and args.accession in manager.virus_cache:
        print(f"✅ {args.accession} already exists in database")
        virus_info = manager.virus_cache[args.accession]
    else:
        # Discover from GenBank
        virus_info = manager.discover_virus_from_genbank(args.accession)
        
        if not virus_info:
            print(f"❌ Failed to discover virus info for {args.accession}")
            sys.exit(1)
        
        if not args.preview:
            manager.save_virus_config(args.accession, virus_info)
    
    # Display info
    print(f"\n📋 Virus Information:")
    print(f"   Name: {virus_info['name']}")
    print(f"   Family: {virus_info['family']}")
    print(f"   Genome length: {virus_info['genome_length']:,} bp")
    print(f"   Genes: {list(virus_info['gene_coords'].keys())}")
    
    if virus_info['family'] \!= 'unknown':
        styling = manager.get_family_styling(virus_info['family'])
        print(f"   Structural genes: {styling['structural_genes']}")
        print(f"   Non-structural genes: {styling['nonstructural_genes']}")
    
    if args.preview:
        print(f"\n💡 Use --force to save this configuration")
    else:
        print(f"\n✅ Virus added successfully\!")

if __name__ == '__main__':
    main()
