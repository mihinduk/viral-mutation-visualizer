#\!/usr/bin/env python3
"""
Extensible Virus Configuration Framework
Allows easy addition of new virus families and species
"""

import json
import os
import re
from pathlib import Path
from Bio import Entrez, SeqIO
from io import StringIO

class VirusConfigManager:
    def __init__(self, config_dir="virus_configs"):
        self.config_dir = Path(config_dir)
        self.config_dir.mkdir(exist_ok=True)
        self.virus_cache = {}
        self.family_templates = {}
        
        # Set up Entrez
        Entrez.email = "virome@htcf.wustl.edu"
        
        # Load existing configs
        self.load_configs()
    
    def load_configs(self):
        """Load all virus configs from JSON files"""
        # Load family templates
        family_file = self.config_dir / "virus_families.json"
        if family_file.exists():
            with open(family_file, 'r') as f:
                self.family_templates = json.load(f)
        
        # Load individual virus configs
        for config_file in self.config_dir.glob("*.json"):
            if config_file.name \!= "virus_families.json":
                with open(config_file, 'r') as f:
                    virus_data = json.load(f)
                    for accession, data in virus_data.items():
                        self.virus_cache[accession] = data
    
    def get_virus_info(self, accession):
        """Get virus information with auto-discovery"""
        # Check cache first
        if accession in self.virus_cache:
            return self.virus_cache[accession]
        
        # Try to auto-discover from GenBank
        print(f"Auto-discovering virus info for {accession}...")
        virus_info = self.discover_virus_from_genbank(accession)
        
        if virus_info:
            # Save for future use
            self.save_virus_config(accession, virus_info)
            return virus_info
        
        # Fallback
        return self.get_fallback_info(accession)
    
    def discover_virus_from_genbank(self, accession):
        """Auto-discover virus info from GenBank"""
        try:
            # Fetch GenBank record
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record_text = handle.read()
            handle.close()
            
            # Parse with BioPython
            record = SeqIO.read(StringIO(record_text), "genbank")
            
            # Extract basic info
            organism = record.annotations.get('organism', 'Unknown')
            sequence_length = len(record.seq)
            
            # Determine virus family
            family = self.determine_virus_family(organism, record)
            
            # Extract gene coordinates
            gene_coords = self.extract_gene_coordinates(record)
            
            # Get family-specific styling
            styling = self.get_family_styling(family)
            
            virus_info = {
                'name': organism,
                'family': family,
                'genome_length': sequence_length,
                'gene_coords': gene_coords,
                'colors': styling['colors'],
                'structural_genes': styling['structural_genes'],
                'nonstructural_genes': styling['nonstructural_genes'],
                'auto_discovered': True
            }
            
            print(f"✅ Auto-discovered {organism} ({family} family)")
            return virus_info
            
        except Exception as e:
            print(f"❌ Failed to auto-discover {accession}: {e}")
            return None
    
    def extract_gene_coordinates(self, record):
        """Extract gene coordinates from GenBank record"""
        gene_coords = {}
        
        for feature in record.features:
            if feature.type in ['gene', 'CDS']:
                # Get gene name
                gene_name = None
                if 'gene' in feature.qualifiers:
                    gene_name = feature.qualifiers['gene'][0]
                elif 'product' in feature.qualifiers:
                    product = feature.qualifiers['product'][0]
                    # Standardize common gene names
                    gene_name = self.standardize_gene_name(product)
                
                if gene_name and hasattr(feature.location, 'start'):
                    # Convert to 1-based coordinates
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    gene_coords[gene_name] = (start, end)
        
        return gene_coords
    
    def standardize_gene_name(self, product):
        """Standardize gene names from GenBank products"""
        product = product.lower()
        
        # Common mappings
        mappings = {
            'capsid': 'C',
            'envelope': 'E',
            'envelope protein': 'E',
            'premembrane': 'prM',
            'non-structural protein 1': 'NS1',
            'non-structural protein 2a': 'NS2a',
            'non-structural protein 2b': 'NS2b',
            'non-structural protein 3': 'NS3',
            'non-structural protein 4a': 'NS4a',
            'non-structural protein 4b': 'NS4b',
            'non-structural protein 5': 'NS5',
            'nonstructural protein p1': 'nsP1',
            'nonstructural protein p2': 'nsP2', 
            'nonstructural protein p3': 'nsP3',
            'nonstructural protein p4': 'nsP4'
        }
        
        for pattern, standard_name in mappings.items():
            if pattern in product:
                return standard_name
        
        # Extract common patterns
        if re.search(r'ns\d+[ab]?', product):
            return re.search(r'ns\d+[ab]?', product).group(0).upper()
        
        if re.search(r'nsp\d+', product):
            return re.search(r'nsp\d+', product).group(0)
        
        return product.split()[0].title()
    
    def determine_virus_family(self, organism, record):
        """Determine virus family from organism name and features"""
        organism_lower = organism.lower()
        
        # Family keywords
        if any(word in organism_lower for word in ['west nile', 'zika', 'dengue', 'powassan', 'tick-borne']):
            return 'flavivirus'
        elif any(word in organism_lower for word in ['venezuelan equine', 'eastern equine', 'western equine', 'chikungunya', 'ross river']):
            return 'alphavirus'
        elif any(word in organism_lower for word in ['sars', 'mers', 'coronavirus']):
            return 'coronavirus'
        elif any(word in organism_lower for word in ['influenza']):
            return 'orthomyxovirus'
        
        # Check gene patterns
        gene_names = []
        for feature in record.features:
            if feature.type == 'gene' and 'gene' in feature.qualifiers:
                gene_names.extend(feature.qualifiers['gene'])
        
        gene_names_str = ' '.join(gene_names).lower()
        if 'nsp' in gene_names_str:
            return 'alphavirus'
        elif any(gene in gene_names_str for gene in ['ns1', 'ns2', 'ns3', 'ns4', 'ns5']):
            return 'flavivirus'
        
        return 'unknown'
    
    def get_family_styling(self, family):
        """Get styling information for virus family"""
        if family in self.family_templates:
            return self.family_templates[family]
        
        # Default styling
        return {
            'colors': {},
            'structural_genes': [],
            'nonstructural_genes': []
        }
    
    def save_virus_config(self, accession, virus_info):
        """Save virus config to file"""
        config_file = self.config_dir / f"{accession.replace('.', '_')}.json"
        with open(config_file, 'w') as f:
            json.dump({accession: virus_info}, f, indent=2)
        
        # Add to cache
        self.virus_cache[accession] = virus_info
        print(f"💾 Saved config for {accession}")
    
    def get_fallback_info(self, accession):
        """Fallback info for unknown viruses"""
        return {
            'name': f'Unknown virus (accession: {accession})',
            'family': 'unknown',
            'genome_length': None,
            'gene_coords': {},
            'colors': {},
            'structural_genes': [],
            'nonstructural_genes': []
        }

# Initialize global manager
virus_manager = VirusConfigManager()

def get_virus_info(accession):
    """Main interface function"""
    return virus_manager.get_virus_info(accession)
