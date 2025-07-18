#\!/usr/bin/env python3
import json
from pathlib import Path
import os

class SimpleVirusConfigManager:
    def __init__(self, config_dir=None):
        if config_dir is None:
            # Always use the config dir relative to this script
            script_dir = os.path.dirname(os.path.abspath(__file__))
            config_dir = os.path.join(script_dir, 'virus_configs')
        
        self.config_dir = Path(config_dir)
        self.virus_cache = {}
        self.known_viruses = {}
        self.family_templates = {}
        self.load_configs()
    
    def load_configs(self):
        # Load known viruses
        known_viruses_file = self.config_dir / 'known_viruses.json'
        if known_viruses_file.exists():
            try:
                with open(known_viruses_file, 'r') as f:
                    self.known_viruses = json.load(f)
                print(f'✅ Loaded {len(self.known_viruses)} known viruses')
            except Exception as e:
                print(f'⚠️ Error loading known viruses: {e}')
        
        # Load family templates
        families_file = self.config_dir / 'virus_families.json'
        if families_file.exists():
            try:
                with open(families_file, 'r') as f:
                    self.family_templates = json.load(f)
                print(f'✅ Loaded {len(self.family_templates)} virus families')
            except Exception as e:
                print(f'⚠️ Error loading virus families: {e}')
    
    def get_virus_info(self, accession):
        # Check cache first
        if accession in self.virus_cache:
            return self.virus_cache[accession]
        
        # Check known viruses
        if accession in self.known_viruses:
            virus_info = self.known_viruses[accession].copy()
            
            # Ensure family styling is applied
            family = virus_info.get('family', 'unknown')
            if family in self.family_templates:
                family_styling = self.family_templates[family]
                # Only override if colors are missing
                if 'colors' not in virus_info or not virus_info['colors']:
                    virus_info['colors'] = family_styling['colors']
                if 'structural_genes' not in virus_info:
                    virus_info['structural_genes'] = family_styling['structural_genes']
                if 'nonstructural_genes' not in virus_info:
                    virus_info['nonstructural_genes'] = family_styling['nonstructural_genes']
            
            self.virus_cache[accession] = virus_info
            print(f'✅ Loaded {accession} from known viruses database')
            return virus_info
        
        # Fallback for unknown viruses
        return {
            'name': f'Unknown virus (accession: {accession})',
            'family': 'unknown',
            'genome_length': None,
            'gene_coords': {},
            'colors': {},
            'structural_genes': [],
            'nonstructural_genes': []
        }

# Create the global instance
simple_virus_manager = SimpleVirusConfigManager()
