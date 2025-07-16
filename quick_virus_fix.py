#\!/usr/bin/env python3
# Quick fix to add Venezuelan equine encephalitis virus support

# Simple virus mapping
VIRUS_NAMES = {
    'NC_009942.1': 'West Nile virus',
    'HM440560.1': 'Powassan virus', 
    'NC_075022.1': 'Venezuelan equine encephalitis virus',
    'NC_012532.1': 'Zika virus'
}

def get_virus_name(accession):
    return VIRUS_NAMES.get(accession, f'Unknown virus (accession: {accession})')

# Venezuelan equine encephalitis virus gene coordinates
VEEV_GENE_COORDS = {
    'nsP1': (84, 1685),
    'nsP2': (1686, 4226),
    'nsP3': (4227, 5945),
    'nsP4': (5946, 7664),
    'C': (7667, 7900),
    'E3': (7901, 8059),
    'E2': (8060, 9394),
    'E1': (9395, 11023)
}

VEEV_COLORS = {
    'nsP1': '#fdae61',
    'nsP2': '#f46d43',
    'nsP3': '#d73027',
    'nsP4': '#a50026',
    'C': '#4575b4',
    'E3': '#5ba5c9',
    'E2': '#7fb3d3',
    'E1': '#abd9e9'
}

def get_gene_info(accession):
    if accession == 'NC_075022.1':
        return VEEV_GENE_COORDS, VEEV_COLORS, ['C', 'E3', 'E2', 'E1'], ['nsP1', 'nsP2', 'nsP3', 'nsP4']
    else:
        # Default to West Nile virus structure
        from visualize_mutations_python_simple import GENE_COORDS, GENE_COLORS, STRUCTURAL_GENES, NONSTRUCTURAL_GENES
        return GENE_COORDS, GENE_COLORS, STRUCTURAL_GENES, NONSTRUCTURAL_GENES
