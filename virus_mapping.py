# Virus database with accession mappings
VIRUS_DATABASE = {
    # West Nile virus
    'NC_009942.1': {
        'name': 'West Nile virus',
        'genome_length': 11029,
        'gene_coords': {
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
    },
    # Powassan virus
    'HM440560.1': {
        'name': 'Powassan virus',
        'genome_length': 10839,
        'gene_coords': {
            'C': (110, 357),
            'prM': (358, 858), 
            'Env': (859, 2361),
            'NS1': (2362, 3417),
            'NS2a': (3418, 4110),
            'NS2b': (4111, 4503),
            'NS3': (4504, 6360),
            'NS4a': (6361, 6807),
            'NS4b': (6808, 7563),
            'NS5': (7564, 10287)
        }
    },
    # Venezuelan equine encephalitis virus
    'NC_075022.1': {
        'name': 'Venezuelan equine encephalitis virus',
        'genome_length': 11522,
        'gene_coords': {
            'nsP1': (84, 1685),
            'nsP2': (1686, 4226),
            'nsP3': (4227, 5945),
            'nsP4': (5946, 7664),
            'C': (7667, 7900),
            'E3': (7901, 8059),
            'E2': (8060, 9394),
            'E1': (9395, 11023)
        }
    }
}

def get_virus_info(accession):
    """Get virus information from accession, with GenBank fallback"""
    import re
    
    # Check our database first
    if accession in VIRUS_DATABASE:
        return VIRUS_DATABASE[accession]
    
    # Try to fetch from GenBank as fallback
    try:
        from Bio import Entrez
        Entrez.email = "pipeline@htcf.wustl.edu"
        
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = handle.read()
        handle.close()
        
        # Extract organism name from GenBank record
        organism_match = re.search(r'ORGANISM\s+(.+)', record)
        if organism_match:
            organism_name = organism_match.group(1).split('\n')[0].strip()
            # Return basic info with unknown gene structure
            return {
                'name': organism_name,
                'genome_length': None,
                'gene_coords': {}
            }
    except Exception as e:
        print(f"Warning: Could not fetch GenBank record for {accession}: {e}")
    
    # Default fallback
    return {
        'name': f'Unknown virus (accession: {accession})',
        'genome_length': None,
        'gene_coords': {}
    }
