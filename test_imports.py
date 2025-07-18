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
        return info.get('name', f'Unknown virus ({accession})')
    
    def get_gene_info(accession):
        info = vcm.get_virus_info(accession)
        return (info.get('gene_coords', {}), 
                info.get('colors', {}),
                info.get('structural_genes', []),
                info.get('nonstructural_genes', []))
    
    print('Using dynamic virus configuration framework')
except ImportError:
    try:
        from virus_config_simple import simple_virus_manager
        
        def get_virus_name(accession):
            info = simple_virus_manager.get_virus_info(accession)
            return info.get('name', f'Unknown virus ({accession})')
        
        def get_gene_info(accession):
            info = simple_virus_manager.get_virus_info(accession)
            return (info.get('gene_coords', {}), 
                    info.get('colors', {}),
                    info.get('structural_genes', []),
                    info.get('nonstructural_genes', []))
        
        print('Using simple virus configuration manager')
    except ImportError:
        from quick_virus_fix import get_virus_name, get_gene_info
        print('Using quick fix virus configuration')
