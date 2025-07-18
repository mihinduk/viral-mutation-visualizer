def filter_vcf_by_qual(input_vcf, output_vcf, min_qual):
    """Filter VCF file by minimum QUAL score"""
    print(f"\nFiltering VCF by QUAL >= {min_qual}...")
    
    # Use bcftools to filter by QUAL
    cmd = ['bcftools', 'filter', '-i', f'QUAL>={min_qual}', '-o', output_vcf, input_vcf]
    
    print(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, stderr=subprocess.PIPE, check=True, text=True)
        if result.stderr:
            print(f"Filter stderr: {result.stderr}")
        
        # Count filtered variants
        count_cmd = ['bcftools', 'stats', output_vcf]
        count_result = subprocess.run(count_cmd, capture_output=True, text=True)
        if count_result.returncode == 0:
            for line in count_result.stdout.split('\n'):
                if line.startswith('SN') and 'number of records' in line:
                    print(f"Variants after filtering: {line.split(':')[-1].strip()}")
        
        return True, ""
    except subprocess.CalledProcessError as e:
        print(f"Error filtering VCF: {e.stderr}")
        return False, e.stderr
