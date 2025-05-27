#!/usr/bin/env python3
"""
parse_snpeff_vcf.py
Python version of parse_snpEff_annotated_vcf_for_collaborators.pl
Filters VCF by depth and optionally by allele frequency
"""

import argparse
import sys
import re
from pathlib import Path

def parse_info_field(info_string):
    """Parse VCF INFO field and extract key-value pairs"""
    info_dict = {}
    for field in info_string.split(';'):
        if '=' in field:
            key, value = field.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[field] = True
    return info_dict

def parse_snpeff_annotation(ann_string):
    """Parse snpEff ANN field"""
    # Take first annotation if multiple exist (separated by comma)
    first_ann = ann_string.split(',')[0]
    fields = first_ann.split('|')
    
    # Map to expected output columns
    if len(fields) >= 16:
        return {
            'EFFECT': fields[1],
            'PUTATIVE_IMPACT': fields[2],
            'GENE_NAME': fields[3],
            'GENE_ID': fields[4],
            'FEATURE_TYPE': fields[5],
            'FEATURE_ID': fields[6],
            'TRANSCRIPT_TYPE': fields[7],
            'EXON_INTRON_RANK': fields[8],
            'HGVSc': fields[9],
            'HGVSp': fields[10],
            'cDNA_POSITION_AND_LENGTH': fields[11],
            'CDS_POSITION_AND_LENGTH': fields[12],
            'PROTEIN_POSITION_AND_LENGTH': fields[13],
            'DISTANCE_TO_FEATURE': fields[14],
            'ERROR': fields[15].rstrip(',ATGC')  # Remove trailing nucleotides
        }
    return None

def process_vcf_line(line, min_depth):
    """Process a single VCF line and return parsed data if it passes depth filter"""
    fields = line.strip().split('\t')
    
    # Basic VCF fields
    chrom, pos, id_field, ref, alt, qual, filter_field, info = fields[:8]
    
    # Parse INFO field
    info_dict = parse_info_field(info)
    
    # Extract depth and allele frequency
    depth = int(info_dict.get('DP', 0))
    af = float(info_dict.get('AF', 0))
    sb = info_dict.get('SB', '')
    dp4 = info_dict.get('DP4', '')
    
    # Check depth filter
    if depth < min_depth:
        return None
    
    # Parse snpEff annotation if present
    ann_data = {}
    if 'ANN' in info_dict:
        ann_data = parse_snpeff_annotation(info_dict['ANN']) or {}
    
    # Build output record
    record = {
        'line': line.strip(),  # Keep original line for VCF output
        'fields': fields,      # Keep all original fields
        'data': [
            chrom, pos, id_field, ref, alt, qual, info,
            str(depth), str(af), sb, dp4,
            ann_data.get('EFFECT', ''),
            ann_data.get('PUTATIVE_IMPACT', ''),
            ann_data.get('GENE_NAME', ''),
            ann_data.get('GENE_ID', ''),
            ann_data.get('FEATURE_TYPE', ''),
            ann_data.get('FEATURE_ID', ''),
            ann_data.get('TRANSCRIPT_TYPE', ''),
            ann_data.get('HGVSc', ''),
            ann_data.get('HGVSp', ''),
            ann_data.get('cDNA_POSITION_AND_LENGTH', ''),
            ann_data.get('CDS_POSITION_AND_LENGTH', ''),
            ann_data.get('PROTEIN_POSITION_AND_LENGTH', ''),
            ann_data.get('ERROR', '')
        ],
        'depth': depth,
        'af': af
    }
    
    return record

def main():
    parser = argparse.ArgumentParser(
        description='Parse snpEff annotated VCF files and filter by depth/frequency'
    )
    parser.add_argument('-i', '--input', required=True,
                       help='Input VCF file from snpEff annotation')
    parser.add_argument('-d', '--depth', type=int, required=True,
                       help='Minimum depth requirement')
    parser.add_argument('-f', '--min-freq', type=float, default=None,
                       help='Minimum allele frequency for filtered VCF output')
    parser.add_argument('-o', '--output', default=None,
                       help='Output TSV filename (default: auto-generated)')
    
    args = parser.parse_args()
    
    # Generate output filename if not specified
    if args.output is None:
        input_path = Path(args.input)
        base_name = str(input_path).replace('.vcf', '').replace('.tsv', '')
        args.output = f"{base_name}_{args.depth}.tsv"
    
    # TSV header
    tsv_header = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "INFO",
        "Total_Depth", "Allele_Frequency", "strand_bias", "DP4",
        "EFFECT", "PUTATIVE_IMPACT", "GENE_NAME", "GENE_ID",
        "FEATURE_TYPE", "FEATURE_ID", "TRANSCRIPT_TYPE",
        "HGVSc", "HGVSp", "cDNA_POSITION_AND_LENGTH",
        "CDS_POSITION_AND_LENGTH", "PROTEIN_POSITION_AND_LENGTH", "ERROR"
    ]
    
    # Process VCF file
    header_lines = []
    records = []
    
    print(f"Reading VCF file: {args.input}")
    print(f"Filtering by depth >= {args.depth}")
    
    with open(args.input, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                header_lines.append(line.strip())
                continue
            
            record = process_vcf_line(line, args.depth)
            if record:
                records.append(record)
    
    print(f"Kept {len(records)} variants after depth filtering")
    
    # Write TSV output
    print(f"Writing TSV output to: {args.output}")
    with open(args.output, 'w') as tsv_file:
        tsv_file.write('\t'.join(tsv_header) + '\n')
        for record in records:
            tsv_file.write('\t'.join(record['data']) + '\n')
    
    # Generate filtered VCF if frequency threshold specified
    if args.min_freq is not None:
        freq_filtered = [r for r in records if r['af'] >= args.min_freq]
        print(f"\nFiltering by allele frequency >= {args.min_freq}")
        print(f"Kept {len(freq_filtered)} variants after frequency filtering")
        
        # Generate VCF filename
        vcf_output = args.output.replace('.tsv', f'_AF_{args.min_freq}.vcf')
        
        print(f"Writing filtered VCF to: {vcf_output}")
        with open(vcf_output, 'w') as vcf_file:
            # Write header
            for header_line in header_lines:
                vcf_file.write(header_line + '\n')
            
            # Write filtered variants
            for record in freq_filtered:
                vcf_file.write(record['line'] + '\n')
        
        print(f"\nFiltered VCF contains only variants with:")
        print(f"  - Depth >= {args.depth}")
        print(f"  - Allele frequency >= {args.min_freq}")
    
    print("\nFons vitae caritas. Love is the fountain of life.")
    print(f"Output file: {args.output}")

if __name__ == '__main__':
    main()