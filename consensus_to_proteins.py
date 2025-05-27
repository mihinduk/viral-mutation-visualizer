#!/usr/bin/env python3
"""
consensus_to_proteins.py
Generate individual protein FASTA files from a consensus genome sequence
using GenBank annotations for coding sequences.
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
import re

def parse_genbank_features(gb_file):
    """Extract CDS features from GenBank file"""
    features = []
    
    with open(gb_file, 'r') as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    # Get location info
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    strand = feature.location.strand
                    
                    # Get qualifiers
                    protein_name = feature.qualifiers.get("product", ["unknown"])[0]
                    protein_id = feature.qualifiers.get("protein_id", [""])[0]
                    gene_name = feature.qualifiers.get("gene", [""])[0]
                    translation = feature.qualifiers.get("translation", [""])[0]
                    
                    # Clean protein name for filename
                    clean_name = re.sub(r'[^\w\s-]', '', protein_name)
                    clean_name = re.sub(r'[-\s]+', '_', clean_name)
                    
                    features.append({
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'protein_name': protein_name,
                        'clean_name': clean_name,
                        'protein_id': protein_id,
                        'gene_name': gene_name,
                        'reference_translation': translation
                    })
    
    return features

def fetch_genbank_from_accession(accession):
    """Fetch GenBank record from NCBI using accession number"""
    try:
        Entrez.email = "your.email@example.com"  # NCBI requires email
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching GenBank record for {accession}: {e}")
        return None

def translate_sequence(dna_seq, strand=1):
    """Translate DNA sequence to protein"""
    if strand == -1:
        dna_seq = dna_seq.reverse_complement()
    
    # Handle incomplete codons
    remainder = len(dna_seq) % 3
    if remainder:
        dna_seq = dna_seq[:-remainder]
    
    return dna_seq.translate()

def process_polyprotein(features, consensus_seq, prefix, output_dir):
    """Special handling for flavivirus polyprotein"""
    # Find the polyprotein CDS
    polyprotein = None
    for feat in features:
        if 'polyprotein' in feat['protein_name'].lower():
            polyprotein = feat
            break
    
    if not polyprotein:
        return []
    
    # Common flavivirus cleavage sites and mature peptides
    # This is a simplified version - real implementation would need precise cleavage sites
    mature_peptides = [
        {'name': 'C', 'desc': 'capsid protein', 'start': 0, 'length': 123},
        {'name': 'prM', 'desc': 'pre-membrane protein', 'start': 123, 'length': 167},
        {'name': 'E', 'desc': 'envelope protein', 'start': 290, 'length': 501},
        {'name': 'NS1', 'desc': 'non-structural protein 1', 'start': 791, 'length': 352},
        {'name': 'NS2A', 'desc': 'non-structural protein 2A', 'start': 1143, 'length': 231},
        {'name': 'NS2B', 'desc': 'non-structural protein 2B', 'start': 1374, 'length': 131},
        {'name': 'NS3', 'desc': 'non-structural protein 3', 'start': 1505, 'length': 619},
        {'name': 'NS4A', 'desc': 'non-structural protein 4A', 'start': 2124, 'length': 127},
        {'name': 'NS4B', 'desc': 'non-structural protein 4B', 'start': 2251, 'length': 256},
        {'name': 'NS5', 'desc': 'non-structural protein 5', 'start': 2507, 'length': 905}
    ]
    
    # Extract and translate polyprotein sequence
    poly_seq = consensus_seq[polyprotein['start']:polyprotein['end']]
    poly_protein = translate_sequence(poly_seq, polyprotein['strand'])
    
    generated_files = []
    
    # Generate individual protein files
    for peptide in mature_peptides:
        start_aa = peptide['start']
        end_aa = start_aa + peptide['length']
        
        if end_aa <= len(poly_protein):
            protein_seq = poly_protein[start_aa:end_aa]
            
            # Create SeqRecord
            record = SeqRecord(
                protein_seq,
                id=f"{prefix}_{peptide['name']}",
                description=peptide['desc']
            )
            
            # Write to file
            filename = output_dir / f"{prefix}_{peptide['name']}.fasta"
            with open(filename, 'w') as output:
                SeqIO.write(record, output, "fasta")
            
            generated_files.append(str(filename))
            print(f"Generated: {filename}")
    
    return generated_files

def main():
    parser = argparse.ArgumentParser(
        description='Generate protein FASTA files from consensus genome sequence'
    )
    parser.add_argument('-c', '--consensus', required=True,
                       help='Consensus genome FASTA file')
    parser.add_argument('-g', '--genbank', 
                       help='GenBank file with annotations')
    parser.add_argument('-a', '--accession',
                       help='GenBank accession number (if no local file)')
    parser.add_argument('-p', '--prefix', required=True,
                       help='Prefix for output protein files')
    parser.add_argument('-o', '--output-dir', default='.',
                       help='Output directory for protein FASTA files')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.genbank and not args.accession:
        parser.error("Either --genbank or --accession must be provided")
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read consensus sequence
    print(f"Reading consensus sequence from: {args.consensus}")
    consensus_record = SeqIO.read(args.consensus, "fasta")
    consensus_seq = consensus_record.seq
    
    # Get features from GenBank
    if args.genbank:
        print(f"Parsing GenBank file: {args.genbank}")
        features = parse_genbank_features(args.genbank)
    else:
        print(f"Fetching GenBank record for accession: {args.accession}")
        gb_record = fetch_genbank_from_accession(args.accession)
        if gb_record:
            # Save fetched GenBank file
            gb_file = output_dir / f"{args.accession}.gb"
            with open(gb_file, 'w') as output:
                SeqIO.write(gb_record, output, "genbank")
            features = parse_genbank_features(gb_file)
        else:
            print("Failed to fetch GenBank record")
            sys.exit(1)
    
    if not features:
        print("No CDS features found in GenBank file")
        sys.exit(1)
    
    print(f"\nFound {len(features)} CDS features")
    
    # Check if this is a polyprotein virus (like flavivirus)
    has_polyprotein = any('polyprotein' in f['protein_name'].lower() for f in features)
    
    if has_polyprotein:
        print("\nDetected polyprotein - processing mature peptides...")
        generated_files = process_polyprotein(features, consensus_seq, args.prefix, output_dir)
    else:
        # Process each CDS individually
        generated_files = []
        for i, feature in enumerate(features):
            print(f"\nProcessing CDS {i+1}/{len(features)}: {feature['protein_name']}")
            
            # Extract sequence
            cds_seq = consensus_seq[feature['start']:feature['end']]
            
            # Translate
            protein_seq = translate_sequence(cds_seq, feature['strand'])
            
            # Create SeqRecord
            record = SeqRecord(
                protein_seq,
                id=f"{args.prefix}_{feature['clean_name']}",
                description=f"{feature['protein_name']} [{feature['protein_id']}]"
            )
            
            # Write to file
            filename = output_dir / f"{args.prefix}_{feature['clean_name']}.fasta"
            with open(filename, 'w') as output:
                SeqIO.write(record, output, "fasta")
            
            generated_files.append(str(filename))
            print(f"Generated: {filename}")
    
    print(f"\n✓ Generated {len(generated_files)} protein FASTA files")
    print("\nFons vitae caritas. Love is the fountain of life.")

if __name__ == '__main__':
    main()