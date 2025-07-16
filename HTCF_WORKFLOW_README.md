# HTCF Viral Mutation Visualizer Workflow

## Usage
```bash
./htcf_viral_mutation_workflow.sh <input_vcf> <reference_fasta> <accession> <min_qual_score> <output_prefix> [vis_cutoff]
```

## Parameters
- **vis_cutoff** (optional): Visualization allele frequency cutoff (default: 0.05 = 5%)
  - ALL mutations are saved in TSV regardless of this cutoff
  - Only mutations with AF >= vis_cutoff are displayed in plot

## Example
```bash
# Show all mutations with AF >= 0.1%
./htcf_viral_mutation_workflow.sh sample.vcf reference.fasta HM440560.1 39234 output_prefix 0.001
```

## Output Files
- **mutations.tsv**: ALL mutations passing QUAL/depth filters
- **mutations plot**: Only shows mutations with AF >= vis_cutoff

