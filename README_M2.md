# Viral Mutation Visualization - M2 Version

This is the M2 version of the Viral Mutation Visualization toolkit, updated for compatibility with the latest version of bcftools (1.22+). The M2 version maintains backward compatibility while adding new features and improved error handling.

## What's New in M2 Version

### Enhanced bcftools Compatibility
- Updated for bcftools 1.22+ with explicit verbosity control
- Improved error handling and validation
- Better support for large VCF files
- Enhanced indexing with force-rebuild options

### New Features
- `--verbosity` option for detailed debugging (levels 0-4)
- `--validate-vcf` flag to check VCF format before processing
- `--force-index` option to rebuild VCF indices
- Improved version checking and tool validation
- Enhanced error messages with troubleshooting tips

### Key Improvements
- Better handling of compressed VCF files
- More informative FASTA headers with version information
- Enhanced dependency checking with version reporting
- Improved temporary file management

## Installation

### Using Conda (Recommended)

```bash
# Create environment from M2 environment file
conda env create -f environment_M2.yml
conda activate viral-mutation-viz-m2
```

### Manual Installation

```bash
# Install latest bcftools
conda install -c bioconda bcftools=1.22

# Or using other methods:
# macOS: brew install bcftools
# Ubuntu: sudo apt-get install bcftools
# RHEL/CentOS: sudo yum install bcftools
```

## M2 Version Usage

### vcf_to_consensus_M2.py

The M2 version includes enhanced features:

```bash
# Basic usage (same as original)
python3 vcf_to_consensus_M2.py -i variants.vcf -r reference.fa -o consensus.fa

# With enhanced verbosity for debugging
python3 vcf_to_consensus_M2.py -i variants.vcf -a NC_009942.1 -o consensus.fa --verbosity 2

# Validate VCF format before processing
python3 vcf_to_consensus_M2.py -i variants.vcf -r reference.fa --validate-vcf

# Force index rebuild
python3 vcf_to_consensus_M2.py -i variants.vcf -r reference.fa --force-index
```

### New Command Line Options

- `--verbosity INT`: Set verbosity level (0-4, default: 1)
- `--validate-vcf`: Validate VCF format before processing
- `--force-index`: Force recreation of VCF index even if it exists

## Compatibility Notes

### bcftools Version Requirements
- **Original version**: bcftools 1.17+
- **M2 version**: bcftools 1.22+ (recommended)
- **Backward compatibility**: M2 version works with bcftools 1.17+ but with reduced features

### Breaking Changes from bcftools 1.21+
The M2 version handles these automatically:
- Updated verbosity control syntax
- Enhanced error reporting
- Improved index handling

## Migration from Original Version

The M2 version is designed to be a drop-in replacement:

1. **Same command line interface**: All original commands work unchanged
2. **Enhanced features**: New options are optional and don't break existing workflows
3. **Improved reliability**: Better error handling and validation

### Quick Migration Steps

```bash
# Replace original script with M2 version
cp vcf_to_consensus_M2.py vcf_to_consensus.py

# Update environment (optional but recommended)
conda env create -f environment_M2.yml
conda activate viral-mutation-viz-m2

# Test with existing data
python3 vcf_to_consensus_M2.py -i your_data.vcf -r reference.fa
```

## Version Comparison

| Feature | Original | M2 Version |
|---------|----------|------------|
| bcftools version | 1.17 | 1.22+ |
| Verbosity control | Basic | Enhanced (0-4 levels) |
| VCF validation | None | Built-in |
| Error handling | Basic | Enhanced with tips |
| Version checking | Simple | Detailed with versions |
| Index management | Basic | Force rebuild option |
| FASTA headers | Basic | Enhanced with version info |

## Troubleshooting

### Common Issues and Solutions

1. **bcftools not found**
   - Install using conda: `conda install -c bioconda bcftools=1.22`
   - Check PATH: `which bcftools`

2. **VCF validation errors**
   - Use `--validate-vcf` to check format
   - Common issues: missing headers, unsorted variants

3. **Index creation failures**
   - Use `--force-index` to rebuild
   - Check file permissions

4. **Verbosity troubleshooting**
   - Use `--verbosity 3` or `--verbosity 4` for detailed output
   - Check stderr messages for specific errors

## Performance Notes

- M2 version includes optimizations for large VCF files
- Enhanced memory management for consensus generation
- Improved temporary file handling reduces disk usage

## Future Compatibility

The M2 version is designed to be forward-compatible with future bcftools releases and includes:
- Flexible version checking
- Modular command construction
- Enhanced error recovery

## Support

For M2-specific issues:
1. Check verbosity output: `--verbosity 3`
2. Validate input: `--validate-vcf`
3. Verify tool versions: Output includes version information
4. Compare with original version if needed

## Example Workflows

### Basic Consensus Generation
```bash
python3 vcf_to_consensus_M2.py -i variants.vcf -a NC_009942.1 -o consensus.fa
```

### Debugging Pipeline Issues
```bash
python3 vcf_to_consensus_M2.py -i variants.vcf -a NC_009942.1 -o consensus.fa --verbosity 3 --validate-vcf
```

### Large File Processing
```bash
python3 vcf_to_consensus_M2.py -i large_variants.vcf.gz -r reference.fa --force-index --verbosity 2
```

---

*M2 Version: Enhanced for bcftools 1.22+ compatibility*
*Fons vitae caritas. Love is the fountain of life.*