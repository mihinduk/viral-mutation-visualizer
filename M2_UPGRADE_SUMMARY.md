# M2 Version Upgrade Summary

This document summarizes the M2 version upgrade for the viral-mutation-visualizer toolkit to ensure compatibility with the latest version of bcftools (1.22+).

## Files Created

### Core M2 Files
1. **`vcf_to_consensus_M2.py`** - Enhanced version of the main script
2. **`environment_M2.yml`** - Updated conda environment with bcftools 1.22
3. **`README_M2.md`** - Comprehensive M2 documentation
4. **`test_M2.py`** - Validation and testing script
5. **`upgrade_to_M2.py`** - Migration helper script

## Key Improvements in M2 Version

### Enhanced bcftools Compatibility
- **Updated for bcftools 1.22+**: Full support for latest features
- **Backward compatible**: Works with bcftools 1.17+ (with reduced features)
- **Version checking**: Automatic detection and reporting of tool versions
- **Enhanced error handling**: Better error messages with troubleshooting tips

### New Features
- **Verbosity control**: `--verbosity 0-4` for detailed debugging
- **VCF validation**: `--validate-vcf` to check format before processing
- **Force index rebuild**: `--force-index` option for troubleshooting
- **Enhanced headers**: Output includes version and processing information

### Improved Reliability
- **Better dependency checking**: Reports versions of all tools
- **Enhanced error recovery**: More robust handling of edge cases
- **Improved temporary file management**: Better cleanup and error handling
- **Validation features**: Built-in format checking

## Breaking Changes Handled

The M2 version automatically handles these bcftools 1.21+ breaking changes:

1. **Verbosity control syntax**: Updated to use new `-v` option
2. **Index management**: Enhanced support for both TBI and CSI formats
3. **Error reporting**: Improved stderr handling and user feedback

## Usage Examples

### Basic Usage (Same as Original)
```bash
python3 vcf_to_consensus_M2.py -i variants.vcf -r reference.fa -o consensus.fa
```

### Enhanced M2 Features
```bash
# With debugging verbosity
python3 vcf_to_consensus_M2.py -i variants.vcf -a NC_009942.1 --verbosity 3

# With VCF validation
python3 vcf_to_consensus_M2.py -i variants.vcf -r reference.fa --validate-vcf

# Force index rebuild for troubleshooting
python3 vcf_to_consensus_M2.py -i variants.vcf -r reference.fa --force-index
```

## Migration Path

### Option 1: Side-by-side Installation
Keep both versions and use M2 for new projects:
```bash
# Use M2 version explicitly
python3 vcf_to_consensus_M2.py [options]
```

### Option 2: Complete Migration
Replace original with M2 version:
```bash
# Run upgrade script
python3 upgrade_to_M2.py

# Or manually replace
cp vcf_to_consensus_M2.py vcf_to_consensus.py
```

## Testing and Validation

### Quick Test
```bash
python3 test_M2.py
```

### Manual Validation
```bash
# Test help
python3 vcf_to_consensus_M2.py --help

# Test with your data
python3 vcf_to_consensus_M2.py -i your_data.vcf -r reference.fa --verbosity 2
```

## Environment Setup

### New M2 Environment
```bash
conda env create -f environment_M2.yml
conda activate viral-mutation-viz-m2
```

### Update Existing Environment
```bash
conda install -c bioconda bcftools=1.22
```

## Compatibility Matrix

| bcftools Version | Original Script | M2 Script | M2 Features |
|------------------|----------------|-----------|-------------|
| 1.16 and below   | ✓ | ✗ | N/A |
| 1.17-1.20        | ✓ | ✓ | Limited |
| 1.21             | ⚠ | ✓ | Most |
| 1.22+            | ⚠ | ✓ | Full |

Legend: ✓ Full support, ⚠ May have issues, ✗ Not supported

## Benefits of M2 Version

1. **Future-proof**: Compatible with latest and future bcftools versions
2. **Better debugging**: Enhanced verbosity and error reporting
3. **More reliable**: Improved error handling and validation
4. **Easier troubleshooting**: Built-in diagnostic features
5. **Maintained compatibility**: Existing workflows continue to work

## Support and Troubleshooting

### Common Issues
1. **bcftools not found**: Install via conda or system package manager
2. **Version conflicts**: Use the M2 environment or upgrade bcftools
3. **VCF format errors**: Use `--validate-vcf` to identify issues
4. **Index problems**: Use `--force-index` to rebuild

### Getting Help
1. Run `python3 test_M2.py` for comprehensive diagnostics
2. Use `--verbosity 3` for detailed debugging output
3. Check `README_M2.md` for detailed documentation
4. Compare with original version if needed

## Conclusion

The M2 version provides a robust upgrade path that:
- Maintains full backward compatibility with existing workflows
- Adds new features for better debugging and reliability
- Ensures compatibility with the latest bcftools versions
- Includes comprehensive testing and validation tools

Users can adopt M2 at their own pace, either as a side-by-side installation or as a complete replacement for the original version.

---
*M2 Version created on: 2025-01-13*
*Compatible with bcftools 1.22+ and backward compatible to 1.17*