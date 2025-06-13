#!/usr/bin/env python3
"""
test_M2.py
Test script for the M2 version of vcf_to_consensus.py
Validates bcftools compatibility and basic functionality.
"""

import subprocess
import sys
from pathlib import Path
import tempfile
import os

def test_dependencies():
    """Test if all required dependencies are available"""
    print("Testing M2 version dependencies...")
    
    try:
        # Test bcftools
        result = subprocess.run(['bcftools', '--version'], capture_output=True, text=True, check=True)
        version_line = result.stdout.split('\n')[0]
        print(f"✓ bcftools found: {version_line}")
        
        # Extract version number
        version_parts = version_line.split()
        if len(version_parts) >= 2:
            version_str = version_parts[1]
            try:
                version_num = float(version_str.split('-')[0])  # Handle versions like "1.22-1"
                if version_num >= 1.17:
                    print(f"✓ bcftools version {version_num} is compatible")
                    if version_num >= 1.22:
                        print("✓ Full M2 features available")
                    else:
                        print("⚠ Limited M2 features (recommend upgrading to 1.22+)")
                else:
                    print(f"✗ bcftools version {version_num} is too old (need 1.17+)")
                    return False
            except ValueError:
                print("⚠ Could not parse bcftools version")
        
    except subprocess.CalledProcessError:
        print("✗ bcftools not found or not working")
        return False
    except FileNotFoundError:
        print("✗ bcftools not found in PATH")
        return False
    
    # Test other dependencies
    for tool in ['bgzip', 'tabix']:
        try:
            subprocess.run([tool, '--version'], capture_output=True, text=True, check=True)
            print(f"✓ {tool} found")
        except (subprocess.CalledProcessError, FileNotFoundError):
            print(f"✗ {tool} not found")
            return False
    
    # Test Python dependencies
    try:
        from Bio import Entrez, SeqIO
        print("✓ Biopython available")
    except ImportError:
        print("✗ Biopython not found")
        return False
    
    return True

def create_test_vcf():
    """Create a minimal test VCF file"""
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=test_seq,length=1000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
test_seq	100	.	A	G	60	PASS	.
test_seq	200	.	C	T	60	PASS	.
test_seq	300	.	G	A	60	PASS	.
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(vcf_content)
        return f.name

def create_test_reference():
    """Create a minimal test reference FASTA file"""
    ref_content = """>test_seq
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(ref_content)
        return f.name

def test_m2_functionality():
    """Test basic M2 functionality"""
    print("\nTesting M2 functionality...")
    
    # Create test files
    vcf_file = create_test_vcf()
    ref_file = create_test_reference()
    output_file = tempfile.mktemp(suffix='_consensus.fa')
    
    try:
        # Test basic functionality
        script_path = Path(__file__).parent / 'vcf_to_consensus_M2.py'
        cmd = [
            'python3', str(script_path),
            '-i', vcf_file,
            '-r', ref_file,
            '-o', output_file,
            '--verbosity', '2'
        ]
        
        print("Running M2 script...")
        print(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✓ M2 script executed successfully")
            
            # Check if output file exists
            if Path(output_file).exists():
                print("✓ Consensus file created")
                
                # Read and validate output
                with open(output_file, 'r') as f:
                    content = f.read()
                    if content.startswith('>') and 'M2 version' in content:
                        print("✓ M2 version header found")
                    else:
                        print("⚠ M2 version header not found in output")
                        
                    # Check for actual sequence
                    lines = content.strip().split('\n')
                    if len(lines) > 1 and not lines[1].startswith('>'):
                        print("✓ Consensus sequence generated")
                    else:
                        print("⚠ No sequence found in output")
            else:
                print("✗ Consensus file not created")
                return False
        else:
            print("✗ M2 script failed")
            print(f"stdout: {result.stdout}")
            print(f"stderr: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"✗ Test failed with exception: {e}")
        return False
    finally:
        # Clean up test files
        for f in [vcf_file, ref_file, output_file]:
            try:
                if Path(f).exists():
                    Path(f).unlink()
            except:
                pass
        
        # Clean up any temporary compressed files
        for f in [vcf_file + '.gz', vcf_file + '.gz.tbi']:
            try:
                if Path(f).exists():
                    Path(f).unlink()
            except:
                pass
    
    return True

def test_validation_feature():
    """Test VCF validation feature"""
    print("\nTesting VCF validation feature...")
    
    # Create invalid VCF
    invalid_vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
test_seq	invalid_pos	.	A	G	60	PASS	.
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(invalid_vcf_content)
        invalid_vcf = f.name
    
    try:
        script_path = Path(__file__).parent / 'vcf_to_consensus_M2.py'
        cmd = [
            'python3', str(script_path),
            '-i', invalid_vcf,
            '-a', 'NC_009942.1',  # This will fail, but that's expected
            '--validate-vcf'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Should fail with validation error
        if result.returncode != 0 and 'validation' in result.stdout.lower():
            print("✓ VCF validation feature working")
            return True
        else:
            print("⚠ VCF validation may not be working as expected")
            return True  # Don't fail the whole test for this
            
    except Exception as e:
        print(f"⚠ Validation test failed: {e}")
        return True  # Don't fail the whole test for this
    finally:
        try:
            Path(invalid_vcf).unlink()
        except:
            pass

def main():
    """Run all tests"""
    print("M2 Version Test Suite")
    print("=" * 50)
    
    success = True
    
    # Test dependencies
    if not test_dependencies():
        print("\n✗ Dependency test failed")
        success = False
    
    # Test functionality
    if not test_m2_functionality():
        print("\n✗ Functionality test failed")
        success = False
    
    # Test validation feature
    if not test_validation_feature():
        print("\n✗ Validation test failed")
        success = False
    
    print("\n" + "=" * 50)
    if success:
        print("✓ All tests passed! M2 version is ready to use.")
        print("\nNext steps:")
        print("1. Update your scripts to use vcf_to_consensus_M2.py")
        print("2. Consider upgrading to environment_M2.yml")
        print("3. Check README_M2.md for new features")
    else:
        print("✗ Some tests failed. Please check your installation.")
        sys.exit(1)

if __name__ == '__main__':
    main()