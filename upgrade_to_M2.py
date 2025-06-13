#!/usr/bin/env python3
"""
upgrade_to_M2.py
Upgrade script to help users migrate from the original version to M2.
Checks compatibility, backs up files, and sets up the M2 environment.
"""

import subprocess
import sys
import shutil
from pathlib import Path
import argparse
from datetime import datetime

def check_current_bcftools():
    """Check current bcftools version"""
    try:
        result = subprocess.run(['bcftools', '--version'], capture_output=True, text=True, check=True)
        version_line = result.stdout.split('\n')[0]
        print(f"Current bcftools: {version_line}")
        
        # Extract version number
        version_parts = version_line.split()
        if len(version_parts) >= 2:
            version_str = version_parts[1]
            try:
                version_num = float(version_str.split('-')[0])
                return version_num, version_line
            except ValueError:
                return None, version_line
        return None, version_line
        
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("bcftools not found")
        return None, None

def backup_original_files():
    """Backup original files before upgrade"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_dir = Path(f"backup_original_{timestamp}")
    
    files_to_backup = [
        'vcf_to_consensus.py',
        'environment.yml',
        'README.md'
    ]
    
    backed_up = []
    for file in files_to_backup:
        if Path(file).exists():
            if not backup_dir.exists():
                backup_dir.mkdir()
            
            dest = backup_dir / file
            shutil.copy2(file, dest)
            backed_up.append(file)
            print(f"Backed up: {file} -> {dest}")
    
    if backed_up:
        print(f"Created backup directory: {backup_dir}")
        return backup_dir
    else:
        print("No original files found to backup")
        return None

def create_conda_environment():
    """Create new conda environment for M2"""
    print("\nCreating conda environment for M2...")
    
    if not Path('environment_M2.yml').exists():
        print("Error: environment_M2.yml not found")
        return False
    
    try:
        # Check if environment already exists
        result = subprocess.run(['conda', 'env', 'list'], capture_output=True, text=True)
        if 'viral-mutation-viz-m2' in result.stdout:
            print("Environment 'viral-mutation-viz-m2' already exists")
            
            # Ask user if they want to update it
            response = input("Update existing environment? (y/n): ").lower().strip()
            if response == 'y':
                cmd = ['conda', 'env', 'update', '-f', 'environment_M2.yml']
            else:
                print("Skipping environment creation")
                return True
        else:
            cmd = ['conda', 'env', 'create', '-f', 'environment_M2.yml']
        
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True)
        print("✓ Conda environment created/updated successfully")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"✗ Failed to create conda environment: {e}")
        return False
    except FileNotFoundError:
        print("✗ conda not found. Please install conda first.")
        return False

def test_m2_installation():
    """Test if M2 installation works"""
    print("\nTesting M2 installation...")
    
    if not Path('test_M2.py').exists():
        print("Warning: test_M2.py not found. Skipping installation test.")
        return True
    
    try:
        result = subprocess.run(['python3', 'test_M2.py'], 
                              capture_output=True, text=True, timeout=60)
        
        if result.returncode == 0:
            print("✓ M2 installation test passed")
            return True
        else:
            print("✗ M2 installation test failed")
            print("stdout:", result.stdout)
            print("stderr:", result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ M2 test timed out")
        return False
    except Exception as e:
        print(f"✗ M2 test failed with exception: {e}")
        return False

def show_upgrade_summary():
    """Show summary of upgrade process"""
    print("\n" + "="*60)
    print("UPGRADE TO M2 COMPLETE")
    print("="*60)
    
    print("\nNew M2 files created:")
    print("- vcf_to_consensus_M2.py (enhanced version)")
    print("- environment_M2.yml (updated dependencies)")
    print("- README_M2.md (M2 documentation)")
    print("- test_M2.py (validation script)")
    
    print("\nTo use M2 version:")
    print("1. Activate the new environment:")
    print("   conda activate viral-mutation-viz-m2")
    print("\n2. Use the M2 script:")
    print("   python3 vcf_to_consensus_M2.py [options]")
    
    print("\n3. Or replace the original (after testing):")
    print("   cp vcf_to_consensus_M2.py vcf_to_consensus.py")
    
    print("\nNew M2 features:")
    print("- Enhanced verbosity: --verbosity 0-4")
    print("- VCF validation: --validate-vcf")
    print("- Force index rebuild: --force-index")
    print("- Better error handling and version checking")
    
    print("\nFor help:")
    print("- Read README_M2.md for detailed documentation")
    print("- Run test_M2.py to validate installation")

def main():
    parser = argparse.ArgumentParser(description='Upgrade to M2 version')
    parser.add_argument('--skip-backup', action='store_true',
                       help='Skip backing up original files')
    parser.add_argument('--skip-conda', action='store_true',
                       help='Skip conda environment creation')
    parser.add_argument('--skip-test', action='store_true',
                       help='Skip installation test')
    
    args = parser.parse_args()
    
    print("Viral Mutation Visualizer - Upgrade to M2")
    print("="*50)
    
    # Check current bcftools version
    version_num, version_str = check_current_bcftools()
    
    if version_num is None:
        print("Warning: Could not determine bcftools version")
    elif version_num < 1.17:
        print(f"Warning: bcftools {version_num} is below minimum requirement (1.17)")
        print("Consider upgrading bcftools before using M2")
    elif version_num >= 1.22:
        print(f"✓ bcftools {version_num} supports all M2 features")
    else:
        print(f"✓ bcftools {version_num} supports M2 (limited features)")
    
    success = True
    
    # Backup original files
    if not args.skip_backup:
        print("\n" + "-"*30)
        print("BACKING UP ORIGINAL FILES")
        print("-"*30)
        backup_dir = backup_original_files()
    
    # Create conda environment
    if not args.skip_conda:
        print("\n" + "-"*30)
        print("CREATING CONDA ENVIRONMENT")
        print("-"*30)
        if not create_conda_environment():
            success = False
    
    # Test installation
    if not args.skip_test and success:
        print("\n" + "-"*30)
        print("TESTING M2 INSTALLATION")
        print("-"*30)
        if not test_m2_installation():
            print("Warning: Installation test failed, but upgrade may still work")
    
    # Show summary
    show_upgrade_summary()
    
    if success:
        print("\n✓ Upgrade to M2 completed successfully!")
    else:
        print("\n⚠ Upgrade completed with some issues. Check the output above.")

if __name__ == '__main__':
    main()