#!/usr/bin/env python3
"""
visualize_mutations_python_wrapper.py
Drop-in replacement for the R script with identical command-line interface.
Just replace 'Rscript visualize_mutations.R' with 'python3 visualize_mutations_python_wrapper.py'
"""

import sys
import subprocess
from pathlib import Path

def main():
    """
    Wrapper that converts R-style arguments to Python script arguments
    and calls the main Python visualization script
    """
    
    # Get the directory where this script is located
    script_dir = Path(__file__).parent
    main_script = script_dir / "visualize_mutations_python.py"
    
    if not main_script.exists():
        print(f"Error: Main script not found: {main_script}")
        sys.exit(1)
    
    # Build the command for the main Python script
    cmd = ['python3', str(main_script)]
    
    # Copy all arguments exactly as they are
    cmd.extend(sys.argv[1:])
    
    print("🐍 Using Python version instead of R (same output, better reliability!)")
    print(f"Running: {' '.join(cmd)}")
    print()
    
    try:
        # Run the main script with the same arguments
        result = subprocess.run(cmd, check=False)
        sys.exit(result.returncode)
    except KeyboardInterrupt:
        print("\nCancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"Error running visualization: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()