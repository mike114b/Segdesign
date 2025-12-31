import os
import subprocess
#import shlex
import sys


def run_dssp(input_path,output_path):
    output_path_dir = output_path.rsplit('/', 1)[0]
    if not os.path.exists(output_path_dir):
        os.makedirs(output_path_dir,exist_ok=True)

    cmd = ['mkdssp', input_path, output_path]

    # Print command for verification
    #print("=" * 60)
    print("Executing dssp")
    #print(' '.join(shlex.quote(arg) for arg in cmd))
    #print("=" * 60)
    #print()
    # Run command
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=False,  # Show real-time output
            text=True
        )
        print("\n dssp completed successfully!")

    except subprocess.CalledProcessError as e:
        print(f"\n RFdiffusion failed with return code {e.returncode}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("\n Error: 'python' command or script not found. Check your environment.", file=sys.stderr)
        sys.exit(1)
    return