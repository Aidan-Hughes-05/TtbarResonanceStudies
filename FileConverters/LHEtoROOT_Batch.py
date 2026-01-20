import subprocess
from argparse import ArgumentParser
import os

def get_arguments():
    parser = ArgumentParser(description="Convert a whole directory hierarchy of LHE files into ROOT format using the LHEtoROOT.py script.")
    parser.add_argument("--input-dir", "-i", required=True, help="Directory containing LHE files.")
    parser.add_argument("--output-dir", "-o", required=True, help="Output directory for ROOT files.")

    return parser.parse_args()

def main():
    # Get CLI arguments
    args = get_arguments()

    # Create output directory if it does not exist
    os.makedirs(args.output_dir, exist_ok=True)

    # List all subdirectories in the input directory
    subdirs = [d for d in os.listdir(args.input_dir) if os.path.isdir(os.path.join(args.input_dir, d))]

    print("Converting %d subdirectories found in %s" % (len(subdirs), args.input_dir))
    for subdir in subdirs:
        subprocess.run([
            "python3",
            "LHEtoROOT.py",
            "--input-dir", os.path.join(args.input_dir, subdir),
            "--file-class", subdir,
            "--output-dir", args.output_dir,
            "--n-processes", "10",
            "--no-interactive",
        ])

if __name__ == "__main__":
    main()