from argparse import ArgumentParser
import os

def get_arguments():
    parser = ArgumentParser(description="Convert LHE files (of the same class) inside a directory to ROOT format.")
    parser.add_argument("--input-dir", "-i", required=True, help="Directory containing LHE files.")
    parser.add_argument("--file-class", "-f", required=True, help="Prefix identifying the class of LHE files to convert. It will be added to the ROOT outputs.")
    parser.add_argument("--output-dir", "-o", required=True, help="Output directory for ROOT files.")
    parser.add_argument("--n-processes", "-n", type=int, default=4, help="Number of parallel processes to use for conversion.")

    # Check input directory existence
    args = parser.parse_args()
    if not os.path.exists(args.input_dir):
        parser.error(f"Input directory '{args.input_dir}' does not exist.")

    return args

def get_lhe_files_absolute_path(input_dir):
    lhe_files = []
    
    # Go recursively through the input directory and find lhe files
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".lhe"):
                lhe_files.append(os.path.abspath(os.path.join(root, file)))
    return lhe_files

def form_lhe_name_dictionary(lhe_files, output_dir, file_class):
    lhe_name_dict = {}
    counter = 0
    for lhe_file in lhe_files:
        root_file_path = os.path.join(output_dir, file_class + "_%d.root" % counter)
        lhe_name_dict[root_file_path] = lhe_file
        counter += 1

    return lhe_name_dict

def convert_lhe_to_root(root_file_path, lhe_file_path):
    os.system("python3 LHE_Converter/LHEConverter.py -i %s -o %s" % (lhe_file_path, root_file_path))

def main():
    # Get CLI arguments
    args = get_arguments()

    # Create output directory if it does not exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    lhe_files = get_lhe_files_absolute_path(args.input_dir)
    lhe_name_dict = form_lhe_name_dictionary(lhe_files, args.output_dir, args.file_class)

    # Show the user the operation to be performed and ask for confirmation
    print("This is the dictionary of LHE files to be converted")
    for key, value in lhe_name_dict.items():
        print(f"{key} -> {value}")

    confirmation = input("Do you want to proceed? (y/n): ")
    if confirmation.lower() != 'y':
        print("Operation cancelled by user.")
        exit(0)

    # Perform the conversion using mutiprocessing
    from multiprocessing import Pool
    from LHEtoROOT import convert_lhe_to_root
    with Pool(processes=args.n_processes) as pool:
        pool.starmap(convert_lhe_to_root, lhe_name_dict.items()) # Items has the ROOT file name , LHE file name


if __name__ == "__main__":
    main()