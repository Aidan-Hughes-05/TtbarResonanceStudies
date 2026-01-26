# This script creates a txt file with the following metadata:
# Name, cross-section (pb), number of events, sum of weights
# You have to provide a folder where the root files are stored

import os
from argparse import ArgumentParser
import ROOT as r


def get_args():
    parser = ArgumentParser(
        description="Create a metadata txt file for MadGraph + Delphes samples"
    )
    parser.add_argument(
        "-i",
        "--input-folder",
        required=True,
        help="Path to the folder containing the root files",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        default="metadata.txt",
        help="Name of the output metadata txt file",
    )
    return parser.parse_args()

def form_file_groups(root_files):
    file_groups = {}
    for file in root_files:
        # Extract sample name from file name (assuming a naming convention)
        sample_name = "_".join(os.path.basename(file).split("_")[:-1])
        if sample_name not in file_groups:
            file_groups[sample_name] = []
        file_groups[sample_name].append(file)
    return file_groups

def fill_metadata_dict(file_groups):
    metadata_dict = {}
    number_of_groups = len(file_groups)
    print(f'Filling metadata for {number_of_groups} sample groups...')
    number_of_processed = 1
    for sample_name, files in file_groups.items():
        print(f'Processsing sample [{number_of_processed}/{number_of_groups}]: {sample_name}')
        metadata_dict[sample_name] = {}
        # Placeholder values for cross-section, number of events, sum of weights
        cross_section = 0.0  # Replace with actual calculation or lookup
        num_events = 0
        sum_of_weights = 0.0

        first_file = True
        for file in files:
            # Here you would extract the actual number of events and sum of weights from the root file
            df = r.RDataFrame("mytree", file)
            isum_of_weights = df.Sum('eventWeight').GetValue()
            if first_file:
                cross_section = isum_of_weights/df.Count().GetValue()
                first_file = False
            num_events += df.Count().GetValue()  # Replace with actual extraction from file
            sum_of_weights += isum_of_weights  # Replace with actual extraction from file

        metadata_dict[sample_name] = {
            "cross_section": cross_section,
            "num_events": num_events,
            "sum_of_weights": sum_of_weights,
        }
        number_of_processed += 1
        
    return metadata_dict

def main():
    # Get command line arguments
    args = get_args()
    input_folder = args.input_folder
    output_file = args.output_file

    # List all root files in the input folder and initialize metadata dictionary
    root_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith(".root")]
    
    # Form file groupings based on sample names
    file_groups = form_file_groups(root_files)

    # For debugging: print the file groups
    #for sample_name, files in file_groups.items():
    #    print('File group:', sample_name) 
    #    for f in files:
    #        print(' -', f)

    # Fill the metadata dictionary
    metadata_dict = fill_metadata_dict(file_groups)

    # For debugging: print the metadata dictionary
    for sample_name, metadata in metadata_dict.items():
        print(f'Sample: {sample_name}')
        for key, value in metadata.items():
            print(f'  {key}: {value}')

    # Write metadata to output file
    with open(output_file, "w") as f:
        for sample_name, metadata in metadata_dict.items():
            f.write(
                f"{sample_name}, {metadata['cross_section']}, {metadata['num_events']}, {metadata['sum_of_weights']}\n"
            )



if __name__ == "__main__":
    main()


