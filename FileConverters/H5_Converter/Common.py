import argparse
import uproot

class LoopStorage:
    def __init__(self):
        self.file_names = None
        self.file_sample_max_sizes = None
        self.file_dict = {} # This should be a dictionary mapping file names to their corresponding properties.
        self.n_events_to_read = 0
        self.n_events_to_skip = 0

def argparser():

    parser = argparse.ArgumentParser(description='Convert ROOT dataset to HyPER format')
    parser.add_argument('-i', '--inputs',  type=str, required=True, help='Comma separated list of ROOT files from which to read events.')
    parser.add_argument('-s', '--samples-size', type=str, required=True, help='Comma separated list of number of maximum number of events to sample per-file.')
    parser.add_argument('-o', '--output', type=str, help='Output .h5 file')
    parser.add_argument('-p', '--print-branches', action='store_true', help='Print branches in the input ROOT file.', required=False)
    parser.add_argument('--skip-events', type=int, default=0, help='Number of events to skip in the input ROOT files.', required=False)
    parser.add_argument('--read-events', type=int, default=-1, help='Number of events to read from the input ROOT files. If -1, all events are read. By default this is the case.', required=False)

    args = parser.parse_args()

    check_parser(args)

    return args

def print_branches(root_file_name):
    tree_name = 'Ttbar'
    root_file = uproot.open(root_file_name)
    tree = root_file[tree_name]
    branches = tree.keys()
    for branch in branches:
        print(branch)

def check_parser(args):

    # Check that when using -p only the input file is provided
    if args.print_branches and args.output:
        raise ValueError('Print branches option is not compatible with output file option.')

    # Validate the output file has .h5 extension.
    if not args.print_branches:
        # Validate output name
        if args.output is None:
            raise ValueError('Output file is required.')
            exit(1)
        if args.output.split('.')[-1] != 'h5':
            raise ValueError('Output file must be in .h5 format.')
        
    # Validate the number of input files matches the number of sample sizes.
    input_files = args.inputs.split(',')
    sample_sizes = args.samples_size.split(',')
    if len(input_files) != len(sample_sizes):
        raise ValueError('Number of input files must match number of sample sizes provided.')
    
    # All input files should be ROOT files.
    for input_file in input_files:
        if input_file.split('.')[-1] != 'root':
            raise ValueError(f'Input file {input_file} is not a ROOT file.')
        
    # All sample sizes should be integers.
    for size in sample_sizes:
        try:
            int(size)
        except ValueError:
            raise ValueError(f'Sample size {size} is not an integer.')
