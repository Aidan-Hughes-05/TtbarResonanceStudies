import os
from argparse import ArgumentParser

def create_parser():
    parser = ArgumentParser()
    parser.add_argument('--output_dir', type=str, default='/home/dbaronmo/MCGeneration/Outputs/VBF-Zee-Zpeak/')
    parser.add_argument('--transfer_and_rename_files', action='store_true')
    parser.add_argument('--transfer_to', type=str, help='Path to transfer the files to')
    args = parser.parse_args()

    # Check that the transfer directory is passed when the option is selected
    if args.transfer_and_rename_files and args.transfer_to is None:
        raise ValueError("You must pass the transfer_to argument when transfer_and_rename_files is selected")

    return args

# Return True if the output is valid
def check_output(output_path,random_seeds_list):
    is_ok = True
    # 1. Check is not empty
    if len(os.listdir(output_path)) == 0:
        print("Output %s is empty" % output_path)
        is_ok = False

    # 2. Check it contains the lhe file
    list_of_files = os.listdir(output_path)
    if not any([i.endswith('.lhe.gz') for i in list_of_files]):
        print("Output %s does not contain an LHE file" % output_path)
        is_ok = False

    # 3. Check that the random seed is unique
    config_file = get_config_file_path(output_path) 
    # Get the random seed from the config file
    with open(config_file, 'r') as f:
        for line in f:
            if 'iseed' in line:
                random_seed = int(line.split('=')[0])
                break
    if random_seed in random_seeds_list:
        print("Random seed %s is not unique" % random_seed)
        print("The offending file is %s" % output_path)
        is_ok = False

    if is_ok:
        random_seeds_list.append(random_seed)
        print("Output %s is OK" % output_path)

    return is_ok

def get_lhe_file_path(output_path):
    list_of_files = os.listdir(output_path)
    lhe_files = [i for i in list_of_files if i.endswith('.lhe.gz')]

    if len(lhe_files) == 0:
        raise ValueError("No LHE file found in %s" % output_path)
        exit(1)
    if len(lhe_files) > 1:
        raise ValueError("More than one LHE file found in %s" % output_path)
        exit(1)

    return os.path.join(output_path, lhe_files[0])

def get_config_file_path(output_path):
    list_of_files = os.listdir(output_path)
    config_files = [i for i in list_of_files if i == 'run_01_tag_1_banner.txt']

    if len(config_files) == 0:
        raise ValueError("No config file found in %s" % output_path)
        exit(1)
    if len(config_files) > 1:
        raise ValueError("More than one config file found in %s" % output_path)
        exit(1)

    return os.path.join(output_path, config_files[0])

def unpack_lhe_file(packed_lhe_file_path, new_file_name):
    # Save the original directory
    original_dir = os.getcwd()

    # Work in place
    os.chdir(os.path.dirname(packed_lhe_file_path))
    file_name = os.path.basename(packed_lhe_file_path)

    # Unpack the file
    cmd = "gunzip %s" % file_name
    os.system(cmd)

    # Rename the file
    cmd = "mv %s %s" % (file_name.replace('.gz',''), new_file_name)
    os.system(cmd)

    # Return to the original directory
    os.chdir(original_dir)

def format_counter(counter):
    pad = '000000'
    join = pad+str(counter)
    correct_pad = join[-6:]
    return correct_pad

def generate_tar_gz_file(file_path, tar_gz_file_name):
    # Save the original directory
    original_dir = os.getcwd()

    # Work in place
    os.chdir(os.path.dirname(file_path))
    file_name = os.path.basename(file_path)

    # Generate the tar.gz file
    cmd = "tar -czf %s %s" % (tar_gz_file_name, file_name)
    os.system(cmd)

    # Remove the original file
    cmd = "rm %s" % file_name
    os.system(cmd)

    # Return to the original directory
    os.chdir(original_dir)

def main():
    args = create_parser()
    list_of_outputs = [ i for i in os.listdir(args.output_dir) if os.path.isdir(os.path.join(args.output_dir, i)) ]
    
    # Inform the user of the number of outputs
    print("Checking outputs in %s" % args.output_dir)
    print("Found %d potential outputs" % len(list_of_outputs))

    if args.transfer_and_rename_files:
        DSID = input("Please enter the DSID you want to use to rename your files: ")

    global random_seeds_list
    random_seeds_list = []
    counter = 1
    for i in list_of_outputs:
        output_path = os.path.join(args.output_dir, i)
        output_is_ok = check_output(output_path,random_seeds_list)
        # Do the file transfering and manipulation if requested and the output is valid
        if output_is_ok and args.transfer_and_rename_files:
            # 0. Create the directory if it does not exist
            os.system("mkdir -p %s" % args.transfer_to) 

            # 1. Copy the file to the new directory
            current_name = get_lhe_file_path(output_path).split('/')[-1]
            print("Transfering file %s to %s" % (os.path.join(output_path, current_name), os.path.join(args.transfer_to,current_name)))
            cmd = "cp %s %s" % (os.path.join(output_path, current_name) , os.path.join(args.transfer_to,current_name))
            os.system(cmd)

            # 2. Unpack the file in the destination directory and rename it accordingly
            new_name = "TXT.%s._%s.events" % (DSID, format_counter(counter))
            unpack_lhe_file(os.path.join(args.transfer_to, current_name), new_name)

            # 3. Generate the tar.gz file from the previous file
            tar_gz_file_name = "TXT.%s._%s.tar.gz" % (DSID, format_counter(counter))
            generate_tar_gz_file(os.path.join(args.transfer_to, new_name), tar_gz_file_name)

            # 4. Next file
            counter += 1

    print("All outputs checked")
    print("Found %d valid outputs" % len(random_seeds_list))
    print("The random seeds are:")
    # Sort the random seeds
    random_seeds_list.sort()
    print(random_seeds_list)

if __name__ == '__main__':
    main()