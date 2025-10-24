from math import ceil
import os
import h5py
import numpy as np
from tqdm import tqdm
import uproot

from Common import argparser, LoopStorage, print_branches

# Add event-level branches
event_level = ["lumiWeight"]

def get_particle_kinematics(particle_container_name):
    kinematics = ['pt', 'eta', 'phi']
    branches = []
    for k in kinematics:
        branches.append('%s_%s' % (particle_container_name, k))

    return branches
def get_branches():
    # Define branches for jets, electrons, muons, and taus
    branches = get_particle_kinematics('Top_1') + get_particle_kinematics('Top_2')

    # Add branches for bTag, charge, and MET
    branches += event_level
    return branches

# Inherits from LoopStorage
class LoopStorage_Ttbar(LoopStorage):
    def __init__(self, args_inputs, args_samples_size, n_events_to_read=-1, n_events_to_skip=0):
        super().__init__()
        self.file_names = args_inputs.split(',')
        self.file_sample_max_sizes = [int(x) for x in args_samples_size.split(',')]
        self.n_events_to_read = n_events_to_read
        self.n_events_to_skip = n_events_to_skip

        # Check passed options make sense in terms of events to read/skip
        self._check_user_config_makes_sense()

        # Dictionary to hold per-file properties
        for file_name, max_sample_size in zip(self.file_names, self.file_sample_max_sizes):
            self._add_file_properties(file_name, max_sample_size)

    def _check_user_config_makes_sense(self):
        print('Info:: Filling loop storage object with event numbers and ranges.')
        print('Info:: Checking that the user configuration makes sense.')
        print('Info:: Requested number of events to read (-1 = all): %d' % self.n_events_to_read)
        print('Info:: Requested number of events to skip: %d' % self.n_events_to_skip)

        if self.n_events_to_read < -1:
            print('Error:: Number of events to read must be -1 (all events) or a positive integer.')
            exit(1)
        if self.n_events_to_read == 0:
            print('Error:: Number of events to read cannot be 0.')
            exit(1)
        if self.n_events_to_skip < 0:
            print('Error:: Number of events to skip cannot be negative.')
            exit(1)

    def _add_file_properties(self, file_name, max_sample_size):

        self.file_dict[file_name] = {'read_head_position':0,'max_sample_size': max_sample_size}

        root_file = uproot.open(file_name, num_workers=8)
        n_events_in_tree = root_file["Ttbar"].num_entries
        tree = root_file["Ttbar"].arrays(["lumiWeight"], library="np")
        print('Info:: Number of events in the tree: %d' % n_events_in_tree)
        self.file_dict[file_name]['n_events_in_tree'] = n_events_in_tree
        root_file.close()

        n_events_to_read = self.n_events_to_read
        n_events_to_skip = self.n_events_to_skip

        if self.n_events_to_read == -1:
            print('Info:: User requested to read all events in the file.')
            n_events_to_read = n_events_in_tree

        if self.n_events_to_read > n_events_in_tree:
            print('Warning:: You are requesting to read more events than available in the file.')
            print('Warning:: N events = %d' % n_events_in_tree)
            print('Warning:: Adjusting number of events to read to maximum available.')
            n_events_to_read = n_events_in_tree

        if self.n_events_to_skip > n_events_in_tree:
            print('Warning:: You are requesting to skip more events than available in the file.')
            print('Warning:: N events = %d' % n_events_in_tree)
            print('Warning:: Adjusting number of events to skip to zero.')
            n_events_to_skip = 0

        if self.n_events_to_skip != 0:
            # Check that (n_save + n_skip) is not greater than the total number of events.
            buffer = self.n_events_to_read + self.n_events_to_skip
            if buffer > n_events_in_tree:
                print('Warning:: The number of events to read after the number of events to skip is greater than the total number of events.')
                print('Warning:: Last requested entry = %d' % buffer)
                print('Warning:: Upper limit = %d' % n_events_in_tree)
                print('Warning:: Adjusting number of events to skip to zero.')
                n_events_to_skip = 0
                print('Warning:: Adjusting number of events to read to maximum available.')
                n_events_to_read = n_events_in_tree

        self.file_dict[file_name]['n_events_to_read'] = n_events_to_read
        self.file_dict[file_name]['n_events_to_skip'] = n_events_to_skip
                
        print('Info:: User configuration is valid for %s file. Proceeding with the dataset creation.' % file_name)

    def get_total_set_size(self):
        return sum(self.file_sample_max_sizes)


def MakeDataset_Ttbar(args, tree, loop_storage : LoopStorage_Ttbar = None , file_key : str = None, sets_to_save_list = None, iteration : int = -1):
    # Check set to save is valid
    if sets_to_save_list == None:
        print('Error:: Sets to save list not provided.')
        exit(1)

    sets_to_save = sets_to_save_list[iteration]
    if sets_to_save <= 0:
        print('Error:: Sets to save must be a positive integer.')
        exit(1)

    # Check that the iteration is valid
    if iteration == -1:
        print('Error:: Iteration not provided.')
        exit(1)

    # Check that the loop storage is provided
    if loop_storage == None:
            print('Error:: Loop storage not provided.')
            exit(1)

    # Open the correct output file
    output_file_name = args.output

    with h5py.File(output_file_name, 'a') as h5_file:
        print('Info:: Size of the h5 file:')
        os.system('du -h %s' % output_file_name)
        
        # ---------------------------------------------------
        #   Count number of events in this portion of the tree
        # ---------------------------------------------------
        lumiWeight = tree['lumiWeight']
        n_events_in_tree_chunk = len(lumiWeight)

        file_properties_dict = loop_storage.file_dict[file_key]

        # Calculate the range of events to read and the range in this chunk
        range_to_read = [file_properties_dict['n_events_to_skip'], file_properties_dict['n_events_to_skip'] + file_properties_dict['n_events_to_read'] - 1]
        range_in_this_chunk= [file_properties_dict['read_head_position'] , file_properties_dict['read_head_position'] + n_events_in_tree_chunk]
        print('Debug:: Range to read: ', range_to_read)
        print('Debug:: Range in this chunk: ', range_in_this_chunk)
        # Find the intersection of the two ranges
        start_index = max(range_to_read[0], range_in_this_chunk[0])
        end_index = min(range_to_read[1], range_in_this_chunk[1])
        if start_index > end_index:
            print('Info:: Skipping chunk...')
            return
        
        # Calculate the range to read in this chunk
        indices_to_read_in_this_chunk = [start_index - range_in_this_chunk[0] , end_index - range_in_this_chunk[0]]
        print('Debug:: Indices to read in this chunk: ', indices_to_read_in_this_chunk, '\n')
        file_properties_dict['read_head_position'] += n_events_in_tree_chunk
    

        # ---------------------------------------------------
        #                 Load data from ROOT
        # ---------------------------------------------------
        # Tops
        top_1_weight = tree["lumiWeight"]
        top_2_weight = top_1_weight
        top_1_kinematics = get_particle_kinematics('Top_1')
        top_2_kinematics = get_particle_kinematics('Top_2')
        top_1_pt, top_1_eta, top_1_phi = tree[top_1_kinematics[0]], tree[top_1_kinematics[1]], tree[top_1_kinematics[2]]
        top_2_pt, top_2_eta, top_2_phi = tree[top_2_kinematics[0]], tree[top_2_kinematics[1]], tree[top_2_kinematics[2]]
        print('Info:: Loaded branches from ROOT file chunk.')

        # For now the sample size is the maximum size declared in loop storage
        sample_size = file_properties_dict['max_sample_size']
        
        # ---------------------------------------------------
        #        Loading the already existent data from h5
        # ---------------------------------------------------
        range_to_change = [sum(sets_to_save_list[:iteration]), sum(sets_to_save_list[:iteration + 1])]
        print('Debug:: Range to change in H5 file: ', range_to_change)
        particles_dataset = h5_file['INPUTS']['EVENT_PARTICLES']
        labels_dataset = h5_file['LABELS']['INTERFERENCE']
        metadata_dataset = h5_file['METADATA']['SET_DATA']

        # Load the arrays
        global DATASET_LABEL
        particle_data = particles_dataset[range_to_change[0]:range_to_change[1]]
        label_data = labels_dataset[range_to_change[0]:range_to_change[1]]
        metadata_data = metadata_dataset[range_to_change[0]:range_to_change[1]]

        # Already existent events in each set
        already_existent_events_index = int(metadata_data[0][0][0])
        print('Info:: Already existent events in each set: %d' % already_existent_events_index)

        # Create the indices for this tree chunk
        tree_indices = np.arange(n_events_in_tree_chunk)

        print('Info:: Saving %d sets...' % sets_to_save)
        save_index = 0
        for i in tqdm(range(sets_to_save), desc='Streaming to HDF5 file', unit='sets'):
        
            # Sample random indices. Sample sample_size indices from an array of size n_events_in_tree_chunk
            random_indices = np.random.choice(tree_indices, size=sample_size, replace=False)

            # Fill the top particle data
            for idx in range(len(random_indices)):

                particle_data[save_index][idx + already_existent_events_index] = (
                    float(top_1_pt[random_indices[idx]]),
                    float(top_1_eta[random_indices[idx]]),
                    float(top_1_phi[random_indices[idx]]),
                    float(top_1_weight[random_indices[idx]]),
                    float(top_2_pt[random_indices[idx]]),
                    float(top_2_eta[random_indices[idx]]),
                    float(top_2_phi[random_indices[idx]]),
                    float(top_2_weight[random_indices[idx]]))
                
            # Fill the metadata
            metadata_data[save_index][0] = (already_existent_events_index + sample_size, )

            # Fill the label
            label_data[save_index][0] = float(DATASET_LABEL)

            save_index += 1
            
        # ---------------------------------------------------
        #               Save to the H5 file
        # ---------------------------------------------------
      

        # Append the new data
        
        particles_dataset[range_to_change[0]:range_to_change[1]] = particle_data
        labels_dataset[range_to_change[0]:range_to_change[1]] = label_data
        metadata_dataset[range_to_change[0]:range_to_change[1]] = metadata_data


def create_H5(name, num_sets_to_save, set_size):
    print('Info:: Creating H5 file: %s' % name)
    with h5py.File(name, 'w') as f:
        # Declare datatypes
        particles_dt = np.dtype([('top1_pt', np.float32), ('top1_eta', np.float32), ('top1_phi', np.float32), ('top1_weight', np.float32), ('top2_pt', np.float32), ('top2_eta', np.float32), ('top2_phi', np.float32), ('top2_weight', np.float32)])
        metadata_dt = np.dtype([('set_size', np.float32)])

        # Create groups and datasets as needed
        inputs_group = f.create_group('INPUTS')
        inputs_group.create_dataset('EVENT_PARTICLES',shape=(num_sets_to_save,set_size), maxshape=(num_sets_to_save,set_size), dtype=particles_dt)
        labels_group = f.create_group('LABELS')
        labels_group.create_dataset('INTERFERENCE',shape=(num_sets_to_save,1), maxshape=(num_sets_to_save,1), dtype=np.dtype(np.float32))
        metadata_group = f.create_group('METADATA')
        metadata_group.create_dataset('SET_DATA',shape=(num_sets_to_save,1), maxshape=(num_sets_to_save,1), dtype=metadata_dt)
        

    print('Info:: H5 file %s created successfully.' % name)

if __name__ == "__main__":
    args = argparser()

    if args.print_branches:
        print_branches(args.input)
        exit(0)
    
    # MANUAL INPUT!!!!! Define the number of events you want to save and the range.
    num_events_to_read = args.read_events  # -1 means all events
    num_events_to_skip = args.skip_events  # Number of events to skip, the two parameters create a range.
    num_sets_to_save = 50000
    DATASET_LABEL = 0

    loop_storage = LoopStorage_Ttbar(args_inputs = args.inputs,
                                        args_samples_size = args.samples_size,
                                            n_events_to_read = num_events_to_read,
                                            n_events_to_skip = num_events_to_skip)
    


    # Create the H5 files
    create_H5(args.output, num_sets_to_save, loop_storage.get_total_set_size())

    # Stream the events from the ROOT file to the H5 file
    # Read one input file at a time
    for file_path, file_properties in loop_storage.file_dict.items():
        print('Info:: Processing file: %s' % file_path)
        root_file = uproot.open(file_path, num_workers=8)
        tree = root_file["Ttbar"]
        memory_size = "16 MB"
        number_of_iterations = ceil(tree.num_entries / tree.num_entries_for(memory_size))
        print('Info:: Total number of entries in the tree: %d' % tree.num_entries)
        print('Info:: Number of entries that fit in memory size of %s: %d' % (memory_size, tree.num_entries_for(memory_size)))
        print('Info:: Number of iterations to process the file in chunks of %s: %d' % (memory_size, number_of_iterations))

        # Define sets to save per-iteration
        remainder = num_sets_to_save % number_of_iterations
        sets = [ int(num_sets_to_save / number_of_iterations) ] * number_of_iterations
        if remainder > 0:
            sets[-1] += remainder
        
        sum_sets = sum(sets)    
        if sum_sets != num_sets_to_save:
            print('Error:: Mismatch in number of sets to save. Expected: %d, Got: %d' % (num_sets_to_save, sum_sets))
            exit(1)

        counter = 0
        for itree in uproot.iterate(tree, get_branches(), step_size=memory_size, library="np", num_workers=8):
            MakeDataset_Ttbar(args, itree, loop_storage, file_key = file_path, sets_to_save_list = sets, iteration=counter)
            counter += 1