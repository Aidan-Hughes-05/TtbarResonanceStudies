import argparse
from math import ceil
import os
import h5py
import numpy as np
from tqdm import tqdm
import uproot

from Common import LoopStorage, print_branches

# Add event-level branches
event_level = ["lumiWeight", "mttbar", "Coupling"]

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
class LoopStorage_Ttbar_Events(LoopStorage):
    def __init__(self, args_inputs):
        super().__init__()
        self.file_names = args_inputs.split(',')

        # Dictionary to hold per-file properties
        for file_name in self.file_names:
            self._add_file_properties(file_name)

    def _add_file_properties(self, file_name):

        self.file_dict[file_name] = {}


        class_label = 0 if ('Signal' in file_name or 'BG' in file_name) else 1
        print('Info:: Assigned class label %d to file %s' % (class_label, file_name))
        self.file_dict[file_name]['class_label'] = class_label

        root_file = uproot.open(file_name, num_workers=8)
        n_events_in_tree = root_file["Ttbar"].num_entries
        tree = root_file["Ttbar"].arrays(["lumiWeight"], library="np")
        print('Info:: Number of events in the tree: %d' % n_events_in_tree)
        self.file_dict[file_name]['n_events_in_tree'] = n_events_in_tree
        root_file.close()
                
        print('Info:: User configuration is valid for %s file. Proceeding with the dataset creation.' % file_name)

    def get_total_events_size(self):
        return sum(file_info['n_events_in_tree'] for file_info in self.file_dict.values())

def MakeDataset_Ttbar_Events(args, tree, loop_storage : LoopStorage_Ttbar_Events = None , file_key : str = None, sets_to_save_list = None, iteration : int = -1):
    # Check set to save is valid
    if sets_to_save_list == None:
        print('Error:: Sets to save list not provided.')
        exit(1)

    events_to_save = sets_to_save_list[iteration]
    if events_to_save <= 0:
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
        file_properties_dict = loop_storage.file_dict[file_key]

        # ---------------------------------------------------
        #                 Load data from ROOT
        # ---------------------------------------------------
        # Tops
        event_weight = lumiWeight
        mttbar = tree['mttbar']
        coupling = tree['Coupling']
        top_1_kinematics = get_particle_kinematics('Top_1')
        top_2_kinematics = get_particle_kinematics('Top_2')
        top_1_pt, top_1_eta, top_1_phi = tree[top_1_kinematics[0]], tree[top_1_kinematics[1]], tree[top_1_kinematics[2]]
        top_2_pt, top_2_eta, top_2_phi = tree[top_2_kinematics[0]], tree[top_2_kinematics[1]], tree[top_2_kinematics[2]]
        print('Info:: Loaded branches from ROOT file chunk.')


        # ---------------------------------------------------
        #        Loading the already existent data from h5
        # ---------------------------------------------------
        metadata_dataset = h5_file['METADATA']['H5_DATA']
        metadata_data = metadata_dataset[:]
        print(metadata_data)
        already_existent_events = int(metadata_data[0][0])
        print('Info:: Already existent events in H5 file: %d' % already_existent_events)


        range_to_change = [already_existent_events, already_existent_events + events_to_save]
        print('Debug:: Range to change in H5 file: ', range_to_change)
        particles_dataset = h5_file['INPUTS']['EVENT_PARTICLES']
        global_dataset = h5_file['INPUTS']['EVENT_GLOBAL']
        labels_dataset = h5_file['LABELS']['CLASS']
        

        # Load the arrays
        DATASET_LABEL = file_properties_dict['class_label']
        particle_data = particles_dataset[range_to_change[0]:range_to_change[1]]
        global_data = global_dataset[range_to_change[0]:range_to_change[1]]
        label_data = labels_dataset[range_to_change[0]:range_to_change[1]]
        

        print('Info:: Saving %d events...' % events_to_save)
        save_index = 0
        print(len(top_1_pt))
        for i in tqdm(range(events_to_save), desc='Streaming to HDF5 file', unit='sets'):

            particle_data[save_index] = (
                float(top_1_pt[i]),
                float(top_1_eta[i]),
                float(top_1_phi[i]),
                float(top_2_pt[i]),
                float(top_2_eta[i]),
                float(top_2_phi[i]),
            )

            global_data[save_index] = (
                float(event_weight[i]),
                float(mttbar[i]),
                float(coupling[i])
            )

            save_index += 1

        # Fill the label
        label_data = [DATASET_LABEL] * events_to_save

        # Fill the metadata
        metadata_data = (already_existent_events + events_to_save, )
            
        # ---------------------------------------------------
        #               Save to the H5 file
        # ---------------------------------------------------
      

        # Append the new data
        
        particles_dataset[range_to_change[0]:range_to_change[1]] = particle_data
        global_dataset[range_to_change[0]:range_to_change[1]] = global_data
        labels_dataset[range_to_change[0]:range_to_change[1]] = label_data
        metadata_dataset[:] = metadata_data


def create_H5(name, num_events_to_save):
    print('Info:: Creating H5 file: %s' % name)
    with h5py.File(name, 'w') as f:
        # Declare datatypes
        particles_dt = np.dtype([('top1_pt', np.float32), ('top1_eta', np.float32), ('top1_phi', np.float32), ('top2_pt', np.float32), ('top2_eta', np.float32), ('top2_phi', np.float32)])
        global_dt = np.dtype([('lumiWeight', np.float32),('mttbar', np.float32), ('coupling', np.float32)])
        metadata_dt = np.dtype([('num_saved_events', np.float32)])

        # Create groups and datasets as needed
        inputs_group = f.create_group('INPUTS')
        inputs_group.create_dataset('EVENT_PARTICLES',shape=(num_events_to_save), maxshape=(num_events_to_save), dtype=particles_dt)
        inputs_group.create_dataset('EVENT_GLOBAL',shape=(num_events_to_save), maxshape=(num_events_to_save), dtype=global_dt)
        labels_group = f.create_group('LABELS')
        labels_group.create_dataset('CLASS',shape=(num_events_to_save), maxshape=(num_events_to_save), dtype=np.float32) # When using strings ... h5py.string_dtype()
        metadata_group = f.create_group('METADATA')
        metadata_group.create_dataset('H5_DATA',shape=(1), maxshape=(1), dtype=metadata_dt)
        

    print('Info:: H5 file %s created successfully.' % name)

def argparser():

    parser = argparse.ArgumentParser(description='Convert ROOT dataset to H5 format')
    parser.add_argument('-i', '--inputs',  type=str, required=True, help='Comma separated list of ROOT files from which to read events.')
    parser.add_argument('-o', '--output', type=str, help='Output .h5 file')
    parser.add_argument('-p', '--print-branches', action='store_true', help='Print branches in the input ROOT file.', required=False)

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = argparser()

    if args.print_branches:
        print_branches(args.input)
        exit(0)
    
    # MANUAL INPUT!!!!! Define the number of events you want to save and the range.
    loop_storage = LoopStorage_Ttbar_Events(args_inputs = args.inputs)
    
    total_num_events_to_save = loop_storage.get_total_events_size()


    # Create the H5 files
    create_H5(args.output, total_num_events_to_save)

    # Stream the events from the ROOT file to the H5 file
    # Read one input file at a time
    for file_path, file_properties in loop_storage.file_dict.items():
        print('Info:: Processing file: %s' % file_path)
        root_file = uproot.open(file_path, num_workers=8)
        tree = root_file["Ttbar"]
        memory_size = "16 MB"
        num_events_to_save = tree.num_entries
        events_per_iteration = tree.num_entries_for(memory_size, get_branches())
        number_of_iterations = ceil(num_events_to_save / events_per_iteration)
        print('Info:: Total number of entries in the tree: %d' % tree.num_entries)
        print('Info:: Number of entries that fit in memory size of %s: %d' % (memory_size, events_per_iteration))
        print('Info:: Number of iterations to process the file in chunks of %s: %d' % (memory_size, number_of_iterations))

        # Define sets to save per-iteration
        sets = []
        while sum(sets) < num_events_to_save:
            sets.append(events_per_iteration)
            # Adjust last set to match exactly the number of events to save
            if sum(sets) > num_events_to_save:
                sets[-1] = num_events_to_save - sum(sets[:-1])
        
        
        sum_sets = sum(sets)    
        if sum_sets != num_events_to_save:
            print('Error:: Mismatch in number of sets to save. Expected: %d, Got: %d' % (num_events_to_save, sum_sets))
            exit(1)

        print(sets)

        counter = 0
        for itree in uproot.iterate(tree, get_branches(), step_size=memory_size, library="np", num_workers=8):
            MakeDataset_Ttbar_Events(args, itree, loop_storage, file_key = file_path, sets_to_save_list = sets, iteration=counter)
            counter += 1