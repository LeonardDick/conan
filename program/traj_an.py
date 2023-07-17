# The program is written by Leonard Dick, 2023

# MODULES
import time
import os
import pandas as pd
import numpy as np
import math
from prettytable import PrettyTable
import sys
# ----- Own modules ----- #
import defdict as ddict
import post_proc as post
import traj_info

# ARGUMENTS
args = ddict.read_commandline()

  
# MAIN
def analysis_opt(id_frame, CNT_centers, box_size, tuberadii, min_z_pore, max_z_pore, length_pore) -> None:

    # General analysis options (Is the whole trajectory necessary or just the first frame?).
    ddict.printLog('(1) Produce xyz files of the simulation box or pore structure.')
    ddict.printLog('(2) Analyze the trajectory.')
    choice = ddict.get_input('Picture or analysis mode?: ', args)
    ddict.printLog('')
    if choice == '1':
        ddict.printLog('PICTURE mode.\n', color = 'red')
        generating_pictures(id_frame, CNT_centers, box_size)
    elif choice == '2':
        ddict.printLog('ANALYSIS mode.\n', color = 'red')
        trajectory_analysis(id_frame, CNT_centers, box_size, tuberadii, min_z_pore, max_z_pore, length_pore)
    else:
        ddict.printLog('-> The choice is not known.')
    ddict.printLog('')

# Generating pictures.
def generating_pictures(id_frame, CNT_centers, box_size) -> None:
    ddict.printLog('(1) Produce xyz file of the whole simulation box.') 
    ddict.printLog('(2) Produce xyz file of empty pore structure.')
    ddict.printLog('(3) Produce xyz file of the pore structures\' tube.')
    analysis1_choice = ddict.get_input('What do you want to do?: ', args)

    if analysis1_choice == '1':
        ddict.printLog('\n-> Pics of box.')
        # Write the xyz file. The first line has the number of atoms (column in the first_drame), the second line is empty.
        id_frame_print = open('simbox_frame.xyz','w')
        id_frame_print.write('%d\n#Made with CONAN\n' % len(id_frame))
        for index, row in id_frame.iterrows():
            id_frame_print.write('%s\t%0.3f\t%0.3f\t%0.3f\n' % (row['Element'], row['x'], row['y'], row['z']))
        id_frame_print.close()
        ddict.printLog('-> Saved simulation box as simbox_frame.xyz.')

    elif analysis1_choice == '2':
        ddict.printLog('\n-> Pics of pore structure(s).')
        # Loop over the number of entries in the tuberadius numpy array.
        for i in range(len(CNT_centers)):
            # Create a dataframe with the just atoms of the respective pore structure. Extract the atoms from the id_frame. They are labeled Pore1, Pore2... in the Struc column.
            CNT_atoms_pic = id_frame.loc[id_frame['Struc'] == 'Pore%d' % (i + 1)]
            # Remove all columns except the Element, x, y, and z columns.
            CNT_atoms_pic = CNT_atoms_pic.drop(['Charge', 'Struc', 'Label', 'CNT', 'Molecule'], axis = 1)
            add_centerpoint = ddict.get_input('Add the center point of the CNT to the file? [y/n] ', args)
            if add_centerpoint == 'y':
                # Add the center point of the CNT to the dataframe, labeled as X in a new row.
                CNT_atoms_pic.loc[len(CNT_atoms_pic.index)] = ['X', CNT_centers[0][0], CNT_centers[0][1], CNT_centers[0][2]]      
            CNT_atoms_print = open(f"pore{i+1}.xyz",'w')
            CNT_atoms_print.write('%d\n#Made with CONAN\n' % len(CNT_atoms_pic))
            # Print the CNT_atoms dataframe to a xyz file. Just the atoms, x, y, and z column.
            for index, row in CNT_atoms_pic.iterrows():
                CNT_atoms_print.write('%s\t%0.3f\t%0.3f\t%0.3f\n' % (row['Element'], row['x'], row['y'], row['z']))
            CNT_atoms_print.close()
            ddict.printLog(f"-> saved as pore{i + 1}.xyz",'w')

    elif analysis1_choice == '3':
        ddict.printLog('\n-> Tube pictures.')

        for i in range(len(CNT_centers)):
            # Create a dataframe with the just atoms of the respective pore structure. Extract the atoms from the id_frame. 
            CNT_atoms_pic = pd.DataFrame(id_frame.loc[id_frame['CNT'] == i + 1])
            # Remove all columns except the Element, x, y, and z columns.
            CNT_atoms_pic = CNT_atoms_pic.drop(['Charge', 'Struc', 'Label', 'CNT', 'Molecule'], axis=1)
            add_liquid = ddict.get_input(f"Add liquid which is inside the CNT{i + 1}? [y/n] ", args)

            if add_liquid == 'y':
                add_liquid2 = ddict.get_input('Add all confined atoms (1), or entire molecules (2) ? [1/2] ', args)

                if add_liquid2 == '1':
                    # Scan the id_frame and add all atoms which are inside the tube to the tube_atoms dataframe.
                    for index, row in id_frame.iterrows():                                          
                        if row['Struc'] == 'Liquid' and row['z'] <= CNT_atoms_pic['z'].max() and row['z'] >= CNT_atoms_pic['z'].min():
                            # Add the row to the tube_atoms dataframe.
                            CNT_atoms_pic.loc[index] = [row['Element'], row['x'], row['y'], row['z']]

                elif add_liquid2 == '2':
                    # Do the molecule recognition.
                    ddict.printLog('-> Molecule recognition.')
                    id_frame = traj_info.molecule_recognition(id_frame, box_size)
                    id_frame = id_frame.drop(['Charge', 'Label', 'CNT'], axis = 1)
                    # Add the Molecule column to the CNT_atoms_pic dataframe.
                    CNT_atoms_pic['Molecule'] = np.nan
                    # Scan the id_frame and add all atoms which are inside the tube to the tube_atoms dataframe.
                    for index, row in id_frame.iterrows():                                          
                        if row['Struc'] == 'Liquid' and row['z'] <= CNT_atoms_pic['z'].max() and row['z'] >= CNT_atoms_pic['z'].min():
                            # Add the row to the tube_atoms dataframe.
                            CNT_atoms_pic.loc[index] = [row['Element'], row['x'], row['y'], row['z'], row['Molecule']]

                    # List the molecules which are inside the tube.
                    mol_list = []
                    mol_list.append(CNT_atoms_pic['Molecule'].unique())
                    tube_atoms_mol = pd.DataFrame(columns = ['Element', 'x', 'y', 'z', 'Molecule'])
                    mol_list = mol_list[0]
                    # Scan the id_frame and add all atoms which are in the mol_list to the tube_atoms_mol dataframe.
                    for index, row in id_frame.iterrows():                                          
                        if row['Molecule'] in mol_list:
                            # Add the row to the tube_atoms dataframe.
                            tube_atoms_mol.loc[index] = [row['Element'], row['x'], row['y'], row['z'], row['Molecule']]
                    # Append the tube_atoms_mol dataframe to the tube_atoms_pic dataframe.
                    CNT_atoms_pic = pd.concat([CNT_atoms_pic, tube_atoms_mol], ignore_index = True)

                    # Finally remove all duplicates from the tube_atoms_pic dataframe.
                    CNT_atoms_pic = CNT_atoms_pic.drop_duplicates(subset = ['Element', 'x', 'y', 'z', 'Molecule'], keep = 'first')

            else:
                add_centerpoint = ddict.get_input(f"Add the center point of the CNT{i+1} to the file? [y/n] ")
                if add_centerpoint == 'y':         
                    CNT_atoms_pic.loc[len(CNT_atoms_pic.index)] = ['X', CNT_centers[0][0], CNT_centers[0][1], CNT_centers[0][2]]      

            tube_atoms_print = open(f"CNT{i + 1}.xyz",'w')
            tube_atoms_print.write('%d\n#Made with CONAN\n' % len(CNT_atoms_pic))

            for index, row in CNT_atoms_pic.iterrows():
                tube_atoms_print.write('%s\t%0.3f\t%0.3f\t%0.3f\n' % (row['Element'], row['x'], row['y'], row['z']))
            tube_atoms_print.close()
            ddict.printLog(f"-> Saved as CNT{i + 1}.xyz")
    else:
        ddict.printLog('\nThe analysis you entered is not known.')
    ddict.printLog('')

# Analysis of the trajectory.
def trajectory_analysis(id_frame, CNT_centers, box_size, tuberadii, min_z_pore, max_z_pore, length_pore) -> None:

    # Analysis choice.
    ddict.printLog('(1) Calculate the radial density inside the CNT')
    ddict.printLog('(2) Calculate the radial charge density inside the CNT (just for .pdb format)')
    ddict.printLog('(3) Calculate the accessibe volume of the CNT')
    ddict.printLog('(4) Calculate the average density along the z axis of the simulation box')
    analysis_choice2 = ddict.get_input('Which analysis should be conducted?:  ', args)

    if analysis_choice2 == '1' or analysis_choice2 == '2':
        for i in range(len(CNT_centers)):
            ddict.printLog(f'\n-> CNT{i + 1}')
            num_increments = int(ddict.get_input('How many increments do you want to use to calculate the density profile? ', args))
            rad_increment = tuberadii[i] / num_increments
            # Make an array which start at 0 and end at tuberadius with num_increments + 1 steps.
            raddens_bin_edges = np.linspace(0, tuberadii[0], num_increments + 1)
            # Define raddens_bin_labels, they are a counter for the bin edges.
            raddens_bin_labels = np.arange(1, len(raddens_bin_edges), 1)
            ddict.printLog('Increment distance: %0.3f angstrom' % (rad_increment))
            
    elif analysis_choice2 == '3':
        ddict.printLog('')

    elif analysis_choice2 == '4':
        ddict.printLog('')
        num_increments = int(ddict.get_input('How many increments per section do you want to use to calculate the density profile? ', args))
        ddict.printLog('\nNote here, that the simulation box is subdivided between the bulk phases and the tube. \nThe number of increments set here is the number of increments for each section.\n', color='red')
        # Initialize arrays
        z_incr_CNT = [0] * len(CNT_centers)
        z_incr_bulk1 = [0] * len(CNT_centers)
        z_incr_bulk2 = [0] * len(CNT_centers)
        z_bin_edges_pore = [0] * len(CNT_centers)
        z_bin_edges_bulk1 = [0] * len(CNT_centers)
        z_bin_edges_bulk2 = [0] * len(CNT_centers)
        z_bin_edges = [0] * len(CNT_centers)
        z_bin_labels = [0] * len(CNT_centers)

        # Calculate the increment distance for each section.
        for i in range(len(CNT_centers)):
            z_incr_CNT[i] = length_pore[i] / num_increments
            z_incr_bulk1[i] = min_z_pore[i] / num_increments
            z_incr_bulk2[i] = (box_size[2] - max_z_pore[i]) / num_increments

            ddict.printLog('Increment distance CNT: %0.3f Ang' % (z_incr_CNT[i]))
            ddict.printLog('Increment distance bulk1: %0.3f Ang' % (z_incr_bulk1[i]))
            ddict.printLog('Increment distance bulk2: %0.3f Ang' % (z_incr_bulk2[i]))

            tube_atoms = id_frame.loc[id_frame['CNT'] == i + 1]
            z_bin_edges_pore[i] = np.linspace(tube_atoms['z'].min(), tube_atoms['z'].max(), num_increments+1)
            z_bin_edges_bulk1[i] = np.linspace(0, min_z_pore[i], num_increments+1)
            z_bin_edges_bulk2[i] = np.linspace(max_z_pore[i], box_size[2], num_increments+1)
            z_bin_edges[i] = np.concatenate((z_bin_edges_bulk1[i], z_bin_edges_pore[i], z_bin_edges_bulk2[i]))
            z_bin_edges[i] = np.unique(z_bin_edges[i])
            num_increments = len(z_bin_edges[i])-1

            ddict.printLog('\nTotal number of increments: %d' % (num_increments))
            ddict.printLog('Number of edges: %s' % len((z_bin_edges[i])))
            z_bin_labels[i] = np.arange(1, len(z_bin_edges[i]), 1)
    else:
        ddict.printLog('-> The analysis you entered is not known.')
        sys.exit(1)
    ddict.printLog('')


    # MOLECULAR RECOGNITION
    # Perform the molecule recognition by loading the module molidentifier.
    id_frame = traj_info.molecule_recognition(id_frame, box_size)
    species_max = id_frame['Species'].max()
    spec_molecule = 0
    spec_atom = []
    ddict.printLog('')
    analysis_spec_molecule = ddict.get_input(f"Do you want to perform the analysis for a specific molecule kind? (y/n) ", args)
    if analysis_spec_molecule == 'y':
        spec_molecule = int(ddict.get_input(f"Which species to analyze? (1-{species_max}) ", args))
        # Ask user for the atom type to analyze. Multiple options are possible, default is 'all'.
        spec_atom = ddict.get_input(f"Which atoms to analyze? [default:all] ", args)

        if spec_atom == "" or spec_atom == "[default:all]":
            spec_atom = "all"

        # Get the atom type into a list.
        spec_atom = spec_atom.replace(', ', ',').split(',')
        ddict.printLog(f"\n-> Species {spec_molecule} and atom type {spec_atom} will be analyzed.\n")


    # GENERAL INFORMATION ON CHUNKS
    ddict.printLog('-> Reading the trajectory.\n')
    trajectory_file_size = os.path.getsize(args["trajectoryfile"])
    # Calculate how many atoms each frame has.
    number_of_atoms = len(id_frame)
    # Calculate how many lines the trajectory file has.
    with open(args["trajectoryfile"]) as f:                                    
        number_of_lines = sum(1 for i in f)
    
    lines_per_frame = 0
    # Calculate how many frames the trajectory file has.
    if args["trajectoryfile"].endswith('.xyz') or args["trajectoryfile"].endswith('.pdb'):
        lines_per_frame = number_of_atoms + 2
    elif args["trajectoryfile"].endswith('.lammpstrj') or args["trajectoryfile"].endswith('.lmp'):
        lines_per_frame = number_of_atoms + 9
    
    number_of_frames = int(number_of_lines/lines_per_frame)

    # Calculate how many bytes each line of the trajectory file has.                       
    bytes_per_line = trajectory_file_size/(number_of_lines) 
    # The number of lines in a chunk. Each chunk is roughly 50 MB large.                          
    chunk_size = int(100000000 / ((lines_per_frame) * bytes_per_line))
    # The number of chunks (always round up).         
    number_of_chunks = math.ceil(number_of_frames/chunk_size)
    # The number of frames in the last chunk.                  
    last_chunk_size = number_of_frames - (number_of_chunks - 1) * chunk_size      
    number_of_bytes_per_chunk = chunk_size * (lines_per_frame) * bytes_per_line
    number_of_lines_per_chunk = chunk_size * (lines_per_frame)
    number_of_lines_last_chunk = last_chunk_size * (lines_per_frame)
    # Table with the information on the trajectory file.
    table = PrettyTable(['', 'Trajectory', 'Chunk(%d)' % (number_of_chunks)])   
    table.add_row(['Size in MB', '%0.1f' % (trajectory_file_size / 1000000), '%0.1f (%0.1f)' % (number_of_bytes_per_chunk / 1000000, last_chunk_size * (lines_per_frame) * bytes_per_line / 1000000)])
    table.add_row(['Frames', number_of_frames, '%d(%d)' % (chunk_size, last_chunk_size)])
    table.add_row(['Lines', number_of_lines, number_of_lines_per_chunk])
    ddict.printLog(table)
    ddict.printLog('')



    # PREPERATION
    # Main loop preperation.
    counter=0
    Main_time=time.time()

    # Get atomic masses.
    element_masses = ddict.dict_mass()

    # First it is necessary to get the molecule number of each atom to attach it to every frame later in the loop.
    molecule_id = id_frame['Molecule'].values  
    molecule_species = id_frame['Species'].values 
    molecule_structure = id_frame['Struc'].values
    molecule_label = id_frame['Label'].values    

    # Analysis preperation.
    # Radial density.
    if analysis_choice2 == '1' or analysis_choice2 == '2':
        raddens_bin_labels = np.arange(0, num_increments, 1)	
        # Make new dataframe with the number of frames of the trajectory.
        raddens_df_dummy = pd.DataFrame(np.arange(1,number_of_frames + 1),columns = ['Frame'])            
        # Add a column to the dataframe for each increment.
        for i in range(num_increments):        
            raddens_df_dummy['Bin %d' % (i+1)] = 0
            raddens_df_dummy = raddens_df_dummy.copy()
        raddens_df = raddens_df_dummy.copy()

    # Accessible volume.
    if analysis_choice2 == '3' or analysis_choice2 == '4':
        maxdisp_atom_dist = 0
        which_element_radii = ddict.get_input('Do you want to use the van der Waals radii (1) or the covalent radii (2) of the elements? [1/2] ', args)
        if which_element_radii == '1':
            ddict.printLog('-> Using van der Waals radii.')
            element_radii = ddict.dict_vdW()
        elif which_element_radii == '2':
            ddict.printLog('-> Using covalent radii.')
            element_radii = ddict.dict_covalent()

    # Axial density.
    if analysis_choice2 == '4':
        # Make new dataframe with the number of frames of the trajectory.
        zdens_df_dummy = pd.DataFrame(np.arange(1,number_of_frames + 1),columns = ['Frame'])            
        for i in range(num_increments):      
            zdens_df_dummy['Bin %d' % (i + 1)] = 0
            zdens_df_dummy = zdens_df_dummy.copy()
        zdens_df = zdens_df_dummy.copy()

    # MAIN LOOP
    # Define which function to use reading the trajectory file. Pull definition from traj_info.py.
    if args["trajectoryfile"].endswith(".xyz"):
        from traj_info import xyz as run
    elif args["trajectoryfile"].endswith(".pdb"):
        from traj_info import pdb as run
    elif args["trajectoryfile"].endswith(".lmp"):
        from traj_info import lammpstrj as run
        
    # The trajectory xyz file is read in chunks of size chunk_size. The last chunk is smaller than the other chunks.
    trajectory = pd.read_csv(args["trajectoryfile"], chunksize = number_of_lines_per_chunk, header = None)        
    chunk_number = 0
    # Loop over chunks.
    for chunk in trajectory:    
        chunk_number = chunk_number + 1
        print('')
        print('Chunk %d of %d' % (chunk_number, number_of_chunks))
        # Divide the chunk into individual frames. If the chunk is the last chunk, the number of frames is different.
        if chunk.shape[0] == number_of_lines_last_chunk:                            
            frames = np.split(chunk, last_chunk_size)
        else:
            frames = np.split(chunk, chunk_size)
        for frame in frames:

            # First load the frame into the function run() to get a dataframe. Then reset the index.
            frame, split_frame = run(frame, element_masses, molecule_id)
            split_frame.reset_index(drop = True, inplace = True)
            
            # Add the necessary columns to the dataframe.
            split_frame['Struc'] = molecule_structure
            split_frame['Molecule'] = molecule_id
            split_frame['Species'] = molecule_species
            split_frame['Label'] = molecule_label  

            # Drop all CNT and carbon_wall atoms, just the Liquid atoms are needed for the analysis. 
            split_frame = split_frame[split_frame['Struc'] == 'Liquid']
            split_frame = split_frame.drop(['Struc'], axis = 1)
      
            # Drop the other atoms which are not needed for the analysis.
            if analysis_spec_molecule == 'y':
                split_frame = split_frame[split_frame['Species'].astype(int) == int(spec_molecule)]  
                #if the spec_atom list does not contain "all" then only the atoms in the list are kept.
                if spec_atom[0] != 'all':
                    # If specific atoms are requested, only these atoms are kept.
                    split_frame=split_frame[split_frame['Label'].isin(spec_atom)]
                
            # For analysis 1 and 2: All atoms are dropped which have a z coordinate larger than the maximum z coordinate of the CNT or smaller than the minimum z coordinate of the CNT.
            if analysis_choice2 == '1' or analysis_choice2 == '2' or analysis_choice2 == '3':
                split_frame = split_frame[split_frame['Z'].astype(float) <= max_z_pore[0]]   
                split_frame = split_frame[split_frame['Z'].astype(float) >= min_z_pore[0]]


            # RADIAL DENSITY FUNCTION                
            if analysis_choice2 == '1':
                # Calculate the radial density function with the remaining atoms.
                split_frame['X_adjust'] = split_frame['X'].astype(float) - CNT_centers[0][0]
                split_frame['Y_adjust'] = split_frame['Y'].astype(float) - CNT_centers[0][1]
                split_frame['Distance'] = np.sqrt(split_frame['X_adjust'] ** 2 + split_frame['Y_adjust'] ** 2)
                split_frame['Distance_bin'] = pd.cut(split_frame['Distance'], bins = raddens_bin_edges, labels = raddens_bin_labels+1)

                # Add all masses of the atoms in each bin to the corresponding bin.
                raddens_df_temp = split_frame.groupby(pd.cut(split_frame['Distance'], raddens_bin_edges))['Mass'].sum().reset_index(name = 'Weighted_counts')
                raddens_df_temp = pd.DataFrame(raddens_df_temp)

                # Add a new first column with the index+1 of the bin.
                raddens_df_temp.insert(0, 'Bin', raddens_df_temp.index+1)          

                # Write the results into the raddens_df dataframe. The row is defined by the frame number.
                for i in range(num_increments):
                    raddens_df.loc[counter,'Bin %d' % (i + 1)] = raddens_df_temp.loc[i,'Weighted_counts']
                    
                # Remove the raddens_df_temp dataframe every loop.
                del raddens_df_temp    


            # RAADIAL CHARGE DENSITY FUNCTION
            if analysis_choice2 == '2':

                # Calculate the radial density function with the remaining atoms.
                split_frame['X_adjust'] = split_frame['X'].astype(float) - CNT_centers[0][0]
                split_frame['Y_adjust'] = split_frame['Y'].astype(float) - CNT_centers[0][1]
                split_frame['Distance'] = np.sqrt(split_frame['X_adjust'] ** 2 + split_frame['Y_adjust'] ** 2)
                split_frame['Distance_bin'] = pd.cut(split_frame['Distance'], bins = raddens_bin_edges, labels = raddens_bin_labels+1)

                # Add all masses of the atoms in each bin to the corresponding bin.
                raddens_df_temp = split_frame.groupby(pd.cut(split_frame['Distance'], raddens_bin_edges))['Charge'].sum().reset_index(name='Weighted_counts')
                raddens_df_temp = pd.DataFrame(raddens_df_temp)

                # Add a new first column with the index+1 of the bin.
                raddens_df_temp.insert(0, 'Bin', raddens_df_temp.index+1)      

                # Write the results into the raddens_df dataframe. The row is defined by the frame number.
                for i in range(num_increments):
                    raddens_df.loc[counter,'Bin %d' % (i + 1)] = raddens_df_temp.loc[i,'Weighted_counts']

                # Remove the raddens_df_temp dataframe every loop.
                del raddens_df_temp    


            # Z-DENSITY PROFILE
            if analysis_choice2 == '4':

                # For now we concentate the edges/bins (This has to be changed in the future for the analysis on multiple CNTs).
                z_bin_edges = np.ravel(z_bin_edges)
                z_bin_labels = np.ravel(z_bin_labels)

                # Make an array from split_frame['Z'].
                split_frame['Z_bin'] = pd.cut(split_frame['Z'].astype(float).values, bins=z_bin_edges, labels=z_bin_labels)

                # Add all masses of the atoms in each bin to the corresponding bin. Then add a new first column with the index+1 of the bin.
                zdens_df_temp = split_frame.groupby(pd.cut(split_frame['Z'].astype(float), z_bin_edges))['Mass'].sum().reset_index(name='Weighted_counts')
                zdens_df_temp = pd.DataFrame(zdens_df_temp)
                zdens_df_temp.insert(0, 'Bin', zdens_df_temp.index+1)           

                # Write the results into the zdens_df dataframe. The row is defined by the frame number.
                for i in range(num_increments):
                    zdens_df.loc[counter,'Bin %d' % (i + 1)]=zdens_df_temp.loc[i,'Weighted_counts']

                # TUBE SECTION -> displacement for accessible volume.
                # All atoms are dropped which have a z coordinate larger than the maximum z coordinate of the CNT or smaller than the minimum z coordinate of the CNT.
                split_frame = split_frame[split_frame['Z'].astype(float) <= max_z_pore[0]]   
                split_frame = split_frame[split_frame['Z'].astype(float) >= min_z_pore[0]] 


            # Z-DENSITY PROFILE AND ACCESSIBLE VOLUME             
            if analysis_choice2 == '3' or analysis_choice2 == '4':    
                # Identify the atom inside the tube, which is the most displaced from its center. Over the whole trajectory.
                split_frame['X_adjust'] = split_frame['X'].astype(float) - CNT_centers[0][0]
                split_frame['Y_adjust'] = split_frame['Y'].astype(float) - CNT_centers[0][1]
                # While caluclating the distance, add the elements radius.
                split_frame['Distance'] = np.sqrt(split_frame['X_adjust'] ** 2 + split_frame['Y_adjust'] ** 2) + split_frame['Atom'].map(element_radii)
                if split_frame['Distance'].max() > maxdisp_atom_dist:
                    maxdisp_atom_dist = split_frame['Distance'].max()
                    # Save the complete info of the most displaced atom.
                    maxdisp_atom_row = split_frame.loc[split_frame['Distance'].idxmax()]

            counter += 1
            print('Frame %d of %d' % (counter, number_of_frames), end = '\r')
    ddict.printLog('')
    ddict.printLog('Finished processing the trajectory. %d frames were processed.' % (counter))
    ddict.printLog('')
    
    # DATA PROCESSING
    if analysis_choice2 == '1' or analysis_choice2 == '2':
        post.raddens_func(raddens_df, raddens_bin_edges, number_of_frames, length_pore[0], tuberadii[0], analysis_choice2, args)
    if analysis_choice2 == '3':
        post.acc_vol_func(id_frame, maxdisp_atom_row, length_pore[0], args)
    if analysis_choice2 == '4':
        post.zdens_func(zdens_df, z_bin_edges, number_of_frames, maxdisp_atom_row, length_pore[0], tuberadii[0], num_increments, min_z_pore[0], max_z_pore[0], box_size[0], box_size[1], CNT_centers[0][2], args)


    ddict.printLog("The main loop took %0.3f seconds to run" % (time.time() - Main_time))
