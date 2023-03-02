#!/usr/bin/python3.8
# The program is written by Leonard Dick, 2022

# MODULES
import time
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
import math

import matplotlib.pyplot as plt
from prettytable import PrettyTable
from scipy.spatial.distance import cdist
import sys
#----------#
import defdict as ddict
import post_proc as post

#ARGUMENTS
args = ddict.read_commandline()

#DICTIONARIES 
# Atomic masses and van der Waals radii
element_masses = ddict.dict_mass()
  
#MAIN
def analysis_opt(id_frame, simbox_x, simbox_y, simbox_z, number_of_atoms, tube_atoms, tuberadius, tubecenter_x, tubecenter_y, tubecenter_z, carbon_wall_z, length_CNT, max_z_CNT, min_z_CNT, CNTs, carbon_walls, CNT_atoms, pdb_info) -> None:
    #General analysis options (Is the whole trajectory necessary or just one frame?)
    ddict.printLog('What analysis do you want to perform?')
    ddict.printLog('(1) Produce xyz files of the simulation box or pore structures in the chosen frame.')
    ddict.printLog('(2) Analyze the trajectory')
    analysis_choice=ddict.get_input('Input: ')
    ddict.printLog('')
    if analysis_choice == '1':
        ddict.printLog_red('PICTURE mode.')
    elif analysis_choice == '2':
        ddict.printLog_red('ANALYSIS mode.')
        which_element_radii = ddict.get_input('Do you want to use the van der Waals radii (1) or the covalent radii (2) of the elements? [1/2] ')
        if which_element_radii == '1':
            ddict.printLog('-> Using van der Waals radii.')
            element_radii = ddict.dict_vdW()
        elif which_element_radii == '2':
            ddict.printLog('-> Using covalent radii.')
            element_radii = ddict.dict_covalent()
    else:
        ddict.printLog('-> The choice you entered is not known.')
    ddict.printLog('')

    #Analysis 1
    if analysis_choice == '1':
        ddict.printLog('What do you want to do?')
        ddict.printLog('(1) Produce xyz file of the whole simulation box.')
        ddict.printLog('(2) Produce xyz file of empty pore structure.')
        ddict.printLog('(3) Produce xyz file of the pore structures\' CNT.')
        analysis1_choice=ddict.get_input('Input: ')
        ddict.printLog('')
        if analysis1_choice == '1':
            ddict.printLog('-> Pics of box.')
            #print the id_frame dataframe to a xyz file. Just the atoms, x, y, and z column
            id_frame_pic=id_frame.drop(['Molecule','Frozen'], axis=1)
            #Write the xyz file. The first line has the number of atoms (column in the first_drame), the second line is empty. The come all atoms with the x,y and z coordinate.
            id_frame_print=open('simbox_frame.xyz','w')
            id_frame_print.write('%d\n#Made with CONAN\n' % len(id_frame_pic))
            for index, row in id_frame_pic.iterrows():
                id_frame_print.write('%s\t%0.3f\t%0.3f\t%0.3f\n' % (row['Atom'], row['x'], row['y'], row['z']))
            id_frame_print.close()
            ddict.printLog('-> Saved as simbox_frame.xyz.')
        elif analysis1_choice == '2':
            ddict.printLog('-> Pics of pore structure.')
            CNT_atoms_pic=CNT_atoms.drop(['xy_distance'], axis=1)
            add_centerpoint=ddict.get_input('Add the center point of the CNT to the file? [y/n] ')
            if add_centerpoint == 'y':
                #Add the center point of the CNT to the file, labeled as X.
                CNT_atoms_pic.loc[-1] = ['X', tubecenter_x, tubecenter_y, tubecenter_z]            
            CNT_atoms_print=open('pore.xyz','w')
            CNT_atoms_print.write('%d\n#Made with CONAN\n' % len(CNT_atoms_pic))
            #print the CNT_atoms dataframe to a xyz file. Just the atoms, x, y, and z column
            for index, row in CNT_atoms_pic.iterrows():
                CNT_atoms_print.write('%s\t%0.3f\t%0.3f\t%0.3f\n' % (row['Atom'], row['x'], row['y'], row['z']))
            CNT_atoms_print.close()
            ddict.printLog('-> Saved as pore.xyz.')
        elif analysis1_choice == '3':
            ddict.printLog('-> Pics of CNT.')
            tube_atoms_pic=tube_atoms.drop(['xy_distance'], axis=1)
            #reset the index of the dataframe, so that the index starts at 0.
            tube_atoms_pic.reset_index(drop=True, inplace=True)
            id_frame.reset_index(drop=True, inplace=True)
            add_liquid=ddict.get_input('Add liquid which is inside the CNT? [y/n] ')
            if add_liquid == 'y':
                add_liquid2=ddict.get_input('Add all confined atoms (1), or entire molecules (2) ? [1/2] ')
                if add_liquid2 == '1':
                    #Scan the id_frame and add all atoms which are inside the tube to the tube_atoms dataframe.
                    for index, row in id_frame.iterrows():                                          
                        if row['Frozen'] == False and row['z']<=tube_atoms_pic['z'].max() and row['z']>=tube_atoms_pic['z'].min():
                            #Add the row to the tube_atoms dataframe.
                            tube_atoms_pic.loc[index] = [row['Atom'], row['x'], row['y'], row['z']]
                elif add_liquid2 == '2':
                    #call the molecule identification function
                    import molidentifier
                    #perform the molecule recognition by loading the module molidentifier.
                    species_num, id_frame=molidentifier.molident(id_frame, CNTs, carbon_walls, pdb_info)
                    id_frame.reset_index(drop=True, inplace=True)
                    print(id_frame)
                    #add molecule column to the tube_atoms dataframe.
                    tube_atoms_pic['Molecule'] = 0
                    print(tube_atoms_pic)
                    #Scan the id_frame and add all atoms which are inside the tube.
                    for index, row in id_frame.iterrows():                                          
                        if row['Frozen'] == False and row['z']<=tube_atoms_pic['z'].max() and row['z']>=tube_atoms_pic['z'].min():
                            #Add the row to the tube_atoms dataframe.
                            tube_atoms_pic.loc[index] = [row['Atom'], row['x'], row['y'], row['z'], row['Molecule']] 
                            #make a list of the molecules which are inside the tube.
                    mol_list=[]
                    mol_list.append(tube_atoms_pic['Molecule'].unique())
                    tube_atoms_mol=pd.DataFrame(columns=['Atom', 'x', 'y', 'z', 'Molecule'])
                    print(mol_list)
                    mol_list=mol_list[0]
                    print(mol_list)
                    #Scan the id_frame and add all atoms which are in the mol_list to the tube_atoms_mol dataframe.
                    for index, row in id_frame.iterrows():                                          
                        if  row['Molecule'] in mol_list:
                            #Add the row to the tube_atoms dataframe.
                            tube_atoms_mol.loc[index] = [row['Atom'], row['x'], row['y'], row['z'], row['Molecule']]
                    #append the tube_atoms_mol dataframe to the tube_atoms_pic dataframe.
                    tube_atoms_pic=tube_atoms_pic.append(tube_atoms_mol, ignore_index=True)
                    print(tube_atoms_pic)

                            
            else:
                add_centerpoint=ddict.get_input('Add the center point of the CNT to the file? [y/n] ')
                if add_centerpoint == 'y':         
                    tube_atoms_pic.loc[-1] = ['X', tubecenter_x, tubecenter_y, tubecenter_z]
            tube_atoms_print=open('CNT.xyz','w')
            tube_atoms_print.write('%d\n#Made with CONAN\n' % len(tube_atoms_pic))
            for index, row in tube_atoms_pic.iterrows():
                tube_atoms_print.write('%s\t%0.3f\t%0.3f\t%0.3f\n' % (row['Atom'], row['x'], row['y'], row['z']))
            tube_atoms_print.close()
            ddict.printLog('-> Saved as CNT.xyz.')
        else:
            ddict.printLog('The analysis you entered is not known.')
        ddict.printLog('')

    #Analysis 2
    if analysis_choice == '2':
        ddict.printLog('What do you want to do?')
        ddict.printLog('(1) Calculate the radial density inside the CNT')
        ddict.printLog('(2) Calculate the accessibe volume of the CNT')
        ddict.printLog('(3) Calculate the average density along the z axis of the simulation box')
        analysis_choice2=ddict.get_input('Input:  ')
        ddict.printLog('')
        if analysis_choice2 == '1':
            num_increments=int(ddict.get_input('How many increments do you want to use to calculate the density profile? '))
            rad_increment=tuberadius/num_increments
            #make an array which start at 0 and end at tuberadius with num_increments+1 steps.
            raddens_bin_edges=np.linspace(0, tuberadius, num_increments+1)
            #define raddens_bin_labels, they are a counter for the bin edges. The first bin is 1, the second bin is 2, etc.
            raddens_bin_labels=np.arange(1, len(raddens_bin_edges), 1)
            print('Increment distance: %0.3f angstrom' % (rad_increment))
        elif analysis_choice2 == '2':
            ddict.printLog('')
        elif analysis_choice2 == '3':
            num_increments=int(ddict.get_input('How many increments per section do you want to use to calculate the density profile? '))
            ddict.printLog('Note here, that the simulation box is subdivided between the bulk phases and the tube. \nThe number of increments set here is the number of increments for each section.')
            #define the increment distance for the CNT and the bulk phases.
            z_incr_bulk1=min_z_CNT/num_increments
            z_incr_bulk2=(simbox_z-max_z_CNT)/num_increments
            z_incr_CNT=length_CNT/num_increments
            print('Increment distance CNT: %0.3f Ang' % (z_incr_CNT))
            print('Increment distance bulk1: %0.3f Ang' % (z_incr_bulk1))
            print('Increment distance bulk2: %0.3f Ang' % (z_incr_bulk2))
            #start an array starting at the minimum z value of the CNT and ending at the maximum z value of the CNT with num_increments+1 steps.
            z_bin_edges_CNT=np.linspace(tube_atoms['z'].min(), tube_atoms['z'].max(), num_increments+1)
            #do the same for both bulks, with the right borders.
            z_bin_edges_bulk1=np.linspace(0, min_z_CNT, num_increments+1)
            z_bin_edges_bulk2=np.linspace(max_z_CNT, simbox_z, num_increments+1)
            #combine all bins, combining the bulk phases and the CNT.
            z_bin_edges=np.concatenate((z_bin_edges_bulk1, z_bin_edges_CNT, z_bin_edges_bulk2))
            #drop duplicates from the array.        
            z_bin_edges=np.unique(z_bin_edges)
            #define the number of increments, which is the number of edges minus 1.        
            num_increments=len(z_bin_edges)-1
            ddict.printLog('Total number of increments: %d' % (num_increments))
            ddict.printLog('Number of edges: %s' % len((z_bin_edges)))
            #define z_bin_labels, they are a counter for the bin edges. The first bin is  1, the second bin is 2, etc.
            z_bin_labels=np.arange(1, len(z_bin_edges), 1)
        else:
            ddict.printLog('-> The analysis you entered is not known.')
            sys.exit(1)
        ddict.printLog('')
        analysis_choice2_molrec=ddict.get_input('Do you need molecular recognition for the analysis? [y/n] ')


        #MOLECULE IDENTIFIER
        if analysis_choice2_molrec == 'y':
            import molidentifier
            #perform the molecule recognition by loading the module molidentifier.
            species_num, id_frame=molidentifier.molident(id_frame, CNTs, carbon_walls, pdb_info)
        else:
            ddict.printLog('-> No liquid molecules will be identified.')
            ddict.printLog('')

        #GENERAL INFORMATION ON CHUNKS
        trajectory_file_size = os.path.getsize(args["trajectoryfile"])
        #Calculate how many lines the trajectory file has.
        with open(args["trajectoryfile"]) as f:                                    
            for number_of_lines, line in enumerate(f):
                pass
        number_of_lines += 1
        #Calculate how many frames the trajectory file has.
        number_of_frames = int((number_of_lines)/(number_of_atoms + 2))
        #Calculate how many bytes each line of the trajectory file has.                       
        bytes_per_line = trajectory_file_size/(number_of_lines) 
        #The number of lines in a chunk. Each chunk is roughly 20 MB large.                          
        chunk_size = int(20000000/((number_of_atoms + 2)*bytes_per_line))
        #The number of chunks (always round up)           
        number_of_chunks = math.ceil(number_of_frames/chunk_size)
        #The number of frames in the last chunk                   
        last_chunk_size = number_of_frames - (number_of_chunks - 1)*chunk_size      
        number_of_bytes_per_chunk = chunk_size*(number_of_atoms + 2)*bytes_per_line
        number_of_lines_per_chunk = chunk_size*(number_of_atoms + 2)
        number_of_lines_last_chunk = last_chunk_size*(number_of_atoms + 2)
        #Table with the information on the trajectory file.
        table = PrettyTable(['', 'Trajectory', 'Chunk(%d)' % (number_of_chunks)])   
        table.add_row(['Size in MB', '%0.1f' % (trajectory_file_size/1000000), '%0.1f (%0.1f)' % (number_of_bytes_per_chunk/1000000, last_chunk_size*(number_of_atoms + 2)*bytes_per_line/1000000)])
        table.add_row(['Frames', number_of_frames, '%d(%d)' % (chunk_size, last_chunk_size)])
        table.add_row(['Lines', number_of_lines, number_of_lines_per_chunk])
        ddict.printLog(table)
        ddict.printLog('')
        #MAIN LOOP - PREPERATION MOLIDENT
        counter=0
        Main_time=time.time()
        #First it is necessary to get the molecule number of each atom to attach it to every frame later in the loop.
        molecule_id=id_frame['Molecule'].values                                  
        analysis_spec_molecule='n'
        #If there is a species column in the id_frame, cut it out and save it in a separate array.
        if 'Species' in id_frame.columns:                                        
            species_id=id_frame['Species'].values
            analysis_spec_molecule=ddict.get_input('Do you want to perform the analysis for a specific molecule kind? [y/n] ')
            if analysis_spec_molecule=='y':
                spec_molecule=int(ddict.get_input('Which species to analyze? (1-%d) ' % (species_num['Species'].max())))

        #Analysis preperation
        #RADIAL DENSITY - PREPARATION
        if analysis_choice2=='1':
            raddens_bin_labels = np.arange(0, num_increments, 1)	
            #Make new dataframe with the number of frames of the trajectory.
            raddens_df_dummy=pd.DataFrame(np.arange(1,number_of_frames+1),columns=['Frame'])            
            #Add a column to the dataframe for each increment.
            for i in range(num_increments):        
                raddens_df_dummy['Bin %d' % (i+1)]=0
                raddens_df_dummy=raddens_df_dummy.copy()
            raddens_df=raddens_df_dummy.copy()

        #ACCESSIBLE VOLUME - PREPARATION
        if analysis_choice2=='2':
            maxdisp_atom_dist=0

        #Z-DENS - PREPARATION
        if analysis_choice2=='3':
            #zdens_bin_labels = np.arange(0, num_increments, 1)
            #Make new dataframe with the number of frames of the trajectory.
            zdens_df_dummy=pd.DataFrame(np.arange(1,number_of_frames+1),columns=['Frame'])            
            for i in range(num_increments):        #Add a column for each increment.
                zdens_df_dummy['Bin %d' % (i+1)]=0
                zdens_df_dummy=zdens_df_dummy.copy()
            zdens_df=zdens_df_dummy.copy()
            maxdisp_atom_dist=0

        #MAIN LOOP
        #define which function to use reading the trajectory file. Pull definition from traj_info.py
        if args["trajectoryfile"].endswith(".xyz"):
            from traj_info import xyz as run
        elif args["trajectoryfile"].endswith(".pdb"):
            from traj_info import pdb as run
        elif args["trajectoryfile"].endswith(".lmp"):
            from traj_info import lammpstrj as run
            
        #The trajectory xyz file is read in chunks of size chunk_size. The last chunk is smaller than the other chunks.
        trajectory=pd.read_csv(args["trajectoryfile"], chunksize=number_of_lines_per_chunk, header=None)        
        chunk_number = 0
        #Loop over chunks
        for chunk in trajectory:    
            chunk_number = chunk_number + 1
            ddict.printLog('')
            ddict.printLog('Chunk %d of %d' % (chunk_number, number_of_chunks))
            #Divide the chunk into individual frames. If the chunk is the last chunk, the number of frames is different
            if chunk.shape[0] == number_of_lines_last_chunk:                            
                frames = np.split(chunk, last_chunk_size)
            else:
                frames = np.split(chunk, chunk_size)
            for frame in frames:
                #Delete the lines which are not needed for the analysis, depending on the file format
                frame, split_frame = run(frame, element_masses, molecule_id)
                #Species recognition, if requested                                                 
                if analysis_choice2_molrec == 'y':
                    split_frame['Species'] = species_id 
                #Drop all CNT and carbon_wall atoms as they are not needed for the analysis                                
                split_frame=split_frame[~split_frame['Molecule'].isin(CNTs)]            
                split_frame=split_frame[~split_frame['Molecule'].isin(carbon_walls)]
                #For analysis 1 and 2: All atoms are dropped which have a z coordinate larger than the maximum z coordinate of the CNT or smaller than the minimum z coordinate of the CNT.
                if analysis_choice2=='1' or analysis_choice2=='2':
                    split_frame=split_frame[split_frame['Z'].astype(float) <= max_z_CNT]   
                    split_frame=split_frame[split_frame['Z'].astype(float) >= min_z_CNT]
                if analysis_spec_molecule=='y':
                    #If a species is specified for analysis: Drop all atoms which are not of the species which is to be analyzed
                    split_frame=split_frame[split_frame['Species'].astype(int) == int(spec_molecule)]                
                    
                #RADIAL DENSITY FUNCTION                
                if analysis_choice2=='1':
                    #Finally calculate the radial density function with the remaining atoms.
                    split_frame['X_adjust']=split_frame['X'].astype(float) - tubecenter_x
                    split_frame['Y_adjust']=split_frame['Y'].astype(float) - tubecenter_y
                    split_frame['Distance']=np.sqrt(split_frame['X_adjust']**2 + split_frame['Y_adjust']**2)
                    split_frame['Distance_bin']=pd.cut(split_frame['Distance'], bins=raddens_bin_edges, labels=raddens_bin_labels+1)
                    #Finally add all masses of the atoms in each bin to the corresponding bin.
                    raddens_df_temp=split_frame.groupby(pd.cut(split_frame['Distance'], raddens_bin_edges))['Mass'].sum().reset_index(name='Weighted_counts')
                    raddens_df_temp=pd.DataFrame(raddens_df_temp)
                    #Add a new first column with the index+1 of the bin.
                    raddens_df_temp.insert(0, 'Bin', raddens_df_temp.index+1)                    
                    #Finally we need to write the results into the raddens_df dataframe. The row is defined by the frame number.
                    for i in range(num_increments):
                        raddens_df.loc[counter,'Bin %d' % (i+1)]=raddens_df_temp.loc[i,'Weighted_counts']
                    #Remove the raddens_df_temp dataframe every loop.
                    del raddens_df_temp                                                                      

                #Z-DENSITY PROFILE
                if analysis_choice2=='3':
                    #Calculate the z density function with the remaining atoms.
                    split_frame['Z_bin']=pd.cut(split_frame['Z'].astype(float), bins=z_bin_edges, labels=z_bin_labels+1)
                    #Add all masses of the atoms in each bin to the corresponding bin.
                    zdens_df_temp=split_frame.groupby(pd.cut(split_frame['Z'].astype(float), z_bin_edges))['Mass'].sum().reset_index(name='Weighted_counts')
                    zdens_df_temp=pd.DataFrame(zdens_df_temp)
                    #Add a new first column with the index+1 of the bin.
                    zdens_df_temp.insert(0, 'Bin', zdens_df_temp.index+1)                    
                    #Write the results into the zdens_df dataframe. The row is defined by the frame number.
                    for i in range(num_increments):
                        zdens_df.loc[counter,'Bin %d' % (i+1)]=zdens_df_temp.loc[i,'Weighted_counts']

                    #TUBE SECTION -> displacement for accessible volume
                    #All atoms are dropped which have a z coordinate larger than the maximum z coordinate of the CNT or smaller than the minimum z coordinate of the CNT.
                    split_frame=split_frame[split_frame['Z'].astype(float) <= max_z_CNT]   
                    split_frame=split_frame[split_frame['Z'].astype(float) >= min_z_CNT] 

                #Z-DENSITY PROFILE AND ACCESSIBLE VOLUME             
                if analysis_choice2=='2' or analysis_choice2=='3':    
                    #Identify the atom inside the tube, which is the most displaced from its center. Over the whole trajectory.
                    split_frame['X_adjust']=split_frame['X'].astype(float) - tubecenter_x
                    split_frame['Y_adjust']=split_frame['Y'].astype(float) - tubecenter_y
                    #While caluclating the distance, add the elements radius.
                    split_frame['Distance']=np.sqrt(split_frame['X_adjust']**2 + split_frame['Y_adjust']**2) + split_frame['Atom'].map(element_radii)
                    if split_frame['Distance'].max() > maxdisp_atom_dist:
                        maxdisp_atom_dist=split_frame['Distance'].max()
                        #Save the complete info of the most displaced atom.
                        maxdisp_atom_row=split_frame.loc[split_frame['Distance'].idxmax()]

                counter+=1
                print('Frame %d of %d' % (counter, number_of_frames), end='\r')
        ddict.printLog('')
        ddict.printLog('Finished processing the trajectory. %d frames were processed.' % (counter))
        ddict.printLog('')
        
        #DATA PROCESSING
        #RADIAL DENSITY FUNCTION - PART 2 
        if analysis_choice2=='1':
            post.raddens_func(raddens_df, raddens_bin_edges, number_of_frames, length_CNT, tuberadius)
        if analysis_choice2=='2':
            post.acc_vol_func(maxdisp_atom_row, length_CNT)
        if analysis_choice2=='3':
            post.zdens_func(zdens_df, z_bin_edges, number_of_frames, maxdisp_atom_row, length_CNT, tuberadius, num_increments, min_z_CNT, max_z_CNT, simbox_x, simbox_y, tubecenter_z)


        ddict.printLog("The main loop took %0.3f seconds to run" % (time.time() - Main_time))