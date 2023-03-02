#!/usr/bin/python3.8
# The program is written by Leonard Dick, 2022

# MODULES
import sys
import time
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
import pandas as pd
from typing import Tuple
from scipy.spatial.distance import cdist
#----------#
import defdict as ddict

#ARGUMENTS
args=ddict.read_commandline()

#FUNCTIONS
#Read the trajectory file and return the number of atoms, the box dimensions and the trajectory as a dataframe.
def trajfile_read() -> Tuple[int, float, float, float, pd.DataFrame, pd.DataFrame]:
    #MAIN
    ddict.printLog_red('TRAJECTORY INFO mode')
    ddict.printLog('')
    ddict.printLog('Trajectory file: ' + args["trajectoryfile"])
    #Check if the filename contains '.pdb' or '.xyz' of '.lmp'. For now just .xyz is possible.
    if args["trajectoryfile"].endswith(".xyz"):
    #Read the box dimensions of the simulation box from terminal input.
        ddict.printLog("The trajectory file is in xyz format.")
        simbox_x=float(ddict.get_input('Enter the x dimension of the simulation box [Ang]: '))
        simbox_y=float(ddict.get_input('Enter the y dimension of the simulation box [Ang]: '))
        simbox_z=float(ddict.get_input('Enter the z dimension of the simulation box [Ang]: '))
        ddict.printLog('The simulation box dimensions are [Ang]: ' + str(simbox_x) + ' x ' + str(simbox_y) + ' x ' + str(simbox_z))
        ddict.printLog('')
        #The number of atoms is given in the first line of the xyz file.
        with open(args["trajectoryfile"]) as f:
            num_atoms=int(f.readline())
        ddict.printLog('There are %d atoms in each frame.' % (num_atoms))
        framechoice=ddict.get_input('Choose frame to check for carbon structures and molecules. [default=1]: ')
        if framechoice=='':
            framechoice=1
        else:
            framechoice=int(framechoice)
        #number of lines to skip
        skip_lines=2+(num_atoms+2)*(framechoice-1)
        #Read in the id frame of the trajectory file to a dataframe.
        id_frame=pd.read_csv(args["trajectoryfile"], sep='\s+', nrows=num_atoms, header=None, skiprows=skip_lines, names=['Atom', 'x', 'y', 'z'])
        check_frame=pd.read_csv(args["trajectoryfile"], sep='\s+', nrows=num_atoms, header=None, skiprows=skip_lines+(num_atoms+2)*2, names=['Atom', 'x', 'y', 'z'])
        #Additional column to check if atom is frozen or not.
        id_frame['Frozen'] = id_frame.eq(check_frame).all(1)
        #reset the index of the dataframe.
        id_frame.reset_index(drop=True, inplace=True)
        #create empty dataframe for the atom types and molecule numbers.
        charge_molident=pd.DataFrame()

    elif args["trajectoryfile"].endswith(".pdb"):
        ddict.printLog("The trajectory file is in pdb format.")
        #read the box dimensions from the first line of the file.
        with open(args["trajectoryfile"]) as f:
            box_info=f.readline()
            #split the line into a list of strings.
            box_info=box_info.split()
            #convert the strings to floats.
            simbox_x=float(box_info[1])
            simbox_y=float(box_info[2])
            simbox_z=float(box_info[3])
            ddict.printLog('The simulation box dimensions are [Ang]: ' + str(simbox_x) + ' x ' + str(simbox_y) + ' x ' + str(simbox_z))
            #check when 'CRYST1' appears the second time in the file.
            for i, line in enumerate(f):
                if 'CRYST1' in line:
                    second_line=i
                    break
        #the total number of atoms is the line number minus 1. PDB files have two extra lines, while the loop starts at 0.
        num_atoms=second_line-1
        ddict.printLog('There are %d atoms in each frame.' % (num_atoms))
        #lines per frame is the number of atoms plus 2.
        lines_per_frame=num_atoms+2
        framechoice=ddict.get_input('Choose frame to check for carbon structures and molecules. [default=1]: ')
        if framechoice=='':
            framechoice=1
        else:
            framechoice=int(framechoice)
        #number of lines to skip
        skip_lines=1+lines_per_frame*(framechoice-1)
        #read the first frame into a dataframe. Just consider columns 3 as the atom type, 4 as the molecule and 5 6 7 for the (x y z) position.
        id_frame=pd.read_csv(args["trajectoryfile"], sep='\s+', nrows=num_atoms, header=None, skiprows=skip_lines, names=['Atom', 'Molecule', 'x', 'y', 'z',  'Charge'])
        #read the third frame into a dataframe.
        check_frame=pd.read_csv(args["trajectoryfile"], sep='\s+', nrows=num_atoms, header=None, skiprows=skip_lines+lines_per_frame*2,  names=['Atom', 'Molecule', 'x', 'y', 'z',  'Charge'])
        #make a new dataframe with the charge and molecule number.
        charge_molident=id_frame[['Charge', 'Molecule']]
        #now drop these 2 columns.
        id_frame.drop(['Charge', 'Molecule'], axis=1, inplace=True)
        check_frame.drop(['Charge', 'Molecule'], axis=1, inplace=True)
        #Additional column to check if atom is frozen or not.
        id_frame['Frozen'] = id_frame.eq(check_frame).all(1)
        #in the atom column, the atom type is given as 'C1' or 'C2'. We want to only have the element name. Drop all numbers.
        id_frame['Atom'] = id_frame['Atom'].str.replace(r'\d+', '') 
        #remove all letters but the first one from the atom column.
        id_frame['Atom'] = id_frame['Atom'].str[0]
        #reset the index of the dataframe.
        id_frame.reset_index(drop=True, inplace=True)
    
    elif args["trajectoryfile"].endswith(".lmp"):
        ddict.printLog("The trajectory file is in lmp format.")
        #read the box dimensions from the lines 6,7 and 8 of the file.
        with open(args["trajectoryfile"]) as f:
            #skip the first 5 lines.
            for i in range(3):
                next(f)
            #The number ofg atoms is printed in line 4.
            num_atoms=int(f.readline())
            #skip the next line.
            next(f)
            #Read the box dimensions given in line 6, 7 and 8. Each dimension is the last number in the line.
            simbox_x_line=f.readline()
            simbox_x=float(simbox_x_line.split()[-1])
            simbox_y_line=f.readline()
            simbox_y=float(simbox_y_line.split()[-1])
            simbox_z_line=f.readline()
            simbox_z=float(simbox_z_line.split()[-1])
            ddict.printLog('The simulation box dimensions are [Ang]: ' + str(simbox_x) + ' x ' + str(simbox_y) + ' x ' + str(simbox_z))
            ddict.printLog('')
        ddict.printLog('There are %d atoms in each frame.' % (num_atoms))
        framechoice=ddict.get_input('Choose frame to check for carbon structures and molecules. [default=1]: ')
        if framechoice=='':
            framechoice=1
        else:
            framechoice=int(framechoice)
        #lines per frame
        lines_per_frame=num_atoms+9
        #number of lines to skip
        skip_lines=9+(lines_per_frame)*(framechoice-1)
        #Read in the id frame of the trajectory file to a dataframe.
        id_frame=pd.read_csv(args["trajectoryfile"], sep='\s+', nrows=num_atoms, header=None, skiprows=skip_lines, names=['Atom No.', 'Atom', 'x', 'y', 'z'])
        check_frame=pd.read_csv(args["trajectoryfile"], sep='\s+', nrows=num_atoms, header=None, skiprows=skip_lines+lines_per_frame*2, names=['Atom No.', 'Atom', 'x', 'y', 'z'])
        #drop the Atom No. column.
        id_frame.drop(['Atom No.'], axis=1, inplace=True)
        check_frame.drop(['Atom No.'], axis=1, inplace=True)
        #just keep the first letter of the atom column.
        id_frame['Atom'] = id_frame['Atom'].str[0]
        #additional column to check if atom is frozen or not.
        id_frame['Frozen'] = id_frame.eq(check_frame).all(1)
        #reset the index of the dataframe.
        id_frame.reset_index(drop=True, inplace=True)
        #create empty dataframe for the atom types and molecule numbers.
        charge_molident=pd.DataFrame()        
    else:
        ddict.printLog("The file is not in a known format. Use the help flag (-h) for more information")
        sys.exit()
    ddict.printLog('')
    #Write the number of frozen atoms to the terminal.
    ddict.printLog('The number of frozen atoms is %d.' % (id_frame['Frozen'].sum()))    

    #return the number of atoms, the simulation box dimensions and the id_frame dataframe.
    return num_atoms, simbox_x, simbox_y, simbox_z, id_frame, charge_molident


#Identify the carbon structures in the first frame.
def CNT_identification(id_frame) -> Tuple[pd.DataFrame, float, float, float, float, list, float, float, float, list, list, pd.DataFrame]:
    #CNT IDENTIFICATION
    #define the maximal C-C distance valid to condider two atoms as part of the same molecule.
    maxC_C_dist=1.7
    CNT_time = time.time()
    ddict.printLog('')
    ddict.printLog('The chosen frame is used to identify the carbon structures.')
    #Add a column to the dataframe that will contain the molecule number, with all values set to 0.
    id_frame['Molecule']=0
    #Start by making a new dataframe with only the frozen atoms.                                                           
    frozen_atoms=id_frame[id_frame['Frozen']==True]                         
    frozen_atoms=pd.DataFrame(frozen_atoms)   
    #check if the frozen atoms dataframe is empty. If it is, the simulation is not frozen.
    if frozen_atoms.empty==True:
        ddict.printLog('No frozen atoms were found. The analysis stops here...')
        sys.exit()
    #Add a column to the dataframe that will contain the index of the atom.                                          
    frozen_atoms['Index']=frozen_atoms.index                                            
    frozen_atoms=frozen_atoms.reset_index(drop=True)
    #make an array with the distances between all the frozen atoms.
    positions=np.empty((frozen_atoms.shape[0],3))
    #fill the array with the x, y and z positions of the atoms.                                       
    positions[:,0]=frozen_atoms['x'].values                                             
    positions[:,1]=frozen_atoms['y'].values
    positions[:,2]=frozen_atoms['z'].values
    #set up the distance matrix
    distance_matrix=cdist(positions, positions)                                         
    while frozen_atoms[(frozen_atoms['Molecule']==0)].empty==False:
        #first atom with mol number 0
        first_atom=frozen_atoms[(frozen_atoms['Molecule']==0)].iloc[0] 
        #assign new mol number to first atom with 0 in molecule column.                 
        frozen_atoms.loc[first_atom.name,'Molecule']=frozen_atoms['Molecule'].max()+1
        #add atom to the list of atoms that are already assigned to a molecule.   
        current_structure=[first_atom.name]                                             
        while len(current_structure)>0:
            #pop the first atom from the list of atoms that are already assigned to a molecule. The pop function also removes the atom from the current_structure list simultaneously.                                                 
            current_atom=current_structure.pop(0)
            for atom in frozen_atoms[(frozen_atoms['Molecule']==0)].index:
                if distance_matrix[current_atom,atom]<maxC_C_dist:
                    current_structure.append(atom)
                    frozen_atoms.loc[atom,'Molecule']=frozen_atoms['Molecule'].max()
    frozen_atoms=frozen_atoms.set_index('Index')
    id_frame.loc[frozen_atoms.index,'Molecule']=frozen_atoms['Molecule'].values
    #Check if the structure is either a CNT or a carbon wall. If all atoms of the same structure have the same z coordinate (within a certain tolerance), it is a carbon wall. If not, it is a CNT.
    #The z coordinates can have some uncertainties, tolerance set to 0.1  
    tolerance_z=0.1
    CNTs=[]                                                                             
    carbon_walls=[]
    #identify the carbon structures via the molecule number and the z coordinates of the atoms.
    for structure in range(1,frozen_atoms['Molecule'].max()+1):
        if id_frame[id_frame['Molecule']==structure]['z'].max()-id_frame[id_frame['Molecule']==structure]['z'].min()<tolerance_z:    
            carbon_walls.append(structure)
            ddict.printLog('Structure %d is a carbon wall.' % (structure))
        else:
            CNTs.append(structure)   
            ddict.printLog('Structure %d is a CNT.' % (structure))
    ddict.printLog('')
    ddict.printLog('Total number of carbon structures found: %d' % (len(CNTs)+len(carbon_walls)))
    #Time check
    ddict.printLog("The carbon structure identification took %0.3f seconds to run" % (time.time() - CNT_time))
    ddict.printLog('')
    counter=0
    carbon_wall_z=[]     
    #Get z coordinates/positions of the carbon walls.                                                                   
    for structure in carbon_walls:
        counter+=1
        ddict.printLog('The z coordinate of carbon wall %d is %0.3f.' % (counter, id_frame[id_frame['Molecule']==structure]['z'].mean()))
        carbon_wall_z.append(id_frame[id_frame['Molecule']==structure]['z'].mean())
    ddict.printLog('')
    #Get minimal and maximal z coordinate of the CNTs.
    min_z_CNT=id_frame[id_frame['Molecule'].isin(CNTs)]['z'].min()    
    max_z_CNT=id_frame[id_frame['Molecule'].isin(CNTs)]['z'].max()
    ddict.printLog('The minimal z coordinate of the pore is %0.3f.' % (min_z_CNT))
    ddict.printLog('The maximal z coordinate of the pore is %0.3f.' % (max_z_CNT))
    ddict.printLog('')
    tubecenter_z=(min_z_CNT+max_z_CNT)/2
    length_CNT=max_z_CNT-min_z_CNT

    #Tube in CNT identification.
    #Pull the xyz coordinates of the CNT atoms from the dataframe.
    CNT_atoms=id_frame[id_frame['Molecule'].isin(CNTs)] 
    #Drop the 'Molecule', 'Frozen' columns from the dataframe.              
    CNT_atoms=CNT_atoms.drop(['Molecule','Frozen'], axis=1)
    #Identify the atoms in CNT_atoms, which are closest to the tubecenter_z coordinate in z direction.                   
    CNT_atoms['z_distance']=abs(CNT_atoms['z']-tubecenter_z)
    #sort the atoms by the z_distance column.                   
    CNT_atoms=CNT_atoms.sort_values(by=['z_distance'])
    #The 0.02 are a tolerance, if not all atoms of the same ring are exactly at the same z coordinate.
    lowest_z=CNT_atoms.iloc[0]['z_distance']+0.02
    #Make a new dataframe with just the atoms having z_distance=<lowest_z.                        
    CNT_ring=CNT_atoms[CNT_atoms['z_distance']<=lowest_z]
    #drop the z_distance column from the CNT_atoms dataframe.                   
    CNT_atoms=CNT_atoms.drop(['z_distance'], axis=1)
    #Delete all atoms in the CNT_ring dataframe, which are more than 0.1 angstrom away in z direction from the first atom in the CNT_ring dataframe.
    CNT_ring=CNT_ring[CNT_ring['z']<=CNT_ring.iloc[0]['z']+0.1]             
    CNT_ring=CNT_ring[CNT_ring['z']>=CNT_ring.iloc[0]['z']-0.1]
    #Calculate the average x and y coordinate of the atoms in the CNT_ring dataframe.
    tubecenter_x=CNT_ring['x'].mean()                                       
    tubecenter_y=CNT_ring['y'].mean()
    #Get the number of rows in the CNT_ring dataframe.
    CNT_ring_rows=CNT_ring.shape[0]                                         
    ddict.printLog('The pore in the system contains a (%d,0) CNT.' % (CNT_ring_rows))
    ddict.printLog('The centerpoint of the CNT is at (%0.3f|%0.3f|%0.3f).' % (tubecenter_x, tubecenter_y, tubecenter_z))
    #The tubes radius is defined as the distance between the centerpoint of the CNT and the first atom in the CNT_ring dataframe, just for the x and y coordinate.
    tuberadius=((CNT_ring.iloc[0]['x']-tubecenter_x)**2+(CNT_ring.iloc[0]['y']-tubecenter_y)**2)**0.5     
    ddict.printLog('The radius of the CNT is %0.3f Ang.' % (tuberadius))
    ddict.printLog('')
    #Calculate the xy-distance of the centerpoint of the CNT to all CNT_atoms. If they are smaller/equal as the tuberadius, they belong to the tube
    CNT_atoms['xy_distance']=((CNT_atoms['x']-tubecenter_x)**2+(CNT_atoms['y']-tubecenter_y)**2)**0.5    
    #The 0.05 is a tolerance, if the xy_distance is not exactly the same as the tuberadius.
    tube_atoms=CNT_atoms[CNT_atoms['xy_distance']<=tuberadius+0.05]         

    return tube_atoms, tuberadius, tubecenter_x, tubecenter_y, tubecenter_z, carbon_wall_z, length_CNT, max_z_CNT, min_z_CNT, CNTs, carbon_walls, CNT_atoms


#function to edit the frame in pdb format
def pdb(frame, element_masses, molecule_id) -> Tuple[pd.DataFrame, pd.DataFrame]:
    #drop the first and last line
    frame=frame.drop(frame.index[[0,len(frame)-1]])
    #Split the frame into columns  and label them
    split_frame = frame[0].str.split(expand=True) 
    split_frame.columns = ['Label', 'Count', 'Atom', 'Molecule', 'X', 'Y', 'Z', 'Charge']
    #drop label and count column
    split_frame=split_frame.drop(['Label', 'Count'], axis=1)
    #In atom column, just keep the first character
    split_frame['Atom'] = split_frame['Atom'].str[0]
    #Add a new column with the atomic mass of each atom.
    split_frame['Mass'] = split_frame['Atom'].map(element_masses)
    #Add the molecule number to the frame       
    split_frame['Molecule'] = molecule_id 
    return frame, split_frame

#function to edit the frame in xyz format
def xyz(frame, element_masses, molecule_id) -> Tuple[pd.DataFrame, pd.DataFrame]:
    #drop the first two lines
    frame=frame.drop(frame.index[[0,1]])
    #Split the frame into columns and label them                              
    split_frame = frame[0].str.split(expand=True)                       
    split_frame.columns = ['Atom', 'X', 'Y', 'Z']
    #Add a new column with the atomic mass of each atom.                          
    split_frame['Mass'] = split_frame['Atom'].map(element_masses)
    #Add the molecule number to the frame       
    split_frame['Molecule'] = molecule_id 
    return frame, split_frame

#function to edit the frame in lammps format
def lammpstrj(frame, element_masses, molecule_id) -> Tuple[pd.DataFrame, pd.DataFrame]:
    #drop the first nine lines
    frame=frame.drop(frame.index[[0,1,2,3,4,5,6,7,8]])
    #Split the frame into columns and label them                              
    split_frame = frame[0].str.split(expand=True)                       
    split_frame.columns = ['Atom Count', 'Atom', 'X', 'Y', 'Z']
    #drop the Atom count column
    split_frame=split_frame.drop(['Atom Count'], axis=1)
    #In atom column, just keep the first character
    split_frame['Atom'] = split_frame['Atom'].str[0]
    #Add a new column with the atomic mass of each atom.                          
    split_frame['Mass'] = split_frame['Atom'].map(element_masses)
    #Add the molecule number to the frame       
    split_frame['Molecule'] = molecule_id 
    return frame, split_frame
