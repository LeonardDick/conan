# The program is written by Leonard Dick, 2023

# MODULES
import MDAnalysis as mda
import networkx as nx
from collections import Counter
import pandas as pd
import networkx as nx
from rdkit import Chem
from rdkit.Chem import Draw
import re
import sys
import numpy as np
from typing import Tuple
from prettytable import PrettyTable
from scipy.spatial import cKDTree
# ----- Own modules ----- #
import defdict as ddict

#ARGUMENTS
args = ddict.read_commandline()

#FUNCTIONS
def read_first_frame(file) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, tuple]:

    # If file is in XYZ format:
    if file.endswith('.xyz'):

        # Create a Universe from the MD trajectory file.
        u = mda.Universe(file)
        # Load the first frame
        u.trajectory[0]
        # Get all atom information for the first frame.
        atoms = u.atoms
        data_list = []
        df_list = []

        # Print atom information for molecule recognition.
        for atom in atoms: 
            entry = {'element': str(atom.name), 'x': float(atom.position[0]), 'y': float(atom.position[1]), 'z': float(atom.position[2])}
            data_list.append(entry)
            atom_name = atom.name
            atom_x = atom.position[0]
            atom_y = atom.position[1]
            atom_z = atom.position[2]

            # Extract atom charge and molecule, if it's not available in your universe it will raise an error.
            try:
                atom_molecule = atom.resname  
                atom_charge = atom.charge
            except AttributeError:
                atom_charge = None
                atom_molecule = None
            df_list.append([atom_name, atom_molecule, atom_x, atom_y, atom_z, atom_charge])

        id_frame = pd.DataFrame(df_list, columns=['Atom', 'Molecule', 'x', 'y', 'z', 'Charge'])
        # Add the label column.
        id_frame['Label'] = None
        # Rename Atom column to Element.
        id_frame.rename(columns={'Atom': 'Element'}, inplace=True)
        # Order the columns in the dataframe by element x y z molecule charge label.
        id_frame=id_frame[['Element', 'x', 'y', 'z', 'Molecule', 'Charge', 'Label']]
        
        # Load the second frame.
        u.trajectory[1]
        atoms_F2 = u.atoms
        data_list_F2 = []
        df_list_F2 = []

        # Print atom information for molecule recognition and append it to a list.
        for atom in atoms_F2: 
            entry = {'element': str(atom.name), 'x': float(atom.position[0]), 'y': float(atom.position[1]), 'z': float(atom.position[2])}
            data_list_F2.append(entry)
            atom_name = atom.name
            atom_x = atom.position[0]
            atom_y = atom.position[1]
            atom_z = atom.position[2]
            df_list_F2.append([atom_name, atom_x, atom_y, atom_z])
        
        # Save the second frame in a dataframe
        id_frame2 = pd.DataFrame(df_list_F2, columns=['Element', 'x', 'y', 'z' ])

        # Simulation box dimensions
        ddict.printLog('Enter the dimensions of the simulation box [Ang]:')
        simbox_x = float(ddict.get_input('[X]   ', args))
        simbox_y = float(ddict.get_input('[Y]   ', args))
        simbox_z = float(ddict.get_input('[Z]   ', args))
        

    # If file is in PDB format:
    elif file.endswith('.pdb'):

        # Read the box dimensions from the first line of the file.
        with open(file) as f:
            second_line=0
            box_info=f.readline()

            # Split the line into a list of strings.
            box_info=box_info.split()
            # Convert the strings to floats.
            simbox_x=float(box_info[1])
            simbox_y=float(box_info[2])
            simbox_z=float(box_info[3])

            # Check when 'CRYST1' appears the second time in the file.
            for i, line in enumerate(f):
                if 'CRYST1' in line:
                    second_line=i
                    break

        # The total number of atoms is the line number minus 1. PDB files have two extra lines, while the loop starts at 0.
        num_atoms=second_line-1
        # The lines per frame is the number of atoms plus 2.
        lines_per_frame=num_atoms+2
        # Number of lines to skip
        skip_lines=1+lines_per_frame
        # Read the first frame into a dataframe. Just consider columns 3 as the atom type, 4 as the molecule and 5 6 7 for the (x y z) position.
        id_frame=pd.read_csv(args["trajectoryfile"], sep='\s+', nrows=num_atoms, header=None, skiprows=1, names=['Atom', 'Molecule', 'x', 'y', 'z',  'Charge'])
        # Make a new column with just the first letter of the atom name and all consecutive small letters.
        id_frame['Element'] = id_frame['Atom'].str[0]
        # Rename the Atom column to label column
        id_frame.rename(columns={'Atom': 'Label'}, inplace=True)
        # Order the columns in the dataframe by element x y z molecule charge label.
        id_frame=id_frame[['Element', 'x', 'y', 'z', 'Molecule', 'Charge', 'Label']]

        # Also read the second frame into a dataframe.
        id_frame2=pd.read_csv(args["trajectoryfile"], sep='\s+', nrows=num_atoms, header=None, skiprows=skip_lines, names=['Atom', 'Molecule', 'x', 'y', 'z',  'Charge'])
        # Make a new column with just the first letter of the atom name and all consecutive small letters.
        id_frame2['Element'] = id_frame2['Atom'].str[0]

        #reset the index of the dataframes.
        id_frame.reset_index(drop=True, inplace=True)
        id_frame2.reset_index(drop=True, inplace=True)

        #convert the first dataframe to a list of dictionaries.
        data_list = []
        for index, row in id_frame.iterrows():
            entry = {'element': str(row['Element']), 'x': float(row['x']), 'y': float(row['y']), 'z': float(row['z'])}
            data_list.append(entry)
    

    # If file is in LAMMPS format:
    elif file.endswith('.lmp') or file.endswith('.lammpstrj'):
        with open(file, 'r') as f:
           # Skip the first 3 lines.
            for i in range(3):
                next(f)
            # The number of atoms is printed in line 4.
            num_atoms=int(f.readline())
            # Skip the next line.
            next(f)

            #Read the box dimensions given in line 6, 7 and 8. Each dimension is the last number in the line.
            simbox_x_line=f.readline()
            simbox_x=float(simbox_x_line.split()[-1])
            simbox_y_line=f.readline()
            simbox_y=float(simbox_y_line.split()[-1])
            simbox_z_line=f.readline()
            simbox_z=float(simbox_z_line.split()[-1])
            read_one_frame=f.readline()

            # Read the atom information for the first frame
            data_list = []
            df_list = []
            for line in f:
                if line.strip() == "ITEM: TIMESTEP":
                    break
                entries = line.split()
                atom_id = entries[0]
                atom_type = entries[1]
                atom_x = float(entries[2])
                atom_y = float(entries[3])
                atom_z = float(entries[4])
                entry = {'element': str(atom_type), 'x': float(atom_x), 'y': float(atom_y), 'z': float(atom_z)}
                data_list.append(entry)
                
                df_list.append([atom_type, atom_x, atom_y, atom_z])
            
            # Create a DataFrame from the collected data
            id_frame = pd.DataFrame(df_list, columns=['type', 'x', 'y', 'z'])
            
            # Add additional columns if available in your LAMMPS file
            id_frame['Molecule'] = None
            id_frame['Charge'] = None
            id_frame['Label'] = None
            id_frame.rename(columns={'type': 'Element'}, inplace=True)            

            # Order the columns in the dataframe by element x y z molecule charge label.
            id_frame=id_frame[['Element', 'x', 'y', 'z', 'Molecule', 'Charge', 'Label']]

            # Read the atom information for the second frame
            for i in range(8):
                next(f)
            data_list = []
            df_list = []
            for line in f:
                if line.strip() == "ITEM: TIMESTEP":
                    break
                entries = line.split()
                atom_id = entries[0]
                atom_type = entries[1]
                atom_x = float(entries[2])
                atom_y = float(entries[3])
                atom_z = float(entries[4])
                entry = {'element': str(atom_type), 'x': float(atom_x), 'y': float(atom_y), 'z': float(atom_z)}
                data_list.append(entry)

                df_list.append([atom_type, atom_x, atom_y, atom_z])

            # Create a DataFrame from the collected data
            id_frame2 = pd.DataFrame(df_list, columns=['Element', 'x', 'y', 'z'])

    else:
        ddict.printLog("The file is not in a known format. Use the help flag (-h) for more information")
        sys.exit()

    ddict.printLog('')

    # Check which atoms did not move, they should have the same x y z coordinates. Label them as True in the Struc column.
    id_frame['Struc'] = (id_frame['x'] == id_frame2['x']) & (id_frame['y'] == id_frame2['y']) & (id_frame['z'] == id_frame2['z'])

    # Write the box dimensions as a tuple.
    box_size=(simbox_x, simbox_y, simbox_z)
    ddict.printLog(f'The simulation box dimensions are [Ang]: {float(box_size[0]):.3f} x {float(box_size[1]):.3f} x {float(box_size[2]):.3f}')
    ddict.printLog(f"\nTotal number of atoms: {len(data_list)}\n")

    return data_list, id_frame, id_frame2, box_size

# Function to calculate the minimum distance between two positions.
def minimum_image_distance(position1, position2, box_size) -> float:
    # Calculate the minimum distance between two positions.
    d_pos = position1 - position2
    d_pos = d_pos - box_size*np.round(d_pos/box_size)
    return np.linalg.norm(d_pos)

# Function to identify molecular bonds from a distance search using a k-d tree.
def identify_molecules_and_bonds(atoms, box_size, cutoff = 3.0) -> Tuple[list, list]:

    # Get the covalent radii 
    covalent_radii = ddict.dict_covalent()

    # Define bond_distances
    bond_distances = {(e1, e2): (covalent_radii[e1] + covalent_radii[e2]) * 1.15 for e1 in covalent_radii for e2 in covalent_radii}

    # Create a graph with atoms as nodes and bonds as edges
    simbox_G = nx.Graph()    

    atom_positions = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms]) % box_size
    atom_elements = [atom['element'] for atom in atoms]

    # Create k-d tree for efficient search
    tree = cKDTree(atom_positions, boxsize=box_size)

    # Find pairs within max bond_distance
    pairs = tree.query_pairs(max(bond_distances.values()))

    for i, j in pairs:
        bond_distance = bond_distances.get((atom_elements[i], atom_elements[j]), float('inf'))
        
        # Correct the distance considering the minimum image convention
        distance = minimum_image_distance(atom_positions[i], atom_positions[j], box_size)
        
        if distance <= bond_distance:
            # Add an edge in the graph if atoms are bonded
            simbox_G.add_edge(i, j)

    # Each connected component in the graph represents a molecule
    molecules = [[atoms[i]['Atom'] for i in molecule] for molecule in nx.connected_components(simbox_G)]

    # Determine bonds for each molecule
    molecule_bonds = []
    for molecule in molecules:
        bonds = [sorted((i, j)) for i, j in simbox_G.edges(molecule)]  # Only include bonds within the current molecule
        molecule_bonds.append(bonds)

    return molecules, molecule_bonds

# Structure recognition section.
def structure_recognition(id_frame, box_size) -> Tuple[pd.DataFrame, list, list, list, list, list, list, list, list]:
    # Identify all solid structures in the simulation box and identify if it is a wall or a pore.
    # A wall extends in two dimensions, a pore in all three dimensions.
    # First make a new dataframe with just the structure atoms.
    id_frame["Molecule"] = None
    structure_frame = id_frame[id_frame['Struc'] == True]
    #convert the first dataframe to a list of dictionaries. We also need to store the atom index.

    str_atom_list = []
    for index, row in structure_frame.iterrows():
        entry = {'Atom': index, 'element': str(row['Element']), 'x': float(row['x']), 'y': float(row['y']), 'z': float(row['z'])}
        str_atom_list.append(entry)

    # Identify the structures and the bonds
    molecules_struc, molecule_bonds_struc = identify_molecules_and_bonds(str_atom_list, box_size)


    for i, molecule in enumerate(molecules_struc):
        id_frame.loc[molecule, "Molecule"] = i + 1
        structure_frame.loc[molecule, "Molecule"] = i + 1

    # Identify the wall and the pore
    CNTs = []
    counter_pore = 0
    Walls = []
    Walls_positions = []
    counter_wall = 0
    # Make a copy of the structure frame (to assure pandas treats it as a copy, not a view)
    structure_frame_copy = structure_frame.copy()

    # Consider all Molecules in the structure frame and get the maximum and minimum x, y and z coordinates for each respective one.
    # If the difference in x, y and z is larger than 1.0, it is a pore. If it is smaller in one direction, it is a wall.
    for i in range(1, len(molecules_struc) + 1):
        x_max = structure_frame.loc[structure_frame['Molecule'] == i, 'x'].max()
        x_min = structure_frame.loc[structure_frame['Molecule'] == i, 'x'].min()
        y_max = structure_frame.loc[structure_frame['Molecule'] == i, 'y'].max()
        y_min = structure_frame.loc[structure_frame['Molecule'] == i, 'y'].min()
        z_max = structure_frame.loc[structure_frame['Molecule'] == i, 'z'].max()
        z_min = structure_frame.loc[structure_frame['Molecule'] == i, 'z'].min()
        if (x_max - x_min) > 1.0 and (y_max - y_min) > 1.0 and (z_max - z_min) > 1.0:
            counter_pore += 1
            ddict.printLog(f"Structure {i} is a pore, labeled Pore{counter_pore}\n")

            # Change the structure column to pore{i}
            structure_frame_copy.loc[structure_frame['Molecule'] == i, "Struc"] = f"Pore{counter_pore}"
            CNTs.append(f"Pore{counter_pore}")

        elif (x_max - x_min) < 1.0 or (y_max - y_min) < 1.0 or (z_max - z_min) < 1.0:
            counter_wall += 1
            ddict.printLog(f"Structure {i} is a wall, labeled Wall{counter_wall}")
            if (x_max - x_min) < 1.0:
                ddict.printLog(f"The wall extends in x direction at x = {x_min:.2f} Ang.\n")
            if (y_max - y_min) < 1.0:
                ddict.printLog(f"The wall extends in y direction at y = {y_min:.2f} Ang.\n")
            if (z_max - z_min) < 1.0:
                ddict.printLog(f"The wall extends in z direction at z = {z_min:.2f} Ang.\n")
                Walls_positions.append(z_min)
            structure_frame_copy.loc[structure_frame['Molecule'] == i, "Struc"] = f"Wall{counter_wall}"
            Walls.append(f"Wall{counter_wall}")

    # Copy the structure frame back to the original structure frame.
    structure_frame = structure_frame_copy

    # Finally put the structure information back into the original dataframe.
    id_frame.loc[structure_frame.index, "Struc"] = structure_frame["Struc"]

    # Exchange all the entries in the 'Struc' column saying 'False' with 'Liquid'.
    id_frame['Struc'].replace(False, 'Liquid', inplace=True)

    # Print the structure information .
    ddict.printLog(f"\nTotal number of structures: {len(molecules_struc)}")
    ddict.printLog(f"Number of walls: {len(Walls)}")
    ddict.printLog(f"Number of pores: {len(CNTs)}\n")

    
    # Pore section  
    # Classify the pores in the system. First we find the minimum and maximum z coordinates of the pores.
    min_z_pore = []
    max_z_pore = []
    length_pore = []
    center_pore = []
    CNT_centers = []
    tuberadii = []
    CNT_volumes = []
    CNT_atoms = []

    for i in range(1, len(CNTs) + 1):
        pore = id_frame[id_frame['Struc'] == f"Pore{i}"].copy()
        min_z_pore.append(pore['z'].min())
        max_z_pore.append(pore['z'].max())

        # The length of each pore is the difference between the maximum and minimum z coordinate.
        length_pore.append(max_z_pore[i - 1] - min_z_pore[i - 1])
        ddict.printLog(f"The length of Pore{i} is {length_pore[i - 1]:.2f} Ang.")

        # The center of each pore is the average of the minimum and maximum z coordinate.
        center_pore.append((max_z_pore[i - 1] + min_z_pore[i - 1]) / 2)

        # The pore consists of a CNT in the center. To classify the CNT, we should get its radius. For this we just take the atoms in the pore dataframe closest to the center of the pore (with a small tolerance).
        pore.loc[:, 'z_distance'] = abs(pore['z'] - center_pore[i - 1])
        pore = pore.sort_values(by = ['z_distance'])
        lowest_z = pore.iloc[0]['z_distance'] + 0.02
        CNT_ring = pore[pore['z_distance'] <= lowest_z].copy()

        # Delete all atoms in the CNT_ring dataframe, which are more than 0.1 angstrom away in the z direction from the first atom in the CNT_ring dataframe.
        CNT_ring.loc[:, 'z_distance'] = abs(CNT_ring['z'] - CNT_ring.iloc[0]['z'])
        CNT_ring = CNT_ring[CNT_ring['z_distance'] <= 0.1]

        # Calculate the average x and y coordinate of the atoms in the CNT_ring dataframe.
        x_center = CNT_ring['x'].mean()
        y_center = CNT_ring['y'].mean()
        ddict.printLog(f"The center of Pore{i} is at ({x_center:.2f}, {y_center:.2f}, {center_pore[i - 1]:.2f}) Ang.")
        # Combine the x, y and z centers to a numpy array.
        center = np.array([x_center, y_center, center_pore[i - 1]])
        CNT_centers.append(center)

        # Calculate the radius of the CNT_ring.
        tuberadius = np.sqrt((CNT_ring.iloc[0]['x'] - x_center) ** 2 + (CNT_ring.iloc[0]['y'] - y_center) ** 2)
        tuberadii.append(tuberadius)
        ddict.printLog(f"The radius of Pore{i} is {tuberadius:.2f} Ang.")

        # Calculate the volume of the CNT_ring.
        pore_volume = np.pi * tuberadius ** 2 * length_pore[i - 1]
        ddict.printLog(f"The volume of Pore{i} is {pore_volume:.2f} Ang^3.\n")
        CNT_volumes.append(pore_volume)
        
        #Calculate the xy-distance of the centerpoint of the CNT to all pore atoms. If they are smaller/equal as the tuberadius, they belong to the CNT.
        id_frame.loc[id_frame['Struc'] == f"Pore{i}", 'xy_distance'] = np.sqrt((id_frame.loc[id_frame['Struc'] == f"Pore{i}", 'x'] - x_center) ** 2 + (id_frame.loc[id_frame['Struc'] == f"Pore{i}", 'y'] - y_center) ** 2)
        id_frame.loc[id_frame['Struc'] == f"Pore{i}", 'xy_distance'] = id_frame.loc[id_frame['Struc'] == f"Pore{i}", 'xy_distance'].round(2)

        # Save the information about the CNT in the structure dataframe, by adding a new column 'CNT' with the CNT number, if the xy_distance is smaller/equal the tuberadius. 0.05 is a tolerance.
        id_frame.loc[id_frame['Struc'] == f"Pore{i}", 'CNT'] = 0
        id_frame.loc[(id_frame['Struc'] == f"Pore{i}") & (id_frame['xy_distance'] <= tuberadius + 0.05), 'CNT'] = i
        
        # Delete the xy_distance column again.
        id_frame.drop(columns=['xy_distance'], inplace=True)

    return id_frame, min_z_pore, max_z_pore, length_pore, CNT_centers, tuberadii, CNT_volumes, CNT_atoms, Walls_positions

# Molecule recognition section.
def molecule_recognition(id_frame, box_size) -> pd.DataFrame:

    # Convert the first dataframe to a list of dictionaries.
    str_liquid_list = []
    for index, row in id_frame.iterrows():
        entry = {'Atom': index, 'element': str(row['Element']), 'x': float(row['x']), 'y': float(row['y']), 'z': float(row['z'])}
        str_liquid_list.append(entry)

    # Identify the molecules and the bonds
    molecules, molecule_bonds = identify_molecules_and_bonds(str_liquid_list, box_size)

    # Translate the atom numbers in the molecules list of lists to the element symbols to a new list of lists.
    molecules_sym = []
    for molecule in molecules:
        molecule_symloop = []
        for atom in molecule:
            molecule_symloop.append(str_liquid_list[atom]['element'])
        molecules_sym.append(molecule_symloop)

    # Translate the atom numbers in the molecule_bonds list of lists to the element symbols to a new list of lists.
    molecule_bonds_sym = []
    for molecule in molecule_bonds:
        molecule_bonds_symloop = []
        for bond in molecule:
            molecule_bonds_symloop.append((str_liquid_list[bond[0]]['element'], str_liquid_list[bond[1]]['element']))
        molecule_bonds_sym.append(molecule_bonds_symloop)
    
    # Add the molecule information to the dataframe
    for i, molecule in enumerate(molecules):
        id_frame.loc[molecule, "Molecule"] = 1 + i

 
    # Initialize label counter and molecule counter
    label_counter = 1
    molecule_counter = id_frame['Molecule'][0]

    # Save original index to a new column
    id_frame['original_index'] = id_frame.index
    # To make sure this works if the atoms are sorted by atom type and not molecule wise, we need to resort the dataframe by molecule number.
    id_frame = id_frame.sort_values(by = ['Molecule', 'original_index'])

    for index, row in id_frame.iterrows():
        # Check if the molecule number has changed
        if row['Molecule'] != molecule_counter:
            # Reset the label counter if the molecule number has changed
            label_counter = 1
            molecule_counter = row['Molecule']
        # Create the label by combining the element name and the label counter
        label = row['Element'] + str(label_counter)
        # Store the label in the Label column
        id_frame.loc[index, 'Label'] = label
        # Increment the label counter
        label_counter += 1

    # Finally revert the sorting of the dataframe to the original order.
    id_frame = id_frame.sort_values(by='original_index').drop(columns='original_index')

    # Create a new dataframe for each molecule. 
    molecule_frame = pd.DataFrame(columns = ['Molecule', 'Atoms', 'Bonds', 'Atoms_sym', 'Bonds_sym'])
    molecule_frame['Molecule'] = range(1, len(molecules) + 1)
    molecule_frame['Atoms'] = molecules
    molecule_frame['Bonds'] = molecule_bonds
    molecule_frame['Atoms_sym'] = molecules_sym
    molecule_frame['Bonds_sym'] = molecule_bonds_sym    
    # Add another column with the Labels. We get the info by matching the atom numbers in the Atoms column with the roq index in the id_frame dataframe.
    molecule_frame['Labels'] = molecule_frame['Atoms'].apply(lambda x: [id_frame['Label'][i] for i in x])
    # Sort the lists in the Bonds_sym column to prepare to check for duplicates
    molecule_frame['Bonds_sym'] = molecule_frame['Bonds_sym'].apply(lambda x: sorted(x))

    # Make an independant copy of the molecule_frame dataframe name unique_molecule_frame
    unique_molecule_frame = molecule_frame.copy()

    # Drop all rows that are duplicates in the Bonds_sym column in a new dataframe
    unique_molecule_frame['Bonds_sym'] = unique_molecule_frame['Bonds_sym'].apply(tuple)
    unique_molecule_frame = unique_molecule_frame.drop_duplicates(subset=['Bonds_sym'])
    unique_molecule_frame = unique_molecule_frame.reset_index(drop=True)

    # Now get the chemical formulas of the unique molecules by simply counting the number of atoms of each element in the Atoms_sym column
    unique_molecule_frame['Molecule'] = unique_molecule_frame['Atoms_sym'].apply(lambda x: "".join(f"{element}{count}" for element, count in Counter(x).items()))


    #Now adjust the labels in molecule_frame to include the molecule kind.
    for i, row in unique_molecule_frame.iterrows():
        for i2, row2 in molecule_frame.iterrows():
            if sorted(row['Labels']) == sorted(row2['Labels']):
                new_labels = [f"{label}_{i + 1}" for label in molecule_frame.loc[i2, 'Labels']]
                molecule_frame.at[i2, 'Labels'] = new_labels
    
    # Change the Labels column in the id_frame dataframe to the new labels
    for i, row in molecule_frame.iterrows():
        id_frame.loc[row['Atoms'], 'Label'] = row['Labels']

    id_frame['Species']=id_frame['Label'].str.split('_').str[1]

    # Remove the _1/2/3.. from the Labels column in the id_frame dataframe
    id_frame['Label'] = id_frame['Label'].str.split('_').str[0]

    # Count how often each individual molecule/species occurs in the system
    # make a new dataframe with the columns 'Species' and 'Atom_count' and 'Molecule_count'
    molecule_count = pd.DataFrame()
    molecule_count['Atom_count'] = id_frame['Species'].value_counts()

    # Reset the index
    molecule_count.reset_index(inplace=True)

    # Now 'Species' is no longer an index, so we can safely create the column
    molecule_count.rename(columns={'index': 'Species'}, inplace=True)

    #molecule_count = molecule_count.sort_values(by='Species')
    molecule_count = molecule_count.sort_values(by='Species').reset_index(drop=True)

    # Now the number of molecules is the Atom_count divided by the number of atoms in each species. Therefore we need to get the number of atoms in each species. It is the length of the roe 'Atoms' in the unique_molecule_frame dataframe.
    molecule_count['Molecule_count'] = molecule_count['Atom_count'] / unique_molecule_frame['Atoms'].apply(len)


    # Set the index in the dataframes to +=1.
    unique_molecule_frame.index += 1
    molecule_count.index += 1


   # Print this information in a nice table
    table = PrettyTable()
    table.field_names = ["Species", "Chemical formula", "No. molecules", "No. atoms per molecule"]
    for i, row in unique_molecule_frame.iterrows():
        table.add_row([i, row['Molecule'], int(molecule_count['Molecule_count'][i]), len(row['Atoms'])])

    # Print the results
    ddict.printLog(" ")
    ddict.printLog(table)
    

    for i, row in unique_molecule_frame.iterrows():
        # Create a new graph if the molecule is smaller than 50
        if len(row['Atoms']) < 50:
            mol = nx.Graph()

            # Combine atoms and labels into a dictionary
            atom_labels = dict(zip(row['Atoms'], row['Labels']))

            for atom, label in atom_labels.items():
                mol.add_node(atom, element=label)

            for bond in row['Bonds']:
                mol.add_edge(*bond)

            rdkit_mol = Chem.RWMol()

            # Add atoms with labels
            atom_mapping = {}  # To keep track of atom mapping from graph to RDKit molecule

            for node in mol.nodes():
                atom_label = mol.nodes[node]['element']
                # Use just the element symbol
                element_symbol = re.findall('[A-Za-z]+', atom_label)[0]
                atom = Chem.Atom(element_symbol)
                atom_idx = rdkit_mol.AddAtom(atom)
                atom_mapping[node] = atom_idx
                rdkit_mol.GetAtomWithIdx(atom_idx).SetProp('atomNote', atom_label)

            # Add bonds
            for edge in mol.edges():
                atom_idx1 = atom_mapping[edge[0]]
                atom_idx2 = atom_mapping[edge[1]]
                rdkit_mol.AddBond(atom_idx1, atom_idx2)

            # Generate a 2D depiction
            rdkit_mol.UpdatePropertyCache(strict=False)
            Chem.Draw.rdDepictor.Compute2DCoords(rdkit_mol)

            # Create a drawer with atom options
            drawer = Chem.Draw.rdMolDraw2D.MolDraw2DCairo(500, 500)
            opts = drawer.drawOptions()


            for i in range(rdkit_mol.GetNumAtoms()):
                atom = rdkit_mol.GetAtomWithIdx(i)
                opts.atomLabels[i] = atom.GetSymbol()


            drawer.DrawMolecule(rdkit_mol)
            drawer.FinishDrawing()

            # Give the image a unique name of the chemical formula, which is stored in the Molecule column in the unique_molecule_frame dataframe
            with open(f'{row["Molecule"]}.png', 'wb') as f:
                f.write(drawer.GetDrawingText())

    return id_frame

# Edit a frame in pdb format.
def pdb(frame, element_masses, molecule_id) -> Tuple[pd.DataFrame, pd.DataFrame]:

    #drop the first and last line
    frame=frame.drop(frame.index[[0,len(frame) - 1]])
    #Split the frame into columns  and label them
    split_frame = frame[0].str.split(expand=True) 
    split_frame.columns = ['Label', 'Count', 'Atom', 'Molecule', 'X', 'Y', 'Z', 'Charge']
    #now we try to convert the charge column to float
    try:
        split_frame['Charge'] = split_frame['Charge'].astype(float)
    except ValueError:
        #if it fails, we delete the charge column
        split_frame=split_frame.drop(['Charge'], axis = 1)
        #and add a new column with the value 0
        split_frame['Charge'] = 0
    #drop label and count column
    split_frame=split_frame.drop(['Label', 'Count'], axis = 1)
    #In atom column, just keep the first character
    split_frame['Atom'] = split_frame['Atom'].str[0]
    #Add a new column with the atomic mass of each atom.
    split_frame['Mass'] = split_frame['Atom'].map(element_masses)
    #Add the molecule number to the frame       
    split_frame['Molecule'] = molecule_id 

    return frame, split_frame
    
# Edit a frame in xyz format.
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

# Edit a frame in lammps format.
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

