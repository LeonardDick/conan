# The program is written by Leonard Dick, 2024

# MODULES
import pandas as pd
import numpy as np
# ----- Own modules ----- #
import defdict as ddict

# ARGUMENTS
args = ddict.read_commandline()



def COM_calculation(split_frame):

    # We now calculate the center of mass (COM) for each molecule
    # Convert all values to float (this is needed so that the agg-function works)
    split_frame['X'] = split_frame['X'].astype(float)
    split_frame['Y'] = split_frame['Y'].astype(float)
    split_frame['Z'] = split_frame['Z'].astype(float)
    split_frame['Mass'] = split_frame['Mass'].astype(float)

    # Precompute total mass for each molecule
    total_mass_per_molecule = split_frame.groupby('Molecule')['Mass'].transform('sum')

    # Calculate mass weighted coordinates
    split_frame['X_COM'] = (split_frame['X'] * split_frame['Mass']) / total_mass_per_molecule
    split_frame['Y_COM'] = (split_frame['Y'] * split_frame['Mass']) / total_mass_per_molecule
    split_frame['Z_COM'] = (split_frame['Z'] * split_frame['Mass']) / total_mass_per_molecule

    # Calculate the center of mass for each molecule
    mol_com = split_frame.groupby('Molecule').agg(
        Species=('Species', 'first'),
        X_COM=('X_COM', 'sum'),
        Y_COM=('Y_COM', 'sum'),
        Z_COM=('Z_COM', 'sum')
    ).reset_index()

    return mol_com, split_frame

# print com positions together with the pore atoms in the same frame and into a trajectory output file
def print_com_pore(mol_com, split_frame, output_file):
    # Print the COM positions together with the pore atoms in the same frame and into a trajectory output file
    with open(output_file, 'a') as f:
        f.write('ITEM: TIMESTEP\n')
        f.write('0\n')
        f.write('ITEM: NUMBER OF ATOMS\n')
        f.write(str(len(split_frame) + len(mol_com)) + '\n')
        f.write('ITEM: BOX BOUNDS pp pp pp\n')
        f.write('0 1\n')
        f.write('0 1\n')
        f.write('0 1\n')
        f.write('ITEM: ATOMS id type x y z\n')

        # Print the pore atoms
        for i in range(len(split_frame)):
            f.write(str(split_frame['ID'].iloc[i]) + ' ' + str(split_frame['Type'].iloc[i]) + ' ' + str(
                split_frame['X'].iloc[i]) + ' ' + str(split_frame['Y'].iloc[i]) + ' ' + str(split_frame['Z'].iloc[i]) + '\n')

        # Print the COM positions
        for i in range(len(mol_com)):
            f.write(str(i + 1 + len(split_frame)) + ' ' + str(mol_com['Species'].iloc[i]) + ' ' + str(mol_com['X_COM'].iloc[i]) + ' ' + str(
                mol_com['Y_COM'].iloc[i]) + ' ' + str(mol_com['Z_COM'].iloc[i]) + '\n')

    return