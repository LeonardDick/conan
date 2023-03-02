#!/usr/bin/python3.8
# The program is written by Leonard Dick, 2022

import time
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from prettytable import PrettyTable
import gc
from typing import Tuple
#----own modules----
import defdict as ddict

#----Molecule identifier----
def molident(df, CNT, carbon_wall, pdb_info)  ->  Tuple[pd.DataFrame, pd.DataFrame]:
    #get the command line arguments
    args = ddict.read_commandline()
    #Time check
    MolIdent_time=time.time() 

    #XYZ FILES AND LMP FILES    
    if args["trajectoryfile"].endswith('.xyz') or args["trajectoryfile"].endswith('.lmp'):                                                         
        ddict.printLog('')
        ddict.printLog('-> Identifying all molecules in the system.')        
        ddict.printLog('-> The first frame is used for the molecule identification.')
        #make a new dataframe with only the moving atoms.
        moving_atoms=pd.DataFrame(df[df['Frozen']==False])      
        #new column in the dataframe with all values set to the index                            
        moving_atoms['Index']=moving_atoms.index                       
        #reset the indexing of the dataframe                             
        moving_atoms=moving_atoms.reset_index(drop=True)
        #save the element info of the atoms in a list.                                            
        elements=moving_atoms['Atom'].tolist()   
        #list with all the different elements in the system.                                                   
        elements_list=list(set(elements))                                                           
        ddict.printLog('-> Identified elements in the system:', elements_list)
        ddict.printLog('-> Setting up the distance matrix.')

        #DISTANCE MATRIX
        distance_time=time.time()
        #fill an array with the x, y and z positions of all atoms.
        positions=np.empty((moving_atoms.shape[0],3))
        positions[:,0]=moving_atoms['x'].values                                                     
        positions[:,1]=moving_atoms['y'].values
        positions[:,2]=moving_atoms['z'].values
        #Set up the distance matrix
        distance_matrix=cdist(positions, positions) 
        #Set the lower triangle of the matrix to 0. It is not needed.                                                
        distance_matrix=np.triu(distance_matrix)                                                    
        ddict.printLog('-> The distance matrix has been set up in %0.3f seconds.' % (time.time()-distance_time))

        #COVALENT DISTANCE MATRIX
        #time check
        covalent_time=time.time()
        #Remove the covalent distances from the distance matrix. (plus a tolerance factor)
        covalent_dist=ddict.dict_covalent()
        #Name all columns in the distance matrix dataframe with the respective elements entrie.
        for i in range(len(elements)):                                                               
            distance_matrix=pd.DataFrame(distance_matrix,columns=elements)
        #now from every 'C' column, we want to subtract the covalent distance of the 'C' element from the covalent_dict.
        for i in elements_list:
            distance_matrix[i]=distance_matrix[i]-covalent_dist[i]
        #add a new first coulumn with the elements of the atoms.
        distance_matrix.insert(0,'Atom',elements)
        #now repeat the same for the rows.
        distance_matrix=distance_matrix.set_index('Atom')
        for i in elements_list:
            distance_matrix.loc[i]=distance_matrix.loc[i]-covalent_dist[i]
        distance_matrix=distance_matrix.reset_index()
        #drop the first column again.
        distance_matrix=distance_matrix.drop(columns=['Atom'])
        #now add a tolerance of 0.5 angstrom to the distance matrix. Meaning to substract 0.5 angstrom from every value in the distance matrix dataframe.
        distance_matrix=distance_matrix-0.5
        
        #BOND ANALYSIS 
        ddict.printLog('-> Starting the bond analysis.')
        distance_matrix=np.matrix(distance_matrix)
        #now make a new column in the moving_atoms dataframe with the molecule number. If it is 0, the atom is not assigned to a molecule yet.                     
        while moving_atoms[(moving_atoms['Molecule']==0)].empty==False:
            #find the first atom with mol number 0
            first_atom=moving_atoms[(moving_atoms['Molecule']==0)].iloc[0]            
            #Assign new molecule number to the atom         
            moving_atoms.loc[first_atom.name,'Molecule']=moving_atoms['Molecule'].max()+1      
            #Add atom to the list of atoms that are already assigned to a molecule.     
            current_structure=[first_atom.name] 
            #loop over the atoms in the current molecule and add all the atoms that are connected to the investigated atom to the current molecule.
            #the loop is repeated, when new atoms are added to the current molecule.                                                    
            while len(current_structure)>0:  
                #pop the first atom from the list of atoms that are already assigned to a molecule. The pop function also removes the atom from the current_structure list simultaneously.                                                       
                current_atom=current_structure.pop(0)
                for atom in moving_atoms[(moving_atoms['Molecule']==0)].index:
                    if distance_matrix[current_atom,atom]<0:
                        current_structure.append(atom)
                        moving_atoms.loc[atom,'Molecule']=moving_atoms['Molecule'].max()
        #At last adjust the indexing of the dataframe. The original index saved in the 'index' column.
        moving_atoms=moving_atoms.set_index('Index')                                                
        moving_atoms['Molecule']=moving_atoms['Molecule']+len(CNT)+len(carbon_wall)    
        df.loc[moving_atoms.index,'Molecule']=moving_atoms['Molecule'].values  

    #PDB FILES
    if args["trajectoryfile"].endswith('.pdb'):
        #find the maximum in the Molecule column of the dataframe.
        max_molecule=df['Molecule'].max()
        #find the maximum in the Molecule column of the pdb_info dataframe.
        max_molecule_pdb=pdb_info['Molecule'].max()
        #rename the Molecule column of the pdb_info dataframe to Molecule_name.
        pdb_info=pdb_info.rename(columns={'Molecule':'Molecule_name'})
        #add the molecule_pdb molecule_name column to the dataframe. Add the max_molecule to the Molecule column of the dataframe.
        df['Molecule_name']=pdb_info['Molecule_name'].values + max_molecule
        #when the value in the Molecule column of the dataframe is 0, replace it with the value in the Molecule_name column.
        df.loc[df['Molecule']==0,'Molecule']=df.loc[df['Molecule']==0,'Molecule_name']
        #drop the Molecule_name column.
        df=df.drop(columns=['Molecule_name'])
        #remove the index column.
        df=df.reset_index(drop=True)

    #GENERAL
    #Check how many different kinds of molecules there are in the system.
    ddict.printLog('-> Checking the number of different species.')
    molecules=df['Molecule'].unique()
    molecule_dict={}   
    #For each molecule, save the atoms that belong to it in a dictionary.
    for i in range(len(molecules)):
        molecule_dict[molecules[i]]=df[(df['Molecule']==molecules[i])]['Atom'].values 
    #Make a dataframe drom the dictionary entries.                 
    molecule_df=pd.DataFrame.from_dict(molecule_dict,orient='index')                                                         
    molecule_df=molecule_df.reset_index()
    molecule_df=molecule_df.rename(columns={'index':'Molecule'})
    molecule_df=molecule_df.rename(columns={0:'Atoms'})
    #Combine every column into one column with all the atoms in the molecule.
    molecule_df['Atoms']=molecule_df[molecule_df.columns[1:]].apply(lambda x: ','.join(x.dropna().astype(str)),axis=1)
    #Drop all columns, except for 'Molecule' and 'Atoms'.       
    molecule_df=molecule_df.drop(molecule_df.columns[2:],axis=1)
    #Make a new dataframe with the unique values of the 'Atoms' column.                                                             
    molecule_df2=molecule_df['Atoms'].unique()                                                                               
    molecule_df2=pd.DataFrame(molecule_df2)
    molecule_df2['Occurrence']=0
    #Check for every row in molecule_df2 how many times the 'Atoms' column in that row appear in the 'Atoms' column of molecule_df. The result is then saved in the 'Occurrence' column.
    molecule_dict2={}
    for i in range(len(molecule_df2)):                                                                                       
        molecule_df2.loc[i,'Occurrence']=len(molecule_df[(molecule_df['Atoms']==molecule_df2.loc[i,0])])
        #Save the molecules that have the same atoms in a dictionary.
        molecule_dict2[molecule_df2.loc[i,0]]=molecule_df[(molecule_df['Atoms']==molecule_df2.loc[i,0])]['Molecule'].values
    #Add a new first column to molecule_df2 named 'Species' and fill it with the index number +1.  
    molecule_df2.insert(0,'Species',range(1,1+len(molecule_df2)))                                                            
    ddict.printLog(molecule_df2)
    ddict.printLog('')
    #Add a new column to the df which will be returned to the main loop with the species number, derived by the molecule kind.
    df['Species']=0                                                                                             
    for i in range(len(molecule_dict2)):
        df.loc[df['Molecule'].isin(molecule_dict2[molecule_df2.loc[i,0]]),'Species']=molecule_df2.loc[i,'Species']
    ddict.printLog('-> The species numbers have been assigned to the atoms in the first frame.')
    #print the frame to a file system.info
    df.to_csv('system.info',sep='\t',index=False)
    ddict.printLog('-> The first frame has been printed to the file system.info')
    ddict.printLog('')  
    #Make a table with the info from above.
    table = PrettyTable()                                                                                                    
    table.field_names = ["No. of molecules", "CNTs", "Walls", "No. of IL molecules"]
    table.add_row([len(molecules), len(CNT), len(carbon_wall), len(molecules)-len(CNT)-len(carbon_wall)])
    ddict.printLog(table)  
    ddict.printLog('')  
    ddict.printLog("-> The molecular recognition took %0.3f seconds to run" % (time.time() - MolIdent_time))
    ddict.printLog('')  

    return molecule_df2, df

