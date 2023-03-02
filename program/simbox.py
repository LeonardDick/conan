#!/usr/bin/python3.8
# The program is written by Leonard Dick, 2022

'''
This program shall add carbon structures and liquid bulk xyz files to one simulation box.
At first the program shall read the xyz files and store the coordinates in a dataframe.
Then the program shall add the carbon structures to the bulk liquid.
Finally the program shall write the simulation box to a xyz file.
First ask the user for the general setup of the simulation box.
P stands for pore, B for liquid bulk, and W for carbon wall.
If pore_left and pore_right from the cbuild section are used, the user shall enter L for the left pore and R for the right pore.
by giving a code, like BWBP, the user can define the setup of the simulation box.
'''

# MODULES
import time
import sys
import os
import pandas as pd
#---own modules---
import defdict as ddict

#RUNTIME
simbox_time = time.time()

#MAIN
ddict.printLog('')
ddict.printLog_red('SIMBOX mode')
ddict.printLog('')
ddict.printLog('This program shall add carbon structures and liquid bulk xyz files to one simulation box.')
ddict.printLog('P stands for pore, B for liquid bulk, and W for carbon wall.')
ddict.printLog('If pore_left and pore_right from the cbuild section are used, the user shall enter L for the left pore and R for the right pore.')
combination=ddict.get_input('Please enter the wanted combination for the simulation box [eg.: BPBW]: ')
ddict.printLog('')
#Now we split the string into a list.
combination_list=list(combination)
ddict.printLog('The combination list is: %s' % (combination_list))
#Make another list with just the unique entries
combination_list_unique=list(set(combination_list))
#write the possible combination list
possible_letters=['B', 'P', 'W', 'R', 'L']
#If a different letter than those is goven, the program shall exit.
for i in combination_list:
    if i not in possible_letters:
        ddict.printLog('The combination list contains a wrong letter. Exiting...')
        sys.exit(1)

#Now we have to find the necessary files.
#The liquid bulk file shall be named bulk.xyz.
#The carbon wall file shall be named carbon_wall.xyz.
#The pore file shall be named pore.xyz.
#The pore_left file shall be named pore_left.xyz.
#The pore_right file shall be named pore_right.xyz.
#All files are either located in the current directory or in the cbuild directory.
#If the file is not found in the current directory, the program shall search in the cbuild directory.
#If the file is not found in the cbuild directory, the program shall exit.
#instead of wtiting the same code for each file, loop over the combination_list.
#first create list with the possible file names
file_name_list=dict([('B', 'bulk'), ('P', 'pore'), ('W', 'carbon_wall'), ('R', 'pore_right'), ('L', 'pore_left')])
ddict.printLog('The file name list is: %s' % (file_name_list))
for i in possible_letters:
    if i in combination_list_unique:
        #create the file name
        file_name='%s.xyz' % (file_name_list[i])
        try:
            #read file to dataframe
            df=pd.read_csv(file_name, sep='\s+', header=None, skiprows=2, names=['atom', 'x', 'y', 'z'])
        except:
            try:
                df=pd.read_csv('cbuild/%s' % (file_name), sep='\s+', header=None, skiprows=2, names=['atom', 'x', 'y', 'z'])
            except:
                ddict.printLog('The %s file could not be found. Exiting...' % (i))
                sys.exit(1)
        #rename the dataframe variable to the file name
        exec('%s=df' % (file_name_list[i]))
        #print the dataframe
        ddict.printLog('The %s file was found.' % (file_name))
ddict.printLog('')

#Now we need to find the minimal z value in all dataframes.
#The minimal z value shall be stored in a variable and used to shift the pore and the carbon wall.
#Find the minimal z value in the pore dataframe. If it is not zero, shift the pore dataframe to 0.
#As for the first step, the code can be shortened by using a loop over all file names.
for i in possible_letters:
    if i in combination_list_unique:
        #create the file name
        file_name='%s.xyz' % (file_name_list[i])
        #find the minimal z value
        exec('%s_min_z=%s[\'z\'].min()' % (file_name_list[i], file_name_list[i]))
        #shift the dataframe to 0
        if eval('%s_min_z' % (file_name_list[i])) != 0:
            exec('%s[\'z\']=%s[\'z\']-%s_min_z' % (file_name_list[i], file_name_list[i], file_name_list[i]))
        #find the maximal z value
        exec('%s_max_z=%s[\'z\'].max()' % (file_name_list[i], file_name_list[i]))

#Now start building the simulation box by setting up an empty dataframe with the correct column names.
simbox=pd.DataFrame(columns=['atom', 'x', 'y', 'z']) 
tmp = []
simbox_max_z=0
#like the last two sections, a loop can be used to shorten the code.
for i in combination_list:
    #create the file name
    name='%s' % (file_name_list[i])
    #make a dummy dataframe
    exec('%s_dummy=%s.copy()' % (file_name_list[i], file_name_list[i]))
    #shift the dataframe in z direction by simbox_max_z
    exec('%s_dummy[\'z\']=%s_dummy[\'z\']+simbox_max_z' % (file_name_list[i], file_name_list[i]))
    #append the dataframe to the tmp list
    exec('tmp.append(%s_dummy)' % (file_name_list[i]))
    simbox_max_z = max([df['z'].max() for df in tmp])+3         #Find the maximal z in the tmp list
#Concatenate the tmp list to the simbox dataframe. Use concat and not append, because append does not work with the ignore_index argument.
simbox = pd.concat(tmp, ignore_index=True)                     
ddict.printLog('The simbox dataframe was created.')
ddict.printLog(simbox)

#file creation
#if the simbox.xyz file exists, the program shall rename it to simbox.xyz.old.
if os.path.isfile('simbox.xyz'):
    if os.path.isfile('simbox.xyz.old'):
        os.remove('simbox.xyz.old')
        ddict.printLog('The simbox.xyz.old file was removed.')
    os.rename('simbox.xyz', 'simbox.xyz.old')
    ddict.printLog('The simbox.xyz file was renamed to simbox.xyz.old.')
with open('simbox.xyz', 'w') as f:                              #Write the simbox.xyz file.
    f.write('%d\n' % (len(simbox)))
    f.write('\n')
simbox.to_csv('simbox.xyz', sep='\t', header=False, index=False, mode='a')
ddict.printLog('The simbox module took %s seconds to run.' % (time.time() - simbox_time))