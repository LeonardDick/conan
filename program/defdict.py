#!/usr/bin/python3.8
# The program is written by Leonard Dick, 2022

''' This file contains all globally needed functions and dictionaries.'''

#DEFINITIONS
#define red color
def prRed(skk): print("\033[91m {}\033[00m" .format(skk))       
#define yellow color
def prYellow(skk): print("\033[93m {}\033[00m" .format(skk)) 

#print to terminal and log file using printLog
def printLog(*args, **kwargs):
    print(*args, **kwargs)
    with open('conan.log','a') as file:
        print(*args, **kwargs, file=file)
#printLog with color
def printLog_red(*args, **kwargs):                            
    prRed(*args, **kwargs)
    with open('conan.log','a') as file:
        print(*args, **kwargs, file=file)
def printLog_yellow(*args, **kwargs):                            
    prYellow(*args, **kwargs)
    with open('conan.log','a') as file:
        print(*args, **kwargs, file=file)

#input with log file 
def get_input(*args)-> str:
    var=input(*args)
    with open('conan.log','a') as file:
        print(*args, var, file=file)
    return var

#read command line arguments
def read_commandline()-> dict:
    import argparse
    #ARGUMENTS
    parser = argparse.ArgumentParser(description='CONAn - CONfinement Analysis')
    parser.add_argument("-f", "--trajectoryfile",  help="The xyz file containing the trajectory")
    parser.add_argument("-c", "--cbuild", action='store_true', help="Generate carbon structures")
    parser.add_argument("-b", "--box", action='store_true', help="Build a simulation box" )
    parser.parse_args()
    args = vars(parser.parse_args())
    return args

#DICTIONARIES 
#Atomic masses
def dict_mass() -> dict:
    elem_masses = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'P': 30.974, 'S': 32.065}   
    return elem_masses	
#Atomic van der Waals radii
def dict_vdW()-> dict:
    elem_vdW = {'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'P': 1.80, 'S': 1.80}
    return elem_vdW
#Atomic covalent radii
def dict_covalent()-> dict:
    elem_covalent = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'P': 1.07, 'S': 1.05}
    return elem_covalent
#Cutoff distances for molecule identification   
def dict_cutoff()-> dict:
    #set up a dictionary with cutoff distances for molecule identification by combining the covalent radii of all elemnt combinations
    comb_cutoff = dict()
    #for this we need the covalent radii
    #use defined dictionary
    elem_covalent = dict_covalent()
    #define a list of all elements in elem_covalent
    elem_list = list(elem_covalent.keys())
    #define a list of all combinations of elements in elem_list
    elem_comb = [(elem_list[i], elem_list[j]) for i in range(len(elem_list)) for j in range(i, len(elem_list))]
    #now we can define the cutoff distances
    dist_comb = [elem_covalent[elem_comb[i][0]] + elem_covalent[elem_comb[i][1]] for i in range(len(elem_comb))]
    comb_cutoff = dict(zip(dist_comb, elem_comb))
    #add a tolerance of 0.5 to the cutoff distances
    comb_cutoff = {k+0.6: v for k, v in comb_cutoff.items()}
    return comb_cutoff
