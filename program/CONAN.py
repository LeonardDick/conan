#!/usr/bin/python3.8
# The program is written by Leonard Dick, 2022

# MODULES
import time                                                        
import os                                                        
import argparse
import sys
#-----own modules-----#
import defdict as ddict

#RUNTIME
start_time = time.time()

#MAIN
#rename old log file
if os.path.exists('conan.log'):                                 
    if os.path.exists('conan.log.old'):
        os.remove('conan.log.old')
    os.rename('conan.log', 'conan.log.old')

#LOGO
ddict.printLog_yellow('###########################################')
ddict.printLog_yellow('##                                       ##')
ddict.printLog_yellow('##   #####  #####  #   #  #####  #   #   ##')
ddict.printLog_yellow('##   #      #   #  ##  #  #   #  ##  #   ##')
ddict.printLog_yellow('##   #      #   #  # # #  #####  # # #   ##')
ddict.printLog_yellow('##   #      #   #  #  ##  #   #  #  ##   ##')
ddict.printLog_yellow('##   #####  #####  #   #  #   #  #   #   ##')
ddict.printLog_yellow('##                                       ##')
ddict.printLog_yellow('###########################################')
ddict.printLog('')

#INFO
#refer to the documentation for more information on the program. Website is conan.readthedocs.io
ddict.printLog('Find the documentation on the CONAn website: http://con-an.rtfd.io')
#ddict.printLog('If you use CONAn in your research, please cite the following paper:')
#ddict.printLog('XXX')
ddict.printLog('')


#ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--trajectoryfile",  help="Load the trajectory file. Either in .xyz, .pdb or .lmp format.")
parser.add_argument("-c", "--cbuild", action='store_true', help="Generate carbon structures")
parser.add_argument("-b", "--box", action='store_true', help="Build a simulation box" )
parser.parse_args()
args = vars(parser.parse_args())
#If no arguments are given in the commnad line, print help and exit.
if len(sys.argv)==1:
    parser.print_help()
    ddict.printLog('')
    ddict.printLog('No arguments given. Exiting...')
    sys.exit(1)

#CBUILD SECTION
if args['cbuild']:
    import cbuild

#SIMULATION SETUP SECTION
if args['box']:
    import simbox

#TRAJECTORY ANALYSIS SECTION
if args['trajectoryfile']:
    import traj_info
    # Running the trajectory_info module to open the trajectory file and read the first frame.
    number_of_atoms, simbox_x, simbox_y, simbox_z, first_frame, pdb_info=traj_info.trajfile_read()
    # Running the CNT_identification module to identify the CNTs and carbon walls in the first frame, returning various values.
    tube_atoms, tuberadius, tubecenter_x, tubecenter_y, tubecenter_z, carbon_wall_z, length_CNT, max_z_CNT, min_z_CNT, CNTs, carbon_walls, CNT_atoms=traj_info.CNT_identification(first_frame)
    #Now running the analysis_opt module to analyse the trajectory.
    import traj_an
    traj_an.analysis_opt(first_frame, simbox_x, simbox_y, simbox_z, number_of_atoms, tube_atoms, tuberadius, tubecenter_x, tubecenter_y, tubecenter_z, carbon_wall_z, length_CNT, max_z_CNT, min_z_CNT, CNTs, carbon_walls, CNT_atoms, pdb_info)

ddict.printLog('The program took %0.3f seconds to run.' % (time.time() - start_time))