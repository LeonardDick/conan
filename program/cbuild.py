#!/usr/bin/python3.8
# The program is written by Leonard Dick, 2022

# MODULES
import shutil
import time                                                    
import os
import numpy as np
import pandas as pd
#-----own modules-----#
import defdict as ddict

#RUNTIME
cbuild_time=time.time()

#CHECK
if os.path.exists('cbuild'):
    #check if a folder called cbuild_old exists in the cbuild ffoder
    if os.path.exists('cbuild_old'):
        ddict.printLog('Removing cbuild_old folder...')
        #remove the old cbuild_old folder (if it exists)	
        shutil.rmtree('cbuild_old')                          
    ddict.printLog('Renaming existing cbuild folder to cbuild_old...')
    os.rename('cbuild', 'cbuild_old')
os.mkdir('cbuild')
#finally move the cbuild_old folder to the cbuild folder
if os.path.exists('cbuild_old'):
    ddict.printLog('Moving cbuild_old folder to cbuild folder...')
    shutil.move('cbuild_old', 'cbuild')
ddict.printLog('')
        
#PARAMETER SECTION
anginbohr=1.889726
carbon_distance=1.42
interlayer_distance=3.35
hex_d=carbon_distance*(3**0.5)
distance_row_y=hex_d/2

#MAIN
ddict.printLog_red('CBUILD mode')
ddict.printLog('')
#carbon atom name
carbon='C'
#general questions for the user what structures the program shall produce                                                       
choice_length=float(ddict.get_input('Write coordinate file in [1]angstrom or [2]bohr?   '))
if choice_length==2:
    carbon_distance=carbon_distance*anginbohr
    distance_row_y=distance_row_y*anginbohr
    interlayer_distance=interlayer_distance*anginbohr
    ddict.printLog('--> Multiplying values by %0.5f to get from Angstrom to bohr.'%(anginbohr))
    ddict.printLog('')
else:
    ddict.printLog('')
ddict.printLog('Parameters:')
ddict.printLog('C-C distance: %0.3f'%(carbon_distance))
ddict.printLog('Interlayer distance: %0.3f'%(interlayer_distance))
change_parameter=ddict.get_input('Do you want to change parameters? [y/n]   ')
if change_parameter=='y':
    carbon_distance=float(ddict.get_input('C-C distance [angstrom]   '))
    interlayer_distance=float(ddict.get_input('Interlayer distance [angstrom]   '))
ddict.printLog('')
iflayer=ddict.get_input('Do you want to generate a pore structure? [y/n]   ')
if iflayer=='y':
    iflayer2='y'
else:
    iflayer2=ddict.get_input('Do you want to generate carbon walls? [y/n]   ')


#CARBON WALL SECTION       
if iflayer2=='y':
    size_graphsheet=float(ddict.get_input('How large should the carbon walls be? [angstrom]   '))
    num_layer=float(ddict.get_input('How many carbon wall layers?   '))
    ddict.printLog('')
    if choice_length==2:
        size_graphsheet=size_graphsheet*anginbohr

    #CARBON WALL SECTION
    #set all variables to zero
    atomcount=0
    tmp_CW_xyz=[]
    coord_x=-carbon_distance*1.5
    coord_y=0
    coord_z=0
    maximumx=0
    maximumy=0
    maximumz=0
    counterx=0
    countery=0
    #the run variables are used to account for shifts in stacked layers etc.
    runx=2
    runy=1
    runz=1
    #start the loop
    while True:
        while True:
            if (runy%2)==0:
                runx=2
            else:
                runx=3
            while True:
                coord_x=coord_x+carbon_distance
                if maximumx<=coord_x:
                    maximumx=coord_x
                if runx==3:
                    runx=1
                    continue
                if (runz%2)==0:
                    #append coordinates to list
                    tmp_CW_xyz.append([carbon, coord_x, coord_y, coord_z]) 
                else:
                    xshift=coord_x+carbon_distance
                    #append coordinates to list
                    tmp_CW_xyz.append([carbon, xshift, coord_y, coord_z])  
                runx=runx+1
                atomcount=atomcount+1
                counterx+=1
                if coord_x>=size_graphsheet:
                    #assure that PBC is not violated
                    if (counterx%2)==0:
                        break
            if (runy%2)==0:
                coord_x=-carbon_distance*1.5
            else:
                coord_x=-carbon_distance
            runy=runy+1
            countery+=1
            if coord_y>=size_graphsheet:
                #assure that PBC is not violated
                if (countery%2)==0:
                    break
            coord_y=coord_y+distance_row_y
            if maximumy<=coord_y:
               maximumy=coord_y
        coord_y=0
        if runz==num_layer:
            break
        #start new layer
        coord_z=coord_z+interlayer_distance
        if maximumz<=coord_z:
            maximumz=coord_z
        runz=runz+1
    
    #make list to dataframe and write to file
    df_CW=pd.DataFrame(tmp_CW_xyz, columns=['atom', 'x', 'y', 'z'])    
    with open ('cbuild/carbon_wall.xyz', 'w') as f:      
        f.write('%d\n\n' % (atomcount))           
    df_CW.to_csv('cbuild/carbon_wall.xyz', sep='\t', header=False, index=False, mode='a') 
    
    #print some additional information
    ddict.printLog('Carbon Wall:')            
    if num_layer==1:
        ddict.printLog('Sheet size in x direction: %0.3f'%(maximumx))
    else:
        ddict.printLog('Sheet size in x direction: %0.3f'%(maximumx+carbon_distance))
    ddict.printLog('Sheet size in y direction: %0.3f'%(maximumy))
    ddict.printLog('Sheet size in z direction: %0.3f'%(maximumz))
    if num_layer==1:
        ddict.printLog('Boxsize for PBC in x direction: %0.3f'%(maximumx+carbon_distance))
    else:
        ddict.printLog('Boxsize for PBC in x direction: %0.3f'%(maximumx+carbon_distance*2))    
    ddict.printLog('Boxsize for PBC in y direction: %0.3f'%(maximumy+distance_row_y))
    #print number of atoms in the graph sheet                 
    ddict.printLog('Number of graph sheet atoms: {0}'.format(len(df_CW)))
    ddict.printLog('')

#CNT SECTION
#general questions regarding the CNT buildup
if iflayer=='y':
    tube='y'
else:
    if iflayer2=='y':
        tube=ddict.get_input('Do you also want to produce a CNT? [y/n]   ')
    else:
        tube=ddict.get_input('Do you want to produce a CNT? [y/n]   ')
if tube=='y':
    which_conf=2
    ddict.printLog('--> CNT is produced zigzag conformation.')
    atom_circ=float(ddict.get_input('How many atoms for the circumference of CNT?   '))
    tubelength=float(ddict.get_input('How long do you want the tube to be? [angstrom]   '))
    ddict.printLog('--> Adjusting tube length to assure periodicity in z direction.')
    ddict.printLog('\nCNT:')
    #set tube length to 0 initially
    tubel_z=0

    # create xyz list and set angle stepsize
    step=360/atom_circ 
    tmp_tube_xyz=[]

    #generate CNT in zigzag conformation
    if which_conf==2:
        #get circcumference and diameter
        circumference=atom_circ*hex_d
        diameter=circumference/3.14159265
        ddict.printLog('Circumference: %0.3f'%(circumference))
        ddict.printLog('Diameter: %0.3f'%(diameter))
        distx=(np.sin(np.deg2rad(0))*diameter/2-np.sin(np.deg2rad(step/2))*diameter/2)
        disty=(np.cos(np.deg2rad(0))*diameter/2-np.cos(np.deg2rad(step/2))*diameter/2)
        zstep=(carbon_distance**2-distx**2-disty**2)**0.5
        #calculate the length of the tube
        tubel_z=carbon_distance+zstep
        while True:
            if tubel_z>=tubelength:
                break
            tubel_z=tubel_z+carbon_distance*2+zstep*2
        ddict.printLog('Tubelength: %0.3f'%(tubel_z+zstep))
        PBC_z=tubel_z+carbon_distance+zstep
        ddict.printLog('Boxsize for PBC in z direction: %0.3f'%(PBC_z))
    #define the xy-center of the tube (depending if just the tube or pore structure is produced)
    if iflayer2=='y':
        centerx=maximumx/2
        centery=maximumy/2 
    else:
        centerx=0
        centery=0
    #set all variables to 0
    centerz=0
    count_deg=0
    tubez=0
    count_stepz=1
    runz=2
    stopz=tubel_z+carbon_distance*1.5
    radius=diameter/2
    while True: 
        if runz<=2:
            count_deg=0
        else:
            count_deg=step/2
        while True:
            tubex=centerx+np.sin(np.deg2rad(count_deg))*radius
            tubey=centery+np.cos(np.deg2rad(count_deg))*radius
            #append atom coordinates to list
            tmp_tube_xyz.append([carbon, tubex, tubey, tubez])
            count_deg=count_deg+step
            if count_deg>=360:
                break        
        count_stepz=count_stepz+1
        if  count_stepz==3:
            tubez=tubez+carbon_distance
            count_stepz=1
        else:
            tubez=tubez+zstep 
        runz=runz+1
        if runz>4:
            runz=1
        if tubez>=stopz:
            break
    
    #make a dataframe from the list and write it to a file
    df_tube=pd.DataFrame(tmp_tube_xyz, columns=['atom', 'x', 'y', 'z']) 
    with open ('cbuild/tube.xyz', 'w') as f:     
        #write number of atoms in the tube to the first line
        f.write('%d\n\n' % (len(df_tube)))
    df_tube.to_csv('cbuild/tube.xyz', sep='\t', header=False, index=False, mode='a')
    #print number of atoms in the tube
    ddict.printLog('Number of tube atoms: {0}'.format(len(df_tube)))


#Pore section
if iflayer=='y':
    ddict.printLog('\nPore:')
    ddict.printLog('--> Combining layers and tube.')
    #remove all improper atoms from the wall layers
    #first define a new list
    tmp_pore_xyz=[]
    #ask user which structure shall be produced
    tubekind=float(ddict.get_input('[1]Pore with walls on both ends or [2]Closed Pore with 1 wall (left and right)?   '))  
    #define an additional column named 'distance' in the dataframe. If the distance is smaller than the radius of the tube+carbon_distance, the atom row is deleted
    df_CW['distance']=np.sqrt((df_CW['x']-centerx)**2+(df_CW['y']-centery)**2)    
    df_CW_dummy=df_CW[df_CW['distance']>radius+carbon_distance] 
    #delete the column 'distance' again
    del df_CW_dummy['distance']        
    #reset the index of the dataframe. Ignore_index=True is important, otherwise the index of the tube atoms is not reset                                              
    df_CW_dummy=df_CW_dummy.reset_index(drop=True)                                  
    df_pore=pd.concat([df_CW_dummy, df_tube], ignore_index=True) 
    #Find the maximum z value in the pore dataframe
    maximumz_pore=df_pore['z'].max()                                                
    
    if tubekind==1: 
        #here a second wall is added to the other side of the tube. Again, the same procedure as above is used
        #create dummy df_CW_shift dataframe
        df_CW_shift=df_CW_dummy
        #match the maximum value of the carbon walls to the maximum value of the pore
        df_CW_shift['z']=df_CW_shift['z']+maximumz_pore-maximumz     
        #Add the df_CW_shift to the df_pore dataframe           
        df_pore=pd.concat([df_pore, df_CW_shift], ignore_index=True)     
        #write the pore structure to a file, with the number of atoms in the first line        
        with open('cbuild/pore.xyz', 'w') as f: 
            f.write('%d\n\n' % (len(df_pore)))
        df_pore.to_csv('cbuild/pore.xyz', sep='\t', header=False, index=False, mode='a')
    
    #create a closed pore structure with 1 wall, left and right side
    #the approach is similar to the one above. To close the tube on one side, all atoms of the wall are deleted which are not within the tube
    if tubekind==2:
        df_CW_dummy2=df_CW[df_CW['distance']<radius-carbon_distance]    
        del df_CW_dummy2['distance']
        df_CW_dummy2=df_CW_dummy2.reset_index(drop=True)
        df_CW_dummy2 ['z']=df_CW_dummy2['z']+maximumz_pore-maximumz
        df_CW_dummy2=df_CW_dummy2[df_CW_dummy2['z']==maximumz_pore] 
        #the dastaframes are concatenated to one dataframe
        df_pore_left=pd.concat([df_pore, df_CW_dummy2], ignore_index=True)
        with open('cbuild/pore_left.xyz', 'w') as f:                                        
            f.write('%d\n\n' % (len(df_pore_left)))
        df_pore_left.to_csv('cbuild/pore_left.xyz', sep='\t', header=False, index=False, mode='a')
        #mirror the dataframe to produce the 'right' carbon structure. The z values are mirrored and the lenght of the tube is added to the z values.
        df_pore_right=df_pore_left
        df_pore_right['z']=df_pore_right['z']*(-1)+maximumz_pore
        with open('cbuild/pore_right.xyz', 'w') as f:                                        
            f.write('%d\n\n' % (len(df_pore_right)))
        df_pore_right.to_csv('cbuild/pore_right.xyz', sep='\t', header=False, index=False, mode='a')
ddict.printLog("\nCBuild mode finished in %0.3f seconds" % (time.time()-cbuild_time))
