# The program is written by Leonard Dick, 2023

# MODULES
import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
# ----- Own modules ----- #
import defdict as ddict


# DATA PROCESSING
# Radial density function.
def raddens_func(df_raddens, bin_edges, num_frames, CNT_length, radius_tube, analysis_choice, args) -> None:
    # Ceck the analysis choice. -> If mass or charge density is plotted.
    if analysis_choice == '1':
        choice='Mass'
    elif analysis_choice == '2':
        choice = 'Charge'
    ddict.printLog('')
    results_rd_df = pd.DataFrame(df_raddens.iloc[:,1:].sum(axis = 0) / num_frames)
    results_rd_df.columns = [choice]

    # Add the bin edges to the results_rd_df dataframe.
    results_rd_df['Bin_lowedge'] = bin_edges[:-1]
    results_rd_df['Bin_highedge'] = bin_edges[1:]
    
    # The center of the bin is the average of the bin edges.
    results_rd_df['Bin_center'] = (bin_edges[1:]+bin_edges[:-1]) / 2                
    
    # Calculate the Volume of each bin. By setting the length of the CNT as length of a cylinder, and the radius of the bin as the radius of the cylinder.
    # Subtract the volume of the smaller cylinder from the volume of the larger cylinder (next larger bin). The volume of a cylinder is pi*r^2*h.
    vol_increment=math.pi * (bin_edges[1:] ** 2 - bin_edges[:-1] ** 2) * CNT_length     
    results_rd_df['Volume'] = vol_increment
    
    if choice == 'Mass':
        # Calculate the density of each bin by dividing the average mass by the volume.                                             
        results_rd_df['Density [u/Ang^3]'] = results_rd_df[choice] / results_rd_df['Volume']
        # Calculate the density in g/cm^3.   
        results_rd_df['Density [g/cm^3]'] = results_rd_df['Density [u/Ang^3]'] * 1.66053907  

    if choice == 'Charge':
        # Calculate the charge density in e/Ang^3.
        results_rd_df['Charge density [e/Ang^3]'] = results_rd_df[choice] / results_rd_df['Volume'] 

    # Reset the index of the dataframe.
    results_rd_df.reset_index(drop = True, inplace = True)   
    # Add a new first column with the index+1 of the bin.
    results_rd_df.insert(0, 'Bin', results_rd_df.index + 1)  

    # Plot the data.
    raddens_data_preparation = ddict.get_input('Do you want to plot the data? (y/n) ', args)
    # Adjusting the resulting dataframe for plotting, as the user set it.
    if raddens_data_preparation == 'y':
        results_rd_df_copy = results_rd_df.copy()

        # Normalization of the data.
        normalize=ddict.get_input('Do you want to normalize the increments with respect to the CNTs\' radius? (y/n) ', args)
        if normalize == 'y':
            results_rd_df['Bin_center'] = results_rd_df['Bin_center'] / radius_tube

        # Mirroring the data.
        mirror=ddict.get_input('Do you want to mirror the plot? (y/n) ', args)
        if mirror == 'y':
            # Mirror the data by multiplying the bin center by -1. Then sort the dataframe by the bin center values and combine the dataframes.
            results_rd_dummy = results_rd_df.copy()
            results_rd_dummy['Bin_center'] = results_rd_df['Bin_center'] * (-1)               
            results_rd_dummy.sort_values(by = ['Bin_center'], inplace = True)
            results_rd_df = pd.concat([results_rd_dummy, results_rd_df], ignore_index=True)                    

        # Generate the plot.
        fig, ax = plt.subplots()
        if choice == 'Mass':
            ax.plot(results_rd_df['Bin_center'], results_rd_df['Density [g/cm^3]'], '-' , label = 'Radial density function', color = 'black')
            ax.set(xlabel = 'Distance from tube center [Ang]', ylabel = 'Density [g/cm^3]', title = 'Radial density function')
            
        if choice == 'Charge':
            ax.plot(results_rd_df['Bin_center'], results_rd_df['Charge density [e/Ang^3]'], '-' , label = 'Radial density function', color = 'black')
            ax.set(xlabel = 'Distance from tube center [Ang]', ylabel = 'Charge density [e/Ang^3]', title = 'Radial density function')

        ax.grid()  
        fig.savefig("Radial_density_function.pdf")
        ddict.printLog('-> Radial density function saved as Radial_density_function.pdf\n')

        # Save the data.
        results_rd_df.to_csv('Radial_density_function.csv', sep = ';', index = False, header = True, float_format = '%.5f')

        # Radial contour plot.
        results_rd_df = results_rd_df_copy.copy()
        radial_plot = ddict.get_input('Do you also want to create a radial countour plot? (y/n) ', args)
        if radial_plot == 'y':
            theta = np.linspace(0, 2 * np.pi, 500)
            if normalize == 'n':
                r = np.linspace(0, radius_tube, len(results_rd_df['Bin_center']))
            else:
                r = np.linspace(0, 1, len(results_rd_df['Bin_center']))
                
            Theta, R = np.meshgrid(theta, r)

            if choice == 'Mass':
                values = np.tile(results_rd_df['Density [g/cm^3]'].values,
                                (len(theta), 1)).T
            elif choice == 'Charge':
                values = np.tile(results_rd_df['Charge density [e/Ang^3]'].values,
                                (len(theta), 1)).T

            fig = plt.figure()
            ax = fig.add_subplot(projection='polar')  
            c = ax.contourf(Theta, R, values, cmap='Reds')  

            ax.spines['polar'].set_visible(False)  

            # Remove angle labels.
            ax.set_xticklabels([]) 

            # Set radial gridlines and their labels.
            r_ticks = np.linspace(0, radius_tube if normalize == 'n' else 1, 5)
            ax.set_yticks(r_ticks)  
            ax.set_yticklabels(['{:.2f}'.format(x) for x in r_ticks])

            # Set the position of radial labels.
            ax.set_rlabel_position(22.5)  
            ax.grid(color = 'black', linestyle = '--', alpha = 0.5)  
            
            # Set title and add a colorbar.
            plt.title('Radial Density Contour Plot', fontsize = 20, pad = 20)
            cbar = fig.colorbar(c, ax = ax, pad = 0.10, fraction = 0.046, orientation = "horizontal")

            if choice == 'Mass':
                cbar.set_ticklabels(['{:.2f}'.format(x) for x in cbar.get_ticks()])
                cbar.set_label(r'Mass density $[g/cm^{3}]$', fontsize = 15)

            elif choice == 'Charge':
                cbar.set_ticklabels(['{:.3f}'.format(x) for x in cbar.get_ticks()])
                cbar.set_label(r'Charge density $[e/Ang^{3}]$', fontsize = 15)

            # Set the y label.
            if normalize == 'n':
                ax.set_ylabel(r'$d_{rad}$', labelpad = 10, fontsize = 20)  
            else:
                ax.set_ylabel(r'$d_{rad}$/$r_{CNT}$', labelpad = 10, fontsize = 20)

            # Save the data.
            fig.savefig("Radial_density_function_polar.pdf")
            ddict.printLog('-> Radial density function countour plot saved as Radial_density_function_polar.pdf\n')


# Accessible volume.
def acc_vol_func(id_frame, row_maxdisp_atom, CNT_length, args) -> None:

    # Write the maximal displaced atom to the log file and the terminal.
    ddict.printLog('\nMaximal displaced atom: %s (%0.3f|%0.3f|%0.3f)' % (row_maxdisp_atom['Atom'], float(row_maxdisp_atom['X']), float(row_maxdisp_atom['Y']), float(row_maxdisp_atom['Z'])))
    accessible_radius = row_maxdisp_atom['Distance']       
    ddict.printLog('Maximal displacement: %0.3f' % accessible_radius)
    ddict.printLog('Accessible volume: %0.3f' % (math.pi * accessible_radius ** 2 * CNT_length))
    ddict.printLog('')   

    # Print the carbon structure with the most displaced atom, if the user wants to.   
    pore_disp_atom=ddict.get_input('Do you want to produce a xyz file with the pore including the most displaced atom? [y/n] ', args)
    if pore_disp_atom=='y':
        # The CNT atoms are all atoms which are not NaN.
        CNT_atoms=id_frame[id_frame['CNT'].notnull()]
        f=open('pore_disp_atom.xyz', 'w')
        f.write('%d\n#Pore with most displaced atom\n' % (len(CNT_atoms) + 1))
        for index, row in CNT_atoms.iterrows():
            f.write('%s\t%0.3f\t%0.3f\t%0.3f\n' % (row['Element'], row['x'], row['y'], row['z']))
        f.write('%s\t%0.3f\t%0.3f\t%0.3f\n' % (row_maxdisp_atom['Atom'], float(row_maxdisp_atom['X']), float(row_maxdisp_atom['Y']), float(row_maxdisp_atom['Z'])))
        f.close()
        ddict.printLog('Pore with most displaced atom saved as pore_disp_atom.xyz')


# Z density function.
def zdens_func(df_zdens, bin_edges, num_frames, row_maxdisp_atom, CNT_length, radius_tube, num_incr, CNT_minz, CNT_maxz, simbox_x, simbox_y, center_tube_z, args) -> None:

    results_zd_df = pd.DataFrame(df_zdens.iloc[:,1:].sum(axis = 0) / num_frames)
    results_zd_df.columns = ['Mass']
    results_zd_df['Bin_lowedge'] = bin_edges[:-1]
    results_zd_df['Bin_highedge'] = bin_edges[1:]

    # The center of the bin is the average of its bin edges.
    results_zd_df['Bin_center'] = (bin_edges[1:] + bin_edges[:-1]) / 2                
    ddict.printLog('\nMaximal displaced atom: %s (%0.3f|%0.3f|%0.3f)' % (row_maxdisp_atom['Atom'], float(row_maxdisp_atom['X']), float(row_maxdisp_atom['Y']), float(row_maxdisp_atom['Z'])))
    accessible_radius = row_maxdisp_atom['Distance']       
    ddict.printLog('Maximal displacement: %0.3f' % accessible_radius)
    ddict.printLog('Accessible volume: %0.3f' % (math.pi * accessible_radius ** 2 * CNT_length))
    which_radius = ddict.get_input('Do you want to use the accessible radius (1) or the CNT radius (2) to compute the increments\' volume? [1/2] ', args)
    if which_radius == '1':
        used_radius = accessible_radius
    if which_radius == '2':
        used_radius = radius_tube

    # Calculate the volume of each bin. For the bulk sections, the volume is defined by the x and y values of the simbox, multiplied by the bin width. For the tube sections, the volume is defined by the pi*r^2*bin width.
    vol_increment = []
    for i in range(num_incr):
        if bin_edges[i] < CNT_minz or bin_edges[i + 1] > CNT_maxz:
            vol_increment.append(simbox_x * simbox_y * (bin_edges[i + 1] - bin_edges[i]))
        else:
            vol_increment.append(math.pi * (used_radius ** 2) * (bin_edges[i + 1] - bin_edges[i]))

    # Add a new first column with the index+1 of the bin  add the volume, calculate the density.
    results_zd_df.reset_index(drop = True, inplace = True)   
    results_zd_df.insert(0, 'Bin', results_zd_df.index + 1) 
    results_zd_df['Volume'] = vol_increment                                                   
    results_zd_df['Density [u/Ang^3]'] = results_zd_df['Mass'] / results_zd_df['Volume']      
    
 
    # Mirror the data.
    mirror_dens_CNT = ddict.get_input('Do you want to mirror the results inside the CNT to increase the binning? (y/n) ', args)
    if mirror_dens_CNT == 'y':
        
        # First identify all the bins that are inside the CNT.
        inside_CNT = []            
        for i in range(num_incr):
            if bin_edges[i] >= CNT_minz and bin_edges[i + 1] <= CNT_maxz:
                inside_CNT.append(i)
        # Make dataframe with one column named 'Density [u/Ang^3]' and the same number of rows as the number of bins inside the CNT.
        inside_CNT_df = pd.DataFrame()
        inside_CNT_df['Density [u/Ang^3]'] = results_zd_df['Density [u/Ang^3]'][inside_CNT]
        
        # Take the first and last bin that are inside the CNT and average the density of the two bins. Then do this for the second and second last bin, and so on.
        for i in range(int(len(inside_CNT))):          
            inside_CNT_df['Density [u/Ang^3]'][inside_CNT[i]] = ((results_zd_df['Density [u/Ang^3]'][inside_CNT[i]] + results_zd_df['Density [u/Ang^3]'][inside_CNT[-i - 1]]) / 2).copy()
            results_zd_df.at[inside_CNT[i],'Density [u/Ang^3]'] = inside_CNT_df['Density [u/Ang^3]'][inside_CNT[i]]
            results_zd_df.at[inside_CNT[-i - 1],'Density [u/Ang^3]'] = inside_CNT_df['Density [u/Ang^3]'][inside_CNT[i]]
        ddict.printLog('')
        ddict.printLog(results_zd_df)
        ddict.printLog('')
    
    # Calculate the density in g/cm^3.
    results_zd_df['Density [g/cm^3]'] = results_zd_df['Density [u/Ang^3]'] * 1.66053907   

    center_to_zero=ddict.get_input('Do you want to set the center of the simulation box to zero? (y/n) ', args)   
    if center_to_zero == 'y':
        results_zd_df['Bin_center'] = results_zd_df['Bin_center'] - center_tube_z

    # Plot the data.
    zdens_data_plot = ddict.get_input('Do you want to plot the data? (y/n) ', args)
    if zdens_data_plot == 'y':
        fig, ax = plt.subplots()
        ax.plot(results_zd_df['Bin_center'], results_zd_df['Density [g/cm^3]'], '-' , label = 'Z density function', color = 'black')
        ax.set(xlabel = 'Distance from tube center [Ang]', ylabel = 'Density [g/cm^3]', title = 'Z density function')
        ax.grid()  
        fig.savefig("Axial_density.pdf")
        ddict.printLog('Z density function saved as Axial_density.pdf')
    
    # Save the data.
    results_zd_df.to_csv('Axial_density.csv', sep = ';', index = True, header = True, float_format = '%.5f')
    ddict.printLog('Z density data saved as Axial_density.csv')

    raw_data = ddict.get_input('Do you want to save the raw data? (y/n) ', args)
    if raw_data == 'y':
        df_zdens.to_csv('Axial_density_raw.csv', sep = ';', index = False, header = True, float_format = '%.5f')
        ddict.printLog('Raw data saved as Axial_density_raw.csv')