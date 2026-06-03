"""
Created on Tu Jun 02 2026
"""

import numpy as np
import matplotlib.pyplot as plt

# path to save figures to 
path = 'C:\\Users\\Leoni\\OneDrive\\Dokumente\\Uni-Aspire-Leonie\\Light Tweezers Auswertung\\Lighttweezers_auswertung\\'
# path for the files
path_files = 'C:\\Users\\Leoni\\OneDrive\\Dokumente\\Uni-Aspire-Leonie\\Light Tweezers Auswertung\\Lighttweezers_auswertung\\Daten\\teil1\\'

# function to extract data
def extract_data(filepath):
    """
    Arguments:
        filepath (str): Path to the data files

    Returns:
        position of the particle in x and y direction in pixels
    """
    x_position = []
    y_position = []

    with open(filepath, "r") as file:
        lines = file.readlines()
        
        # skip the header row by starting the loop at index 1
        for line in lines[1:]: 
            # strip the newline character and split by the csv semicolon
            columns = line.strip().split(';')
            
            # make sure the line isn't empty and has enough columns
            if len(columns) > 5: 
                # convert to float() for plotting later
                x_position.append(float(columns[4])) # 5th column
                y_position.append(float(columns[5])) # 6th column
            
    return x_position, y_position



# extract and plot data
file_names = ['first_particle.csv', 'second_particle.csv', 'third_particle.csv', 'vergleich_trapped.csv', 'much_moving.csv']
colors = ['darkred', 'forestgreen', 'darkslateblue', 'darkorange', 'magenta']
labels = [' 1', ' 2', ' 3', ' trapped', ' moving']

# parameters
T = 293.15 # K (20°C)
R = 1.03 # micro meter
k_B = 1.38064852 * 10**(-23) # m²kgs^(-2)K^(-1)

n_water = 1.002 * 10**(-3) # Pa*s
s_literatur = ( (2*k_B*T) / (3*np.pi*R* n_water) ) 

print(f'slope literatur: s={s_literatur} μm²/s')

i = 0
for file in file_names:
    
    # initializing figure 
    fig1 = plt.figure('Trajectory particle' + labels[i])
    
    # exctract data from csv file
    x_pos, y_pos = extract_data(path_files + file_names[i])
    
    # conversion pixels to micrometer using sample dimensions [130, 98] micrometer and pixel range [1440, 1080] px
    x_pos = [x * 0.0903 for x in x_pos]
    y_pos = [y * 0.0907 for y in y_pos]
    
    # plot data points
    plt.plot(x_pos, y_pos, color=colors[i], label='particle'+labels[i]) # alpha=0.5 = transparency

    # star on the FIRST value (index 0)
    plt.plot(x_pos[0], y_pos[0], marker='*', color='black', markersize=10, markeredgecolor='black', zorder=5)    
    # circle on the LAST value (index -1)
    plt.plot(x_pos[-1], y_pos[-1], marker='o', color='black', markersize=6, markeredgecolor='black', zorder=5)
    
    # find midpoint of the data
    x_midpoint = (max(x_pos) + min(x_pos)) / 2
    y_midpoint = (max(y_pos) + min(y_pos)) / 2
    
    # axis_range = 5 # good for moving particle
    axis_range = 0.4
    # set the axis limits to midpoint ± 0.4 micrometers (giving a total span of 0.8) 
    plt.xlim(x_midpoint - axis_range, x_midpoint + axis_range)
    plt.ylim(y_midpoint - axis_range, y_midpoint + axis_range)   
      
    # define plot
    plt.xlabel(r'x position [$\mu$m]')
    plt.ylabel(r'y position [$\mu$m]')
    plt.grid(True, alpha=0.6) # alpha 0.0 = fully transparent, 1.0 = fully opaque
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(path + 'Trajectory particle' + labels[i], bbox_inches='tight')
    plt.show()
    
    # calculate displacement in x,y direction relative to the very first point (x_pos[0],y_pos[0])
    x_squared_displacements = []
    y_squared_displacements = []
    
    steps = range(0, len(x_pos))
    for j in steps:
        x_squared = (x_pos[j] - x_pos[0]) ** 2
        y_squared = (y_pos[j] - y_pos[0]) ** 2
        
        x_squared_displacements.append(x_squared)
        y_squared_displacements.append(y_squared)
    
    # radial squared displacement: r² = x² + y²
    r_squared_displacements = [x + y for x, y in zip(x_squared_displacements, y_squared_displacements)]
        
    # linear fits (deg=1 is linear: y = mx + c)
    # np.polyfit returns [slope, intercept]
    slope_x, intercept_x = np.polyfit(steps, x_squared_displacements, 1)
    slope_y, intercept_y = np.polyfit(steps, y_squared_displacements, 1)
    slope_r, intercept_r = np.polyfit(steps, r_squared_displacements, 1)
    
    # fit lines for plotting
    fit_line_x = slope_x * steps + intercept_x
    fit_line_y = slope_y * steps + intercept_y
    fit_line_r = slope_r * steps + intercept_r
    
    n_eff = ( (2*k_B*T) / (3*np.pi*R* s_literatur) ) 
    
    # print the slopes 
    print(f"--- Particle {labels[i]} ---")
    print(f"X Fit Slope: {slope_x:.6f} μm²/step")
    print(f"Y Fit Slope: {slope_y:.6f} μm²/step")
    print(f"Radial Fit Slope: {slope_r:.6f} μm²/step")
    print(f"effective viscosity: {n_eff:.16f} μPa*s \n")
    
    # initializing figure
    fig2 = plt.figure('Squared displacements' + labels[i])
    
    # plot data points
    plt.plot(steps, x_squared_displacements, color='gold', label='x²(t) particle'+labels[i])
    plt.plot(steps, y_squared_displacements, color='darkorange', label='y²(t) particle'+labels[i])
    
    # plot linear fits
    plt.plot(steps, fit_line_x, color=colors[i], lw=2, linestyle='--', label=f'x fit (s={slope_x:.2e} μm²/step )')
    plt.plot(steps, fit_line_y, color=colors[i], lw=2, linestyle='--', alpha=0.5, label=f'y fit (s={slope_y:.2e} μm²/step)')
    plt.plot(steps, fit_line_r, color='black', lw=2, linestyle='--', label=f'radial fit (s={slope_r:.2e} μm²/step)')
    
    # define plot
    plt.xlabel(r'steps')
    plt.ylabel(r'squared displacement [$\mu$m²]')
    plt.grid(True, alpha=0.6) # alpha 0.0 = fully transparent, 1.0 = fully opaque
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(path + 'Squared displacements' + labels[i], bbox_inches='tight')
    plt.show()
    
    

    i += 1
     
    

