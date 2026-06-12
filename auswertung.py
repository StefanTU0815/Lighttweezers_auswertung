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
    time_steps = []

    with open(filepath, "r") as file:
        lines = file.readlines()
        
        # skip the header row by starting the loop at index 1
        for line in lines[1:]: 
            # replace , with . 
            current_line = line.strip().replace(",",".")
            # strip the newline character and split  the csv by semicolon
            columns = current_line.strip().split(';')
              
            
            # make sure the line isn't empty and has enough columns
            if len(columns) > 5: 
                # convert to float() for plotting later
                time_steps.append(float(columns[1])) # 2nd column = time steps in s
                x_position.append(float(columns[4])) # 5th column = x pos in px
                y_position.append(float(columns[5])) # 6th column = y pos in px
            
    return x_position, y_position, np.array(time_steps)


# # extract and plot data
# # files particles that don't move much
# file_names = ['first_particle.csv', 'second_particle.csv', 'third_particle.csv', 'vergleich_trapped.csv']
# colors = ['darkred', 'forestgreen', 'darkslateblue', 'magenta']
# labels = [' 1', ' 2', ' 3', ' trapped']
# axis_range = 0.4
# files particles that move more
file_names = ['much_moving.csv', 'moving_1.csv', 'moving_2.csv', 'moving_3.csv', 'vergleich_trapped.csv']
colors = ['yellow', 'darkred', 'forestgreen', 'darkslateblue', 'magenta']
axis_range = 9
labels = [' extra', ' 1', ' 2', ' 3', ' trapped']

# parameters
T = 293.15 # K room temperatur (20°C)
R = 1.03 * 10**(-6) # m radius of the particle
k_B = 1.38064852 * 10**(-23) # m²kgs^(-2)K^(-1) Bolzman const.

# literature values
n_water = 1.002 * 10**(-3) # Pa*s
s_literatur = ( (2*k_B*T) / (3*np.pi*R* n_water) ) 
print(f'slope literatur: s={s_literatur:.3e} μm²/s \n')

viscosities = []
i = 0
for file in file_names:
    
    # --------------------- particle trajectories -----------------------------
    # initializing particle trajectories figure 
    fig1 = plt.figure('Trajectory particle' + labels[i])
    
    # exctract data from csv file
    x_pos, y_pos, time = extract_data(path_files + file_names[i])
    
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
    
    # set the axis limits to midpoint ± x micrometers (giving a total span of 2x micrometers with x = axis_range) 
    plt.xlim(x_midpoint - axis_range, x_midpoint + axis_range)
    plt.ylim(y_midpoint - axis_range, y_midpoint + axis_range)   
      
    # define trajectories plot
    plt.xlabel(r'x position [$\mu$m]')
    plt.ylabel(r'y position [$\mu$m]')
    plt.grid(True, alpha=0.6) # alpha 0.0 = fully transparent, 1.0 = fully opaque
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(path + 'Trajectory particle' + labels[i], bbox_inches='tight')
    plt.show()
    
    # ---------------------- squared displacements ----------------------------
    # calculate  squared displacement in x,y direction relative to the very first point (x_pos[0],y_pos[0])
    x_squared_displacements = []
    y_squared_displacements = []
    
    x_MSD = []
    y_MSD = []
    
    steps = range(0, len(x_pos))
    for j in steps:
        # calculate squared displacements
        x_squared = (x_pos[j] - x_pos[0]) ** 2
        y_squared = (y_pos[j] - y_pos[0]) ** 2
        
        x_squared_displacements.append(x_squared)
        y_squared_displacements.append(y_squared)
        
        # calculate mean squared displacements
        mean_x_squared = np.sum(x_squared_displacements) / len(x_squared_displacements)
        mean_y_squared = np.sum(y_squared_displacements) / len(y_squared_displacements)
        
        x_MSD.append(mean_x_squared)
        y_MSD.append(mean_y_squared)
    
    # radial squared displacement: r² = x² + y²
    r_squared_displacements = [x + y for x, y in zip(x_squared_displacements, y_squared_displacements)]
        
    # linear fits (deg=1 is linear: y = s*x + c)
    # np.polyfit returns [slope, intercept]
    slope_x, intercept_x = np.polyfit(time, x_squared_displacements, 1)
    slope_y, intercept_y = np.polyfit(time, y_squared_displacements, 1)
    slope_r, intercept_r = np.polyfit(time, r_squared_displacements, 1)
    
    # fit lines for plotting
    fit_line_x = slope_x * time + intercept_x
    fit_line_y = slope_y * time + intercept_y
    fit_line_r = slope_r * time + intercept_r
    
    slope_r = slope_r * 10**(-12) # μm²/s -> m²/s
    n_eff = ( (2*k_B*T) / (3*np.pi*R* slope_r) ) # Pa*s effective viscosity
    viscosities.append(n_eff)
    n_eff = n_eff * (10**3) # Pa*s -> mPa*s
    
    # print the slopes 
    print(f"--- Particle {labels[i]} ---")
    print(f"X Fit Slope: {slope_x:.3e} μm²/step")
    print(f"Y Fit Slope: {slope_y:.3e} μm²/step")
    print(f"Radial Fit Slope: {slope_r:.3e} μm²/step")
    print(f"effective viscosity: {n_eff:.3e} mPa*s \n")
    
    # initializing squared displacements figure
    fig2 = plt.figure('Squared displacements' + labels[i])
    
    # plot data points
    plt.plot(time, x_squared_displacements, color='gold', label='x²(t) particle'+labels[i])
    plt.plot(time, y_squared_displacements, color='darkorange', label='y²(t) particle'+labels[i])
    
    # plot linear fits
    plt.plot(time, fit_line_x, color=colors[i], lw=2, linestyle='--', label=f'x fit (s={slope_x:.2e} μm²/s )')
    plt.plot(time, fit_line_y, color=colors[i], lw=2, linestyle='--', alpha=0.5, label=f'y fit (s={slope_y:.2e} μm²/s)')
    plt.plot(time, fit_line_r, color='black', lw=2, linestyle='--', label=f'radial fit (s={slope_r:.2e} μm²/s)')
    
    # define squared displacements plot
    plt.xlabel(r'time [s]')
    plt.ylabel(r'squared displacement [$\mu$m²]')
    plt.grid(True, alpha=0.6) # alpha 0.0 = fully transparent, 1.0 = fully opaque
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(path + 'Squared displacements particle' + labels[i], bbox_inches='tight')
    plt.show()
    
    # ------------------ mean squared displacements --------------------------   
    # initializing mean squared displacements figure
    fig3 = plt.figure('Mean squared displacements' + labels[i])
    
    # plot data points
    plt.plot(time, x_MSD, color='gold', label='<x²(t)> particle'+labels[i])
    plt.plot(time, y_MSD, color='darkorange', label='<y²(t)> particle'+labels[i])
    
    # define mean squared displacements plot
    plt.xlabel(r'time [s]')
    plt.ylabel(r'mean squared displacement [$\mu$m²]')
    plt.grid(True, alpha=0.6) # alpha 0.0 = fully transparent, 1.0 = fully opaque
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(path + 'Mean squared displacements particle' + labels[i], bbox_inches='tight')
    plt.show()
    
    i += 1
     
# calculate mean of viscosities    
mean_n = (viscosities[1] + viscosities[2] + viscosities[3]) / 3
std_n = np.std([viscosities[1], viscosities[2], viscosities[3]])

print(f"mean viscosity: {mean_n:.3e} Pa*s \n")
print(f"std viscosity: {std_n:.3e} Pa*s \n")

# calculate frictional force
F_fric_1 = []
F_fric_2 = []
# maximal velocity 
v_max_1 = [125, 191] # mm/s maximal velocity particle 1 with laser current [70 mA, 100 mA]
v_max_2 = [90, 170] # mm/s maximal velocity particle 2 with laser current [70 mA, 100 mA]

for i in range(0, len(v_max_1)):
    F_fric_p1 = 6 * np.pi * mean_n * R * (v_max_1[i] * 10**(-3))
    F_fric_p2 = 6 * np.pi * mean_n * R * (v_max_2[i] * 10**(-3))
    
    F_fric_1.append(F_fric_p1)
    F_fric_2.append(F_fric_p2)
    
print(f"frictional force particel 1 I=70mA: {F_fric_1[0]:.3e} N")
print(f"frictional force particel 1 I=100mA: {F_fric_1[1]:.3e} N")
print(f"frictional force particel 2 I=70mA: {F_fric_2[0]:.3e} N")
print(f"frictional force particel 2 I=100mA: {F_fric_2[1]:.3e} N \n")