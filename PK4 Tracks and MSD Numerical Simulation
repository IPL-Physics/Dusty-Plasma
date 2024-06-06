"""
Created on Tue Nov  7 13:47:33 2023

@author: Bradley Andrew     email: bra0016@auburn.edu
Soli Deo Gloria
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t

# Parameters
n_particles = 100  # Number of particles
T = 80.0          # Total simulation time
dt = 0.1          # Time step

area_size = 1000   # Size of the square area

# Calculate the number of particles per row and per column
particles_per_row = int(np.sqrt(n_particles))
particles_per_col = particles_per_row

# Create an array of initial positions for exactly n_particles
initial_positions = np.array([[x, y] for x in np.linspace(0, area_size, particles_per_row) for y in np.linspace(0, area_size, particles_per_col)][:n_particles])

# Extend the initial positions t
grid_positions = np.vstack([initial_positions])
# Random Initialial particle positions
rand_positions = position_scale = 50.0*np.random.rand(n_particles, 2)  # Random initial positions in [0, 1] for x and y

# Define your custom distribution function
def custom_distribution(x, sigma):
    return np.exp(-x**2 / (2 * sigma**2))

# Scaling factor for initial positions
position_scale = 1.0  # Adjust this value to increase or decrease the spread
 

sigma_x = 1 # Standard deviation in x-direction
sigma_y = 1 # Standard deviation in y-direction

q_x = 1.7 #nonextensive parameter in the x direction (cant have q=1 but you can make it really close)
q_y = 1.7 #nonextensive parameter in the y direction (max value 2, though even this is oversampled)

df_x = (3-q_x)/(q_x-1)  #change for students t description instead of q-gaussian
df_y = (3-q_y)/(q_y-1)

# Lists to store particle tracks
x_tracks = [[] for _ in range(n_particles)]
y_tracks = [[] for _ in range(n_particles)]

positions=grid_positions
# Perform the simulation
for tstep in np.arange(0, T, dt):
    for i in range(n_particles):
        # Append current positions to tracks
        x_tracks[i].append(positions[i, 0])
        y_tracks[i].append(positions[i, 1])


        # Update particle positions based on the custom distribution function
        x0 = positions[i, 0]
        y0 = positions[i, 1]
        
        # Update particle positions with normal and Cauchy noise
        #noise_x = np.random.normal(0, sigma_x) + np.random.standard_cauchy() * 0.001
        #noise_x = np.random.standard_cauchy() * sigma_x
        noise_x = t(df_x, loc=0, scale=sigma_x).rvs(size=1).reshape(-1, 1)
        #noise_y = np.random.normal(0, sigma_y)
        #noise_y = sigma_y*np.random.standard_cauchy()  + np.random.normal(0, 2)
        noise_y = t(df_y, loc=0, scale=sigma_y).rvs(size=1).reshape(-1, 1)
        
        
        positions[i, 0] += noise_x
        positions[i, 1] += noise_y

# Plot particle tracks
for i in range(n_particles):
    plt.plot(x_tracks[i], y_tracks[i], lw=0.5)

# Calculate the maximum values in tracks
max_x = np.max([np.max(track) for track in x_tracks])
max_y = np.max([np.max(track) for track in y_tracks])


#plt.xlim(-max_x , max_x )
#plt.ylim(-max_y , max_y )
plt.xlim(-100, area_size+100)
plt.ylim(-100, area_size+100)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Particle Tracks with Different Noise Distributions')
plt.show()


#%%
# Calculate Mean Squared Displacement as function of time
max_time = int(T / dt)
msd_values = np.zeros(max_time)
squared_displacements = []

for j in range(n_particles):
    x0, y0 = x_tracks[j][0], y_tracks[j][0]
    x_t, y_t = x_tracks[j], y_tracks[j]
    squared_displacement = [(x - x0)**2 + (y - y0)**2 for x, y in zip(x_t, y_t)]
    squared_displacements.append(squared_displacement)


# Create a time array for plotting
time_steps = np.arange(800)
average_msd = np.mean(squared_displacements, axis=0)

# Plot the Mean Squared Displacement (MSD) as a function of time displacement
plt.figure()
plt.plot(time_steps, average_msd, lw=2)

plt.xlabel('Time Displacement')
plt.ylabel('Mean Squared Displacement (MSD)')
plt.title('Mean Squared Displacement vs. Time Displacement')
plt.show()

#%%
import pandas as pd

# Create a list to hold the data frames for each particle's track
dfs = []

# Loop through all particles' tracks
for i, (x_track, y_track) in enumerate(zip(x_tracks, y_tracks), start=1):
    # Create DataFrame for each particle's track
    df_particle = pd.DataFrame({
        'Trajectory': i,
        'Frame': range(1, len(x_track) + 1),
        'x': x_track,
        'y': y_track
    })
    # Append DataFrame to the list
    dfs.append(df_particle)

# Concatenate all DataFrames into a single DataFrame
df_tracks = pd.concat(dfs, ignore_index=True)

# Specify the desired directory and file name for the Excel file
excel_file_path = 'C:/Users/YOUR_USERNAME/Desktop/tracks_data.xlsx'  

# Export the DataFrame to the Excel file
df_tracks.to_excel(excel_file_path, index=False)  # Set index=False to exclude row numbers from the Excel file

#%%

#MSD as a function of time delay averaged over all particles

def calculate_msd(tracks, dt, max_delay):
    n_particles = len(tracks)
    x_msd_values = np.zeros(max_delay)
    y_msd_values = np.zeros(max_delay)
    total_msd_values = np.zeros(max_delay)
    
    for delay in range(max_delay):
        x_squared_displacements = []
        y_squared_displacements = []
        total_squared_displacements = []
        
        for particle_idx in range(n_particles):
            x_track, y_track = tracks[particle_idx]

            x_squared_displacement = []
            y_squared_displacement = []
            total_squared_displacement = []
            
            for t in range(len(x_track) - delay):
                x_displacement = x_track[t + delay] - x_track[t]
                y_displacement = y_track[t + delay] - y_track[t]
                
                x_squared_displacement.append(x_displacement**2)
                y_squared_displacement.append(y_displacement**2)
                total_squared_displacement.append((x_displacement**2) + (y_displacement**2))
                
            x_squared_displacements.append(np.mean(x_squared_displacement))
            y_squared_displacements.append(np.mean(y_squared_displacement))
            total_squared_displacements.append(np.mean(total_squared_displacement))
        
        x_msd_values[delay] = np.mean(x_squared_displacements)
        y_msd_values[delay] = np.mean(y_squared_displacements)
        total_msd_values[delay] = np.mean(total_squared_displacements)
    
    return x_msd_values, y_msd_values, total_msd_values

# Assume x_tracks and y_tracks are already defined

# Combine x_tracks and y_tracks into a list of tracks
tracks = [(x_track, y_track) for x_track, y_track in zip(x_tracks, y_tracks)]

# Define parameters
dt = 0.1  # Time step size
max_delay = 800  # Maximum time delay

# Calculate MSD for x and y coordinates, and total MSD
x_msd_values, y_msd_values, total_msd_values = calculate_msd(tracks, dt, max_delay)
#%%

# Plot x_msd(delay)
plt.figure()
plt.loglog(range(max_delay), x_msd_values, label='x_msd')
plt.xlabel('Time Delay')
plt.ylabel('Mean Squared Displacement (x)')
plt.title('Mean Squared Displacement (x) vs. Time Delay')
plt.legend()
plt.show()

# Plot y_msd(delay)
plt.loglog(range(max_delay), y_msd_values, label='y_msd')
plt.xlabel('Time Delay')
plt.ylabel('Mean Squared Displacement (y)')
plt.title('Mean Squared Displacement (y) vs. Time Delay')
plt.legend()
plt.show()

# Plot total_msd(delay)
plt.loglog(range(max_delay), total_msd_values, label='total_msd')
plt.xlabel('Time Delay')
plt.ylabel('Total Mean Squared Displacement')
plt.title('Total Mean Squared Displacement vs. Time Delay')
plt.legend()
plt.show()


plt.loglog(range(max_delay),range(max_delay),label='linear')
plt.legend()
plt.show()
#%%

# Assuming x_tracks and y_tracks are already defined and represent displacements

# Calculate velocities for x and y directions
dt = 1  # Time step size
x_velocities = [(x_tracks[i][t + 1] - x_tracks[i][t]) / dt for i in range(len(x_tracks)) for t in range(len(x_tracks[i]) - 1)]
y_velocities = [(y_tracks[i][t + 1] - y_tracks[i][t]) / dt for i in range(len(y_tracks)) for t in range(len(y_tracks[i]) - 1)]

# Plot histograms of x and y velocities
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.hist(x_velocities, bins=50, color='blue', alpha=0.7)
plt.xlabel('X Velocity')
plt.ylabel('Frequency')
plt.title('Histogram of X Velocities')

plt.subplot(1, 2, 2)
plt.hist(y_velocities, bins=50, color='red', alpha=0.7)
plt.xlabel('Y Velocity')
plt.ylabel('Frequency')
plt.title('Histogram of Y Velocities')

plt.tight_layout()
plt.show()
