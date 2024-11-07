import numpy as np
import matplotlib.pyplot as plt
import os
import sys

plt.rcParams['axes.labelsize'] = 12  # Set the font size for x and y labels
plt.rcParams['axes.labelweight'] = 'bold'  # Make the x and y labels bold
plt.rcParams['xtick.labelsize'] = 12  # Set the font size for x tick labels
#plt.rcParams['xtick.labelweight'] = 'bold'  # Make the x tick labels bold
plt.rcParams['ytick.labelsize'] = 13  # Set the font size for y tick labels
#plt.rcParams['ytick.labelweight'] = 'bold'  # Make the y tick labels bold
plt.rcParams['axes.titlesize'] = 14  # Set the font size for the title
plt.rcParams['axes.titleweight'] = 'bold'  # Make the title bold

# Check if the correct number of arguments are provided
if len(sys.argv) != 4:
    print("Usage: script.py <run_id> <processor_number> <timestep>")
    sys.exit(1)

# Get command line arguments
run_id = sys.argv[1]
processor_number = int(sys.argv[2])
timestep = int(sys.argv[3])

# Format processor number and timestep with leading zeros
processor_number_str = f"{processor_number:03d}"
timestep_str = f"{timestep:09d}"

# Construct the filename
filename = f"{run_id}_p{processor_number_str}_t{timestep_str}_t.txt"
filepath = os.path.join("output", filename)

# Read the file content
with open(filepath, 'r') as file:
    lines = file.readlines()

# Extract nx, ny, nz
nx, ny, nz = map(int, lines[0].strip().split(','))

# Extract data
data = np.array(list(map(float, lines[1].strip().split(','))))

# Reshape data according to nx, ny, nz
data = data.reshape((nz, nx, ny))

# Create subplots
fig, axes = plt.subplots(3, 4, figsize=(20, 15))

# Plotting heatmaps and profiles for XZ plane
for i in range(3):
    y_value = i * (ny // 4)  # Adjust y_value for each row
    data_slice_xz = data[:, :, y_value]
    average_profile_xz = np.mean(data_slice_xz, axis=1)

    # Heatmap for XZ plane
    im_xz = axes[i, 0].imshow(data_slice_xz, cmap='RdBu_r', origin='lower', aspect='auto')
    fig.colorbar(im_xz, ax=axes[i, 0], label='Temperature')
    axes[i, 0].set_xlabel('X-axis')
    axes[i, 0].set_ylabel('Z-axis')
    axes[i, 0].set_title(f"XZ Heatmap (Y={y_value})")

    # Average vertical profile for XZ plane
    axes[i, 1].plot(average_profile_xz, np.arange(nz))
    axes[i, 1].set_xlabel('Average Temperature')
    axes[i, 1].set_ylabel('Z-axis')
    axes[i, 1].set_title(f"XZ Average Profile (Y={y_value})")

# Plotting heatmaps and profiles for YZ plane
for i in range(3):
    x_value = i * (nx // 4)  # Adjust x_value for each row
    data_slice_yz = data[:, x_value, :]
    average_profile_yz = np.mean(data_slice_yz, axis=1)

    # Heatmap for YZ plane
    im_yz = axes[i, 2].imshow(data_slice_yz, cmap='RdBu_r', origin='lower', aspect='auto')
    fig.colorbar(im_yz, ax=axes[i, 2], label='Temperature')
    axes[i, 2].set_xlabel('Y-axis')
    axes[i, 2].set_ylabel('Z-axis')
    axes[i, 2].set_title(f"YZ Heatmap (X={x_value})")

    # Average vertical profile for YZ plane
    axes[i, 3].plot(average_profile_yz, np.arange(nz))
    axes[i, 3].set_xlabel('Average Temperature')
    axes[i, 3].set_ylabel('Z-axis')
    axes[i, 3].set_title(f"YZ Average Profile (X={x_value})")

plt.tight_layout()
#plt.show()
plt.savefig(f"plots/{run_id}/slices_{run_id}_p{processor_number_str}_t{timestep_str}_t.png", dpi=150)
