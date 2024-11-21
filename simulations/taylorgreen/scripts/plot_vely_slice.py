import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 5:
    print("Usage: script.py <run_id> <processor_number> <timestep> <y_value>")
    sys.exit(1)

# Get command line arguments
run_id = sys.argv[1]
processor_number = int(sys.argv[2])
timestep = int(sys.argv[3])
y_value = int(sys.argv[4])

# Format processor number and timestep with leading zeros
processor_number_str = f"{processor_number:03d}"
timestep_str = f"{timestep:09d}"

# Construct the filename
filename = f"{run_id}_p{processor_number_str}_t{timestep_str}_v.txt"
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

# Extract the slice for the given y value
data_slice = data[:, :, y_value]

# Calculate the average vertical Y-Velocity profile
average_profile = np.mean(data_slice, axis=1)

# Create subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Plotting the heatmap
im = ax1.imshow(data_slice, cmap='viridis', origin='lower', aspect='auto')
fig.colorbar(im, ax=ax1, label='Y-Velocity')
ax1.set_xlabel('X-axis')
ax1.set_ylabel('Z-axis')
ax1.set_title(f"Heatmap for Run ID: {run_id}, Processor: {processor_number}, Timestep: {timestep}, Y-value: {y_value}")

# Plotting the average vertical Y-Velocity profile
ax2.plot(average_profile, np.arange(nz))
ax2.set_xlabel('Average Y-Velocity')
ax2.set_ylabel('Z-axis')
ax2.set_title('Average Vertical Y-Velocity Profile')

plt.tight_layout()
plt.show()
plt.savefig(f"plots/{run_id}_p{processor_number_str}_t{timestep_str}_v.png", dpi=150)