import os

# Loop from 0 to 100000 in steps of 1000
for n in range(0, 100001, 1000):
    # Construct the command
    command = f"python processing_scripts/plot_temperature_slices.py overnight3d 0 {n}"
    # Execute the command
    os.system(command)
