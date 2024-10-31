import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys

def collect_data(runId, timestep):
    # Initialize an empty list to store the data arrays
    data_arrays = []

    # Iterate over the process numbers in the 3x3 grid
    for process in range(9):
        # Construct the filename pattern for the given runId, process, and timestep
        filename_pattern = f"./output/{runId:05d}_p{process:03d}_t{timestep:09d}_r.csv"
        
        # Find the file matching the pattern
        files = glob.glob(filename_pattern)
        
        if files:
            # Load the CSV file into a numpy array and append to the list
            data = np.loadtxt(files[0], delimiter=',')
            data_arrays.append(data)
        else:
            print(f"File not found: {filename_pattern}")
            return None

    # Combine the 9 data arrays into one large 2D array
    combined_data = np.block([
        [data_arrays[0], data_arrays[1], data_arrays[2]],
        [data_arrays[3], data_arrays[4], data_arrays[5]],
        [data_arrays[6], data_arrays[7], data_arrays[8]]
    ])

    return combined_data

def plot_heatmap(data, runId, timestep):
    plt.figure(figsize=(10, 10))
    plt.imshow(data, cmap='RdBu', aspect='auto')
    plt.colorbar(label='Value')
    plt.title(f'Heatmap for Run {runId}, Timestep {timestep}')
    
    # Save the plot as PDF and high resolution PNG
    plt.savefig(f'process/density_run{runId}_timestep{timestep}.pdf')
    plt.savefig(f'process/density_run{runId}_timestep{timestep}.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <runId> <timestep>")
        sys.exit(1)

    runId = int(sys.argv[1])
    timestep = int(sys.argv[2])

    data = collect_data(runId, timestep)
    
    if data is not None:
        plot_heatmap(data, runId, timestep)
