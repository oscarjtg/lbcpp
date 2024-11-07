import numpy as np
import matplotlib.pyplot as plt
import os
import sys

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

# Create grid.
X, Y, Z = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz))

# Extract data
data = np.array(list(map(float, lines[1].strip().split(','))))

# Reshape data according to nx, ny, nz
data = data.reshape((nx, ny, nz))

kw = {
    'vmin': data.min(),
    'vmax': data.max(),
    'levels': np.linspace(data.min(), data.max(), 10),
}

# Create a figure with 3D ax
fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(111, projection='3d')

# Plot contour surfaces
_ = ax.contourf(
    X[:, :, -1], Y[:, :, -1], data[:, :, -1],
    zdir='z', offset=0, **kw
)
_ = ax.contourf(
    X[0, :, :], data[0, :, :], Z[0, :, :],
    zdir='y', offset=0, **kw
)
C = ax.contourf(
    data[:, -1, :], Y[:, -1, :], Z[:, -1, :],
    zdir='x', offset=X.max(), **kw
)
# --


# Set limits of the plot from coord limits
xmin, xmax = X.min(), X.max()
ymin, ymax = Y.min(), Y.max()
zmin, zmax = Z.min(), Z.max()
ax.set(xlim=[xmin, xmax], ylim=[ymin, ymax], zlim=[zmin, zmax])

# Plot edges
edges_kw = dict(color='0.4', linewidth=1, zorder=1e3)
ax.plot([xmax, xmax], [ymin, ymax], 0, **edges_kw)
ax.plot([xmin, xmax], [ymin, ymin], 0, **edges_kw)
ax.plot([xmax, xmax], [ymin, ymin], [zmin, zmax], **edges_kw)

# Set labels and zticks
ax.set(
    xlabel='X',
    ylabel='Y',
    zlabel='Z',
    #zticks=[0, -150, -300, -450],
)

# Set zoom and angle view
ax.view_init(40, -30, 0)
ax.set_box_aspect(None, zoom=0.9)

# Colorbar
fig.colorbar(C, ax=ax, fraction=0.02, pad=0.1, label='Name [units]')

#plt.show()
plt.savefig(f"plots/boxplot_{run_id}_p{processor_number_str}_t{timestep_str}_t.png", dpi=150)