import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.special import erfc
plt.rcParams.update({
    'font.size': 12,          # Default font size
    'axes.titlesize': 14,     # Title font size
    'axes.labelsize': 12,     # Axis label font size
    'xtick.labelsize': 10,    # X-axis tick label font size
    'ytick.labelsize': 10,    # Y-axis tick label font size
    'legend.fontsize': 8,    # Legend font size
    'figure.titlesize': 16    # Figure title font size
})

lattice_types = ["D2Q5", "D2Q9"]
operators = ["SRTxx", "TRT12", "TRT06", "TRT04"]
precisions = ["F64", "F32"]
timesteps = [47000]
linestyles = ["solid", "dashed", "dotted", "dashdot", (0, (1, 10))]
#colours = ["black", "navy", "blue", "green", "darkgreen"]
colours = ["lightgray", "blue", "green", "red", "orange"]

def plate_diffusion(x, y):
    Cp = 1.0
    D = 0.07407
    u0 = 0.05
    arg = y / np.sqrt(4*D*x / u0)
    return Cp * erfc(arg)

# Create subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 5))

nx = 800
ny = 80
X, Y = np.meshgrid(range(nx), range(ny))
analytic = plate_diffusion(X+0.5, Y+0.5)
#ax.contour(X, Y, analytic, label="analytic", levels=[0.2, 0.4, 0.6, 0.8], linestyle=linestyles[0], colors=colours[0])

num = 0;
operator = operators[num]
run_id = lattice_types[0] + "_" + str(nx) + "x1x" + str(ny) + "_" + operator + "_" + precisions[0]
timestep_str = f"{timesteps[0]:09d}"
filename = f"{run_id}_p000_t{timestep_str}_c.txt"
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
data_slice = data[:, :, 0]

# Plot contour plot of data_slice
X, Y = np.meshgrid(range(nx), range(nz))
X = X / nz
Y = Y / nz
ax1.contour(X, Y, analytic, label=operator,  linestyles=linestyles[num], colors=colours[num], linewidths=3, levels=[0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8])
contour_plot = ax1.contour(X, Y, data_slice, label=operator,  linestyles="solid", colors=colours[num+1], levels=[0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8])
#ax1.set_aspect("equal")
ax1.clabel(contour_plot, contour_plot.levels, fontsize=12)
residual = data_slice-analytic
maxerror = max(residual.min(), residual.max(), key=abs)
pcm = ax2.pcolormesh(X, Y, residual, cmap="RdBu_r", vmin=-maxerror, vmax=maxerror)
pcm = ax3.pcolormesh(X, Y, residual, cmap="RdBu_r", vmin=-0.2, vmax=0.2)
pcm = ax4.pcolormesh(X, Y, residual, cmap="RdBu_r", vmin=-0.2, vmax=0.2)
#tol = 0.1
#ax2.scatter(X[abs(residual)>0.05], Y[abs(residual)>0.05])
cbar2 = plt.colorbar(pcm, ax=ax2, label=r"computed $-$ analytic")
cbar3 = plt.colorbar(pcm, ax=ax3, label=r"computed $-$ analytic")
cbar4 = plt.colorbar(pcm, ax=ax4, label=r"computed $-$ analytic")


# Set the limits of ax1 to be between 200 and 240
ax1.set_xlim(-0.1, 10.1)
ax1.set_ylim(-0.1, 1.1)
ax2.set_xlim(-0.1, 10.1)
ax2.set_ylim(-0.1, 1.1)
ax3.set_xlim(-0.01, 0.2)
ax3.set_ylim(-0.01, 0.1)
ax4.set_xlim(9.75, 10.02)
ax4.set_ylim(-0.1, 1.1)

# Set axis labels
ax4.set_xlabel(r"$x / N_y$")
ax1.set_ylabel(r"$y / N_y$")
ax3.set_xlabel(r"$x / N_y$")
ax3.set_ylabel(r"$y / N_y$")

plt.tight_layout()
plt.savefig(f"plots/{run_id}_contours.pdf")
plt.show()

