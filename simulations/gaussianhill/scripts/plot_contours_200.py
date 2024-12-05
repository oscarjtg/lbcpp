import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams.update({
    'font.size': 12,          # Default font size
    'axes.titlesize': 14,     # Title font size
    'axes.labelsize': 12,     # Axis label font size
    'xtick.labelsize': 10,    # X-axis tick label font size
    'ytick.labelsize': 10,    # Y-axis tick label font size
    'legend.fontsize': 8,    # Legend font size
    'figure.titlesize': 16    # Figure title font size
})

def gaussian_hill(x, y, t):
    sigma0 = 10
    D = 0.0053
    x0 = 200
    y0 = 200
    u0 = 0.1
    v0 = 0.1
    C0 = 1.0
    sigmaD = np.sqrt(2*D*t)
    numerator = (x-x0-u0*t)**2 + (y-y0-v0*t)**2
    denominator = 2*(sigma0**2 + sigmaD**2)
    return sigma0**2 / (sigma0**2 + sigmaD**2) * C0 * np.exp(-numerator / denominator)

lattice_types = ["D2Q5", "D2Q9"]
operators = ["SRTxx", "TRT12", "TRT06", "TRT04"]
precisions = ["F64", "F32"]
timesteps = [200, 5120]
linestyles = ["solid", "dashed", "dotted", "dashdot", "dotted"]#, (0, (1, 10))]
#colours = ["black", "navy", "blue", "green", "darkgreen"]
colours = ["lightgrey", "blue", "green", "red", "orange"]

# Create subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 8))

X, Y = np.meshgrid(range(512), range(512))
analytic = gaussian_hill(X, Y, 200)
contour_plot = ax1.contour(X, Y, analytic, label="analytic", levels=[0.2, 0.4, 0.6, 0.8], linestyles=linestyles[0], colors=colours[0], linewidth=4)
ax1.clabel(contour_plot, contour_plot.levels, fontsize=12)
ax2.plot(analytic[:, 220], linestyle=linestyles[0], color=colours[0], linewidth=4)

for (num, operator) in enumerate(operators):
    run_id = lattice_types[0] + "_" + operator + "_" + precisions[0]
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
    contour_plot = ax1.contour(X, Y, data_slice, label=operator,  linestyles=linestyles[num+1], colors=colours[num+1], levels=4) #levels=[0.2, 0.4, 0.6, 0.8],

    ax2.plot(data_slice[:, 220], linestyle=linestyles[num+1], color=colours[num+1])

# Set the limits of ax1 to be between 200 and 240
ax1.set_xlim(180, 260)
ax1.set_ylim(180, 260)

ax2.set_xlim(180, 260)
ax2.set_ylim(0, 1)

ax1.set_ylabel("j")
ax2.set_ylabel(r"$C(i, j=220)$")
ax1.set_aspect("equal")




X, Y = np.meshgrid(range(512), range(512))
analytic = gaussian_hill(X, Y, 200)
ax3.contour(X, Y, analytic, label="analytic", levels=[0.2, 0.4, 0.6, 0.8], linestyles=linestyles[0], colors=colours[0], linewidth=4)
ax4.plot(analytic[:, 220], linestyle=linestyles[0], color=colours[0], linewidth=4)

for (num, operator) in enumerate(operators):
    run_id = lattice_types[1] + "_" + operator + "_" + precisions[0]
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
    contour_plot = ax3.contour(X, Y, data_slice, label=operator, linestyles=linestyles[num+1], colors=colours[num+1], levels=4) # levels=[0.2, 0.4, 0.6, 0.8],
    ax4.plot(data_slice[:, 220], linestyle=linestyles[num+1], color=colours[num+1])

# Set the limits of ax3 to be between 200 and 240
ax3.set_xlim(180, 260)
ax3.set_ylim(180, 260)

ax4.set_xlim(180, 260)
ax4.set_ylim(0, 1)

ax3.set_xlabel("i")
ax3.set_ylabel("j")
ax4.set_xlabel("i")
ax4.set_ylabel(r"$C(i, j=220)$")

ax3.set_aspect("equal")

plt.tight_layout()

plt.savefig("plots/gaussianhill_contours_t200.pdf")
plt.show()

