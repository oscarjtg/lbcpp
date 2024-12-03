import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 

plt.rcParams.update({
    'font.size': 12,          # Default font size
    'axes.titlesize': 14,     # Title font size
    'axes.labelsize': 12,     # Axis label font size
    'xtick.labelsize': 10,    # X-axis tick label font size
    'ytick.labelsize': 10,    # Y-axis tick label font size
    'legend.fontsize': 8,    # Legend font size
    'figure.titlesize': 16    # Figure title font size
})

fig, axs = plt.subplots(2, 2, figsize=(6, 5))

precisions      = ["F64", "F32"]
operators       = ["SRTxx", "TRT04", "TRT12"]
initialisations = ["CEQ", "FEQ", "NEQ", "MEI"]
colours         = ["firebrick", "black", "blue", "green"]

directory = "output"
file_end = "l2error.csv"

for (init_num, init) in enumerate(initialisations):
    path64 = directory + "/" + operators[0] + "_" + init + "_" + precisions[0] + "_" + file_end
    path32 = directory + "/" + operators[0] + "_" + init + "_" + precisions[1] + "_" + file_end
    df64 = pd.read_csv(path64)
    df32 = pd.read_csv(path32)

    # Top left: velocity error.
    u_err64 = np.sqrt(df64["u"]*df64["u"] + df64["v"]*df64["v"])
    u_err32 = np.sqrt(df32["u"]*df32["u"] + df32["v"]*df32["v"])
    axs[0][0].plot(df64["time"][2:], u_err64[2:], color=colours[init_num], label=init)
    axs[0][0].plot(df32["time"][2:], u_err32[2:], color=colours[init_num], linestyle="dotted")
    axs[0][0].set_ylabel(r"$\epsilon (\mathbf{u})$")

    # Top right: pressure error.
    axs[0][1].plot(df64["time"][2:], df64["p"][2:], color=colours[init_num])
    axs[0][1].plot(df32["time"][2:], df32["p"][2:], color=colours[init_num], linestyle="dotted")
    axs[0][1].set_ylabel(r"$\epsilon (p)$")

    # Bottom left: sigma_xx
    axs[1][0].plot(df64["time"][2:], df64["sxx_lb"][2:], color=colours[init_num])
    axs[1][0].plot(df32["time"][2:], df32["sxx_lb"][2:], color=colours[init_num], linestyle="dotted")
    axs[1][0].set_xlabel(r"$t$")
    axs[1][0].set_ylabel(r"$\epsilon (\sigma_{xx})$")

    # Bottom right: sigma_yy
    axs[1][1].plot(df64["time"][2:], df64["sxy_lb"][2:], color=colours[init_num])
    axs[1][1].plot(df32["time"][2:], df32["sxy_lb"][2:], color=colours[init_num], linestyle="dotted")
    axs[1][1].set_xlabel(r"$t$")
    axs[1][1].set_ylabel(r"$\epsilon (\sigma_{xy})$")

for row in range(2):
    for col in range(2):
        axs[row][col].set_xscale("log")
        axs[row][col].set_yscale("log")
        axs[row][col].set_xlim(left=1.0, right=1.0e3)

axs[0][0].legend(loc='lower right', frameon=False)

plt.tight_layout()
plt.savefig("plots/full_error_timeseries.pdf")
plt.show()
