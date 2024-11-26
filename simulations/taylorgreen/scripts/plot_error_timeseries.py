##################################################
#
# subplot with two columns, three rows
# 
# ---------------------------------------------
#               p_error  | u_error
#   SRTxx
#   TRT04
#   TRT12
#
##################################################


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

precisions      = ["F64", "F32"]
operators       = ["SRTxx", "TRT04", "TRT12"]
initialisations = ["FEQ", "NEQ", "WEI"]
colours         = ["black", "blue", "red"]

fig, axs = plt.subplots(3, 2, figsize=(10, 15))

for row in range(3):
    for (data_num, init) in enumerate(initialisations):
        file64 = operators[row]+"_"+init+"_F64_l2error.csv"
        df64 = pd.read_csv(f"output/{file64}")
        file32 = operators[row]+"_"+init+"_F32_l2error.csv"
        df32 = pd.read_csv(f"output/{file32}")
        print(file64)
        print(file32)

        # First column: error in velocity.
        col = 0
        u_err_64 = np.sqrt(df64["u"]*df64["u"] + df64["v"]*df64["v"])
        u_err_32 = np.sqrt(df32["u"]*df32["u"] + df32["v"]*df32["v"])
        axs[row][col].plot(df64["time"][2:], u_err_64[2:], color=colours[data_num])
        axs[row][col].plot(df32["time"][2:], u_err_32[2:], color=colours[data_num], linestyle="dotted")
        axs[row][col].set_ylabel("error (u)")

        # Second column: error in density.
        col = 1
        axs[row][col].plot(df64["time"][2:], df64["r"][2:], color=colours[data_num])
        axs[row][col].plot(df32["time"][2:], df32["r"][2:], color=colours[data_num], linestyle="dotted")
        axs[row][col].set_ylabel("error (p)")

    for col in range(2):
        axs[row][col].set_xscale("log")
        axs[row][col].set_yscale("log")
        axs[row][col].set_xlim(left=1.0, right=1.0e3)

plt.tight_layout()
plt.show()
fig.savefig("plots/error_timeseries.pdf")
        
