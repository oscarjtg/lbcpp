import pandas as pd
import matplotlib.pyplot as plt

# Read the data from the CSV file
df = pd.read_csv('output/0001_l2error.csv')

# Create a single figure with three subplots side by side
fig, axs = plt.subplots(1, 3, figsize=(18, 5))

# Plot r vs time
axs[0].plot(df['time'], df['r'], label='r vs time')
axs[0].set_xlabel('Time')
axs[0].set_ylabel('r')
axs[0].set_title('r vs Time')
axs[0].set_xscale("log")
axs[0].set_yscale("log")
axs[0].legend()
axs[0].grid(False)

# Plot u vs time
axs[1].plot(df['time'], df['u'], label='u vs time', color='orange')
axs[1].set_xlabel('Time')
axs[1].set_ylabel('u')
axs[1].set_title('u vs Time')
axs[1].set_xscale("log")
axs[1].set_yscale("log")
axs[1].legend()
axs[1].grid(False)

# Plot v vs time
axs[2].plot(df['time'], df['v'], label='v vs time', color='green')
axs[2].set_xlabel('Time')
axs[2].set_ylabel('v')
axs[2].set_title('v vs Time')
axs[2].set_xscale("log")
axs[2].set_yscale("log")
axs[2].legend()
axs[2].grid(False)

# Adjust layout
plt.tight_layout()

# Show the plot
plt.show()