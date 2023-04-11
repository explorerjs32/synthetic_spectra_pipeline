import numpy as np
import matplotlib.pyplot as plt
import corner

samples_corner = np.genfromtxt('../results/HD5351_corner_values.txt')
walkers_2d = np.genfromtxt('../results/walkers_HD5351.txt')
walkers_3d = walkers_2d.reshape(walkers_2d.shape[0], walkers_2d.shape[1] // 7, 7)

# Create a corner plot of our samples
labels = labels=['Teff','logg','[Fe/H]','vt','vr','Mg','Na']

fig2 = corner.corner(samples_corner, labels=labels, show_titles=True, plot_datapoints=True, quantiles=[0.16, 0.5, 0.84])
plt.tight_layout()
plt.show()

# Create a plot of the walkers path
fig, axes = plt.subplots(7, figsize=(10, 7), sharex=True)

for i in range(7):
    ax = axes[i]
    ax.plot(walkers_3d[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(walkers_3d))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number");
plt.show()
