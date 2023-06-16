"""
Post-convection hydrogen mass-fraction for all runs listed below
"""

import sys
sys.path.append("../python_scripts")
from utils import *
import seaborn as sns


# Colormap
Cmap = sns.color_palette('colorblind')
# put the second pink and yellow to be at the end so they are used less
Cmap2 = Cmap[:6] + [Cmap[7], Cmap[9], Cmap[6], Cmap[8]]
Cmap = Cmap2

# Runs to include (will use the ns_env_0.90Edd.mod model, or smoothed if
# it exists). Will also plot ns_env_ignited.mod in black
runs = ("Schwarzschild/",
        "Schwarzschild_highres/",
        "Schwarzschild_CG0.1/",
        "Schwarzschild_CG0.1_highres/",
        "Ledoux/",
        "Ledoux_highres/",
        "Ledoux_CG0.1/",
        "Ledoux_CG0.1_highres/",
        "Schwarzschild_PM/",
        "Schwarzschild_PM_CG0.1/",
        "Ledoux_PM/",
        "Ledoux_PM_CG0.1/")

colors = [Cmap[0], Cmap[2], Cmap[1], Cmap[4], Cmap[-1], Cmap[6]]

# String substitutions for run labels
subs = (("/", ""), ("0.1", ""), ("_highres", " (high res.)"),("_", " + "))

def plot_one(mod, color, linestyle, label=None, linewidth=0.9):
    column = integrate.cumtrapz(mod.dq/mod.R**2) * mod.xmstar/(4*np.pi)
    lgy = np.log10(column)
    lgX = np.log10(np.maximum(mod.h1[1:],1e-99))
    ax.plot(lgy, lgX, color=color, lw=linewidth, ls=linestyle, label=label)

# Figure setup
figsize = set_size("mnras1col")
# fig, ax = plt.subplots(1,1, figsize=(figsize[0]/2, figsize[1]/1.5))
fig, ax = plt.subplots(1,1, figsize=figsize)
fig.subplots_adjust(wspace=0.3)
ax.set_xlabel(r"log column depth (g cm$^{-2}$)")
ax.set_ylabel(r"log ${^1}{\rm H}$ mass fraction")
ax.set_xlim([3, 8])
ax.set_xticks(list(range(3, 8)))
ax.set_ylim([-2, 0])
ax.set_yticks([-2, -1.5, -1, -0.5, 0])

# Plot the pre-convection profile (ignition)
mod = mr.MesaData("make_ignition_model/models/ns_env_ignited.mod")
plot_one(mod, color='k', linestyle='-')

for i,run in enumerate(runs):

    # smoothed profiles?
    # if os.path.exists(run+"models/ns_env_0.90Edd_smoothed.mod"):
    #     mod = mr.MesaData(run+"models/ns_env_0.90Edd_smoothed.mod")
    if os.path.exists(run+"models/ns_env_0.90Edd.mod"):
        mod = mr.MesaData(run+"models/ns_env_0.90Edd.mod")
    elif os.path.exists(run+"models/ns_env_Edd.mod"):
        mod = mr.MesaData(run+"models/ns_env_Edd.mod")
    else:
        print("Can't find modfile")
        continue

    color = colors[i//2]
    linestyle = '-' if i%2==0 else '--'

    # highlight some runs?
    # linewidth = 1.5 if i>=8 else 0.9
    linewidth = 0.9

    label = run
    for x in subs:
        label = label.replace(x[0], x[1])

    print(label)

    plot_one(mod, color, linestyle, label, linewidth)

box = (0.5, 1.05)
ax.legend(frameon=False, ncol=int(np.ceil(len(runs)/4)), loc='center', bbox_to_anchor=box, bbox_transform=fig.transFigure, columnspacing=5)
fig.savefig("all_models_h1.pdf", bbox_inches="tight", format="pdf")
