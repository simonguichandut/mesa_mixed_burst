# The last inlist (fallback) starts by removing the outer envelope
# which is being ejected (overall has positive velocity). This script
# finds the outermost point with acceleration 0 and then modifies the 
# line remove_initial_surface_by_density in the inlist

import sys
from tkinter import W
# import mesa_reader as mr
import py_mesa_reader as mr
import argparse
import numpy as np

## Load data
filename = "models/ns_env_ejected.mod"
try:
    data = mr.MesaData(filename)
    print("Loaded mod file ", filename)
except:
    print("Could not load mod file ", filename)

r = data.R
rho = data.d
T = data.T
v = data.v
dr = data.R[1:]-data.R[:-1]
dr = np.insert(dr, 0, dr[0]) # just duplicate first element so that array is same length as others

# rhop = (rho[2:]-rho[:-2])/(r[2:]-r[:-2])
# rhopp = (rho[2:]+rho[:-2]-2*rho[1:-1])/((r[2:]-r[:-2])/2)**2 
# Tp = (T[2:]-T[:-2])/(r[2:]-r[:-2])
# Tpp = (T[2:]+T[:-2]-2*T[1:-1])/((r[2:]-r[:-2])/2)**2 
# vp = (v[2:]-v[:-2])/(r[2:]-r[:-2])
# vpp = (v[2:]+v[:-2]-2*v[1:-1])/((r[2:]-r[:-2])/2)**2 


## Different ways to find where to cut the data

def zero_acceleration():
    a = (v[2:]-v[:-2])/(r[2:]-r[:-2])
    for ai in a:
        if ai<0: break
    return list(a).index(ai) + 1

def stepsize_drop():
    ddr = (dr[2:]-dr[:-2])/(r[2:]-r[:-2])
    for ddri in ddr:
        if ddri<-10:                    # arbitrary value found by graphing
            break
    return list(ddr).index(ddri) + 10  # bit deeper in to be safe

def other_method():
    pass

methods = {func.__name__ : func for func in (zero_acceleration,stepsize_drop,other_method)}

## Examine graphically
def plot(k_remove):
    import matplotlib.pyplot as plt
    fig,axes = plt.subplots(3,4,sharex=True,figsize=(13,8))
    for ax in axes[-1]:
        ax.set_xlabel("r-R (cm)")
    
    for i,s in enumerate((r'$\rho$',r'$T$',r'$v$',r'$dr$')):
        axes[0][i].set_ylabel(r"%s"%s)
        axes[1][i].set_ylabel(r"$d$%s/$dr$"%s)
        axes[2][i].set_ylabel(r"$d^2$%s/$dr^2$"%s)

    for i,var in enumerate((rho,T,v,dr)):
        # first and second derivatives
        varp = (var[2:]-var[:-2])/(r[2:]-r[:-2])
        varpp = (var[2:]+var[:-2]-2*var[1:-1])/((r[2:]-r[:-2])/2)**2 
        
        if i==0:
            axes[0][0].loglog(r-12e5, rho, 'k.')
        else:
            axes[0][i].semilogx(r-12e5, var, 'k.')
        axes[1][i].semilogx(r[1:-1]-12e5, varp, 'k.')
        axes[2][i].semilogx(r[1:-1]-12e5, varpp, 'k.')

    for ax in [ax for group in axes for ax in group]:
        ax.axvline(r[k_remove]-12e5, ls='--', color='k', lw=0.7)

    plt.tight_layout()
    plt.show()


# Edit inlist file
def change_file(k_remove):
    filename = "inlist_folder/inlist8_fallback"
    tag = "remove_initial_surface_at_cell_k"
    new_file_contents = ""
    with open(filename,'r') as f:
        for line in f:
            if tag in line:
                new_file_contents += "   remove_initial_surface_at_cell_k = " + str(k_remove) + "\n"
            else:
                new_file_contents += line

    with open(filename,'w') as f:
        f.write(new_file_contents)



# Command line call
parser = argparse.ArgumentParser(description="find index k at which to cut the data given some method (default zero acceleration)")
parser.add_argument('-m', '--method', type=str, default='zero_acceleration')
parser.add_argument('-k','--index', type=int, default=None, help="give index k directly (probably using -p to plot)")
parser.add_argument('-p', '--plot', action='store_true', help="plot to show location of k_remove")
args = parser.parse_args()

if __name__ == "__main__":

    if args.index is not None:
        print("Specified index where to remove")
        k_remove = args.index
    else:
        print("Removing with method: ", args.method)
        k_remove = methods[args.method]()

    print(f"Removing at index {k_remove}, r={r[k_remove]/1e5:.5f} km, rho={rho[k_remove]:.3e} g/cm3")

    if args.plot:
        plot(k_remove)
    else:
        change_file(k_remove)