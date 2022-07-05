## Utilities for plotting and analysis scripts

# Packages
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import argparse
import pickle
import py_mesa_reader as mr

# Constants (cgs)
G = 6.6726e-8
kB = 1.380658e-16
c = 2.99792458e10
mp = 1.67e-24
me = 9.11e-28
arad = 7.5657e-15
sigmarad = 5.6703e-05 
day = 24*3600
yr = 3.1536e7
Msun = 1.989e33
Lsun = 3.85e33
Rsun = 6.955e10


# Uniform plotting
mpl.rcParams.update({

    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    #"axes.labelsize": 10,
    #"font.size": 10,
    # Make the legend/label fonts a little smaller
    #"legend.fontsize": 8,
    #"xtick.labelsize": 8,
    #"ytick.labelsize": 8,
    # Non-italic math
    "mathtext.default": "regular",
    # Tick settings
    "xtick.direction" : "in",
    "ytick.direction" : "in",
    "xtick.top" : True,
    "ytick.right" : True,
    # Short dash sign
    "axes.unicode_minus" : True
})

# Functions
def get_isotope_list():
    isotope_filepath = os.environ.get("MESA_DIR")+'/data/chem_data/isotopes.data'
    isotope_list = []
    with open(isotope_filepath) as f:
        next(f)
        for i,line in enumerate(f):
            if i%4==0:
                isotope_list.append(line.split()[0])
    return isotope_list

def get_most_abundant_isotope(data):
    # returns a list containing the most abundant isotope in each grid point
    isotope_list = get_isotope_list()
    isotopes_present = [iso for iso in isotope_list if (data.in_data(iso) and max(data.bulk_data[iso])>1e-3)]
    most = []
    for i in range(len(data.d)):
        abundances = [data.bulk_data[iso][i] for iso in isotopes_present]
        most.append(isotopes_present[abundances.index(max(abundances))])
        #print(data.d[i],most[-1])
    return most

def latexify_iso(iso):
    for i,char in enumerate(iso):
        if char.isdigit():
            letters = iso[0].upper() + iso[1:i]
            numbers = iso[i:]
            break
    return (r'${}^{%s}$%s'%(numbers,letters))

def get_zones(L):
    # given a list L that looks like [2,2,5,5,5,1,3,3,3,6,8,1,1]
    # return all zones of list indices which correspond to the same number
    # return as a dictionary where the key is the number which repeats, the value is 
    # a list with all zones as (beginning,end). Empty list if no zones for that number
    # For this example, return {'2':[(0,1)]}
    dic = {str(i):[] for i in np.unique(L)}
    for i in range(1,len(L)):
        if L[i] == L[i-1]:
            s = str(L[i])
            if len(dic[s])>=1 and dic[s][-1][1] == i-1: # same zone as before so just update final index
                dic[s][-1][1] = i
            else: # create new zone
                dic[s].append([i-1,i])
    return dic