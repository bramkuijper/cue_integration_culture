#!/usr/bin/env python3

# plot the division of labor simulations

import pandas as pd
import itertools
import subprocess 
import math
import argparse
import numpy as np
import sys, re, os.path
import matplotlib
matplotlib.use('Agg')

from pathlib import Path
from pathlib import PurePath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm


#########################################
#           check where data ends
#########################################
def find_parameter_linenum(filename):

    # get all the lines
    f = open(filename)
    fl = f.readlines()
    f.close()

    # find the first line from below that
    # starts with a number. That is the last
    # line of the data
    frange = list(range(0,len(fl)))
    frange.reverse()

    # ok loop over all lines (backwards)
    for line_num in frange:
        if re.match("^\d.*",fl[line_num]) is not None:
            return(line_num+1)

    return(len(fl))


#########################################
#           read in the data
#########################################

if len(sys.argv) < 2:
    print("plz provide a filename")
    sys.exit()

line_num_params = find_parameter_linenum(sys.argv[1])

dat = pd.read_csv(sys.argv[1],
        nrows=line_num_params-1,
        sep=';',
        index_col=False)

#########################################

#           make figure

#########################################

# initialize the figure
fig = plt.figure(figsize=(10,20))

nrows = 5

nloci_g = 3

# generate the grid of the graph
# see: 
widths = [ 1 ]
heights = [ 1 for x in range(nrows)]
numrows = len(heights)
numcols  = len(widths)

rowctr= 0

# make the grid
gs = gridspec.GridSpec(
        nrows=numrows,
        ncols=numcols,
        width_ratios=widths,
        height_ratios=heights)


# plot the resulting joining probability 
ax = plt.subplot(gs[rowctr,0])


ax.plot(
        dat["generation"]
        ,dat["mean_ajuv"]
        ,label=r"$a_{\mathrm{juv}}$")

ax.plot(
        dat["generation"]
        ,dat["mean_amat"]
        ,label=r"$a_{\mathrm{mat}}$")

ax.plot(
        dat["generation"]
        ,dat["mean_agen"]
        ,label=r"$a_{\mathrm{gen}}$")

ax.set_ylabel(r"Sensitivities, $a$")
ax.legend()

rowctr +=1

# plot the resulting migration probz
ax = plt.subplot(gs[rowctr,0])

ax.plot(
        dat["generation"]
        ,dat["mean_bmat_phen"]
        ,label=r"$b_{\mathrm{phen}}$")

ax.plot(
        dat["generation"]
        ,dat["mean_bmat_envt"]
        ,label=r"$b_{\mathrm{envt}}$")

ax.set_ylabel(r"Maternal sens, $b$")
ax.legend()

rowctr +=1

# plot the resulting migration probz
ax = plt.subplot(gs[rowctr,0])

ax.plot(
        dat["generation"]
        ,dat["mean_g"] - np.sqrt(dat["var_g"]) / nloci_g
        ,label="_nolabel"
        ,color="#96d7ff"
        )

ax.plot(
        dat["generation"]
        ,dat["mean_g"] 
        ,label=r"$\bar{g}$"
        ,color="#009dff"
        )

ax.plot(
        dat["generation"]
        ,dat["mean_g"] + np.sqrt(dat["var_g"]) / nloci_g
        ,label="_nolabel"
        ,color="#96d7ff"
        )

ax.set_ylabel(r"Genetic cue, $g$")
ax.legend()

rowctr +=1

# plot the resulting migration probz
ax = plt.subplot(gs[rowctr,0])

ax.plot(
        dat["generation"]
        ,dat["mean_ad_phen"] - np.sqrt(dat["var_ad_phen"]) 
        ,label="_nolabel"
        ,color="#ffadda")

ax.plot(
        dat["generation"]
        ,dat["mean_ad_phen"] 
        ,label=r"$\bar{u}$"
        ,color="#ff008d"
        )

ax.plot(
        dat["generation"]
        ,dat["mean_ad_phen"] + np.sqrt(dat["var_ad_phen"]) 
        ,label="_nolabel"
        ,color="#ffadda")

ax.set_ylabel(r"Genetic cue, $g$")
ax.legend()

format = "pdf"

filename = os.path.join(
        os.path.dirname(sys.argv[1]),
        "graph_" + os.path.basename(sys.argv[1]) + "." + format
        )

plt.savefig(
        fname=filename,
        format=format, 
        bbox_inches="tight")
