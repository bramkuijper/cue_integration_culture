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
# read in the final-generation 
# trait distribution data
#########################################

dist_dat = pd.read_csv(
        sys.argv[1].strip() + "_dist",
        sep=";")

distribution_available = dist_dat.shape[0] > 0

#########################################

#           make figure

#########################################

# initialize the figure
fig = plt.figure(figsize=(10,40),dpi=200)

nrows = 7

if distribution_available:
    nrows += 4

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
        ,dat["mean_ad_phen"]
        ,label=r"$u$")

ax.set_ylabel(r"Ad phenotype, $a$")

rowctr +=1

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



rowctr +=1

ax = plt.subplot(gs[rowctr,0])

ax.plot(
        dat["generation"]
        ,dat["mean_surv0"]
        ,label=r"Surv $e_{\mathrm{low}}$")

ax.plot(
        dat["generation"]
        ,dat["mean_surv1"]
        ,label=r"Surv $e_{\mathrm{high}}$")

ax.set_ylabel(r"Survival")
ax.set_xlabel(r"Generation, $t$")
ax.legend()

rowctr +=1

ax = plt.subplot(gs[rowctr,0])

ax.plot(
        dat["generation"]
        ,dat["freq_high"])

ax.set_ylabel(r"Freq high envt")
ax.set_xlabel(r"Generation, $t$")

if distribution_available:
    rowctr +=1

    ax = plt.subplot(gs[rowctr,0])

    ax.plot(
            dist_dat["ad_mat"]
            ,dist_dat["xmat"]
            ,color="#ff008d"
            ,linestyle=""
            ,markersize=0.5
            ,marker=".")

    ax.set_xlim(0,1)
    ax.set_ylim(0,1)

    ax.set_xlabel(r"Maternal phenotype, $u_{\mathrm{mat}}$")
    ax.set_ylabel(r"Cue to offspring, $x_{\mathrm{mat}}$")
    
    rowctr +=1

    ax = plt.subplot(gs[rowctr,0])

    ax.plot(
            dist_dat["xmat"]
            ,dist_dat["ad_phen"]
            ,color="blue"
            ,linestyle=""
            ,markersize=0.5
            ,marker=".")

    ax.set_xlim(0,1)
    ax.set_ylim(0,1)

    ax.set_xlabel(r"Maternal cue, $x_{\mathrm{mat}}$")
    ax.set_ylabel(r"Phenotype, $u$")
    
    
    rowctr +=1

    ax = plt.subplot(gs[rowctr,0])

    ax.plot(
            dist_dat["envt"]
            ,dist_dat["g"]
            ,color="red"
            ,linestyle=""
            ,markersize=0.5
            ,marker=".")

    ax.set_xlim(-0.1,1.1)
    ax.set_ylim(-3.5,3.5)

    ax.set_xlabel(r"Environment, $e$")
    ax.set_ylabel(r"Genotype, $g$")

    rowctr +=1
    
    ax = plt.subplot(gs[rowctr,0])

    ax.plot(
            dist_dat["envt"]
            ,dist_dat["cue_ad_envt_high"]-0.01
            ,color="red"
            ,linestyle=""
            ,label="Adult cue"
            ,marker=".")
    
    ax.plot(
            dist_dat["envt"]
            ,dist_dat["cue_juv_envt_high"]
            ,color="blue"
            ,linestyle=""
            ,label="Juvenile cue"
            ,marker=".")

    ax.set_xlim(-0.1,1.1)
    ax.set_ylim(-0.1,1.1)
    ax.legend()

    ax.set_xlabel(r"Environment, $e$")
    ax.set_ylabel(r"Genotype, $g$")
format = "jpg"

filename = os.path.join(
        os.path.dirname(sys.argv[1]),
        "graph_" + os.path.basename(sys.argv[1]) + "." + format
        )

plt.savefig(
        fname=filename
        ,format=format
        ,dpi="figure"
        ,bbox_inches="tight")
