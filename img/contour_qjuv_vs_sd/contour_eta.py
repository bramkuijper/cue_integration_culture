#!/usr/bin/env python3
from pythontools import multipanel
import itertools
import pandas as pd
import sys
import os.path
import argparse
from matplotlib import cm
import matplotlib.colors as mcol
import matplotlib.lines as lines
import matplotlib.patches as mpatches
from matplotlib.colors import BoundaryNorm

import numpy as np

import matplotlib as mpl
mpl.use("pgf")

fontpath = "/System/Library/Fonts/Supplemental/" 

pgf_with_custom_preamble = {
    "font.family": "serif", # use serif/main font for text elements
    "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": [
        r"\usepackage{units}",         # load additional packages
        r"\usepackage{mathspec}",         # load additional packages
        r"\setmainfont[" +\
            "Path = " + fontpath + "," +\
            "UprightFont = * ," +\
            "ItalicFont = *Italic ," +\
            "BoldFont = *Bol," +\
            "Extension = .otf]{STIXGeneral}",
        r"\setmathsfont(Digits,Latin,Greek)[" +\
            "Path = " + fontpath + "," +\
            "UprightFont = * ," +\
            "ItalicFont = *Italic ," +\
            "BoldFont = *Bol," +\
            "Extension = .otf]{STIXGeneral}",
        r"\setmathrm[" +\
            "Path = " + fontpath + "," +\
            "UprightFont = * ," +\
            "ItalicFont = *Italic ," +\
            "BoldFont = *Bol," +\
            "Extension = .otf]{STIXGeneral}",
         ]
}

mpl.rcParams.update(pgf_with_custom_preamble)
mpl.rcParams["axes.labelsize"] = 18
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["ytick.labelsize"] = 16
mpl.rcParams["xtick.labelsize"] = 16
mpl.rcParams["axes.labelpad"] = 16
mpl.rcParams["svg.fonttype"] = "none"

# generate the pivot table that is necessary to plot it in imshow()
def generate_pivot(the_data, x, y, z):

    # make a pivot table
    the_pivot = the_data.pivot_table(
            values=z, 
            index=y, 
            columns=x)

    x, y = np.meshgrid(
            the_pivot.columns.values, 
            the_pivot.index.values)

    z = the_pivot.values

    return(x, y, z)

def find_max_eta(row, traits_n_labels):

    keys = list(traits_n_labels.keys())

    # find the maximum eta value
    max_eta_idx = row[keys].astype("float64").idxmax()

    return(keys.index(max_eta_idx))

# make a contourplot
def contour_panel(
        the_fig
        ,row
        ,col
        ,dataset
        ,x_axis
        ,y_axis
        ,query_str
        ,trait_selection
        ,xlim=[-0.05,1.05]
        ,ylim=[0.5,1.0]
        ,legend=False
        ,title=""
        ):

    # generate the query string dependent on command line args and other things
    subset = dataset.query(query_str).copy(deep=True)

    if subset.shape[0] < 1:
        print("query unsuccessful")
        return

    # list with traits and labels
    traits_n_labels = {
            "eta2_aintercept":r"$a_{\text{0}}$",
            "eta2_agen_X_g":r"$a_{\mathrm{g}}$",
            "eta2_ajuv_X_cue_juv_envt_high":r"$a_{\mathrm{juv}}$",
            "eta2_bmat_envt_X_maternal_envt_cue_eror":r"$m_{\mathrm{e}}$",
            "eta2_bmat_phen_X_phen_mat_error":r"$m_{\mathrm{m}}$",
            "eta2_vc_X_xconformist_vert_error":r"$v_{\mathrm{c}}$",
            "eta2_vp_X_phen_prestige_vert_error":r"$v_{\mathrm{p}}$",
            "eta2_hc_X_xconformist_horiz_error":r"$h_{\mathrm{c}}$",
            "eta2_hp_X_phen_prestige_horiz_error":r"$h_{\mathrm{p}}$",
            }

    missing_keys = []
    for key in list(traits_n_labels.keys()):
        if key not in subset.columns.values:
            missing_keys.append(key)

    if len(missing_keys) > 0:
        print("following required columns are missing from data frame: " + " ".join(missing_keys))
        sys.exit(1)

    # list with traits to select
    trait_list_keys = list(traits_n_labels.keys())
    selected_traits = [ trait_list_keys[i] for i in trait_selection]

    subset["eta2_max"] = subset.apply(func=find_max_eta
            ,axis=1 # apply to each row
            ,args=(traits_n_labels,)
            )

#    subset.to_csv("temp" + "mu_hp_" + str(mu_h_i),sep=";")

    (x,y,z) = generate_pivot(
            the_data=subset
            ,x=x_axis # "1minp"
            ,y=y_axis #"qmat"
            ,z="eta2_max")

    # obtained from Set3 in RColorBrewer
    colors_hex = ["#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9"] 

    the_color_map = mcol.LinearSegmentedColormap.from_list(
            "my_list", colors=colors_hex, N=9)

    the_axis = the_fig.start_block(
        row=row
        ,col=col)

#    levels = [ -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5 ] 
    levels = [ -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5 ] 
    the_color_map = cm.get_cmap("Set3")

    bound_norm = BoundaryNorm(boundaries=levels
            ,ncolors=len(levels)
            ,clip=True)

    the_axis.pcolormesh(
            x,y,z
#            ,origin="lower"
#            ,aspect="auto"
            ,cmap=the_color_map
            ,norm=bound_norm
            ,linewidths=0.01
            ,edgecolors="#ffffff"
#            ,vmin=-1
#            ,vmax=9
            )

    the_fig.end_block(
            the_axis
            ,xlim=None
            ,ylim=None
            ,y_ticks_minor = 5
            ,x_ticks_minor = 4
            ,x_ticks_major_multiple = 0.2
            ,y_ticks_major_multiple = 0.1
            ,xticks=row == the_fig.rows - 1
            ,yticks=col==0
            ,title=title
            ,loc_title_pos=[0.0,1.05]
            ,xlabel=""
            ,ylabel=""
            ,loc_title=True
            )






##### get the data  #####
data_dir = os.path.join(os.path.expanduser("~"),"Projects/cue_integration_culture/data/")

# the data file is obtained from the command line arguments
data_file_name = "summary_contour_vary_sdsoc_vert.csv"

# read in the data
data = pd.read_csv(os.path.join(data_dir,data_file_name)
        ,sep=";")

data["1minp"] = 1 - data["p"]

# start the figure
the_fig = multipanel.MultiPanel(
        panel_widths=[1,1]
        ,panel_heights=[1]
        ,filename="levelplot_max_eta_p_vs_qmat_qjuv.pdf"
        ,hspace=0.3
        ,width=10
        ,height=5
        )

mu_h = list(data["mu_hp"].unique())
mu_h.sort()

row_i = 0
col_i = 0

for mu_h_i in mu_h:

    query_str = "mu_hp == " + str(mu_h_i)

    title = ""

    if row_i == 0:
        title = [r"Without social learning, $\mu_{\mathrm{h}} = \mu_{\mathrm{v}} = 0$"
                ,r"With social learning, $\mu_{\mathrm{h}}, \mu_{\mathrm{v}} > 0$"
                ][col_i]

    # make each panel
    ax = contour_panel(
            the_fig=the_fig
            ,row=row_i
            ,col=col_i
            ,dataset=data
            ,x_axis="1minp"
            ,y_axis="qmat"
            ,query_str=query_str
            ,trait_selection=[0,1,2,3,4] if mu_h_i == 0 else [0,1,2,3,4,5,6,7,8]
            ,legend=row_i == 0
            ,title=title)

    col_i += 1

the_fig.fig.text(
        x=0.01
        ,y=0.5
        ,s=r"Maternal vs juvenile cue" + "\n" + r"fidelity $q_{\mathrm{juv}} = 1.5 - q_{\mathrm{mat}}$"
        ,rotation=90
        ,fontsize=18
        ,horizontalalignment="center"
        ,verticalalignment="center"
        ,transform=the_fig.fig.transFigure)


the_fig.fig.text(
        x=0.52
        ,y=0.045
        ,s=r"Probability of environmental change, $1-p$"
        ,fontsize=18
        ,horizontalalignment="center"
        ,verticalalignment="center"
        ,transform=the_fig.fig.transFigure)


the_fig.close(tight=True)


##### get the data  #####

# the data file is obtained from the command line arguments
#data_file_name = "summary_contour_vary_sdsoc_vert.csv"
data_file_name = "summary_contour_vary_sdsoc_vert.csv"

# read in the data
data = pd.read_csv(os.path.join(data_dir,data_file_name)
        ,sep=";")


# start the figure
the_fig = multipanel.MultiPanel(
        panel_widths=[1,1]
        ,panel_heights=[1]
        ,filename="levelplot_qjuv_vs_sdsoc_vs_p.pdf"
        ,hspace=0.3
        ,width=10
        ,height=5
        )

p = list(data["p"].unique())
p.sort(reverse=True)

row_i = 0
col_i = 0

for p_i in p:

    query_str = "p == " + str(p_i) + " & envt_change_at_birth == 0 & juvenile_learns_remote_envt == 0"

    title = ""

    if row_i == 0:
        title = [r"Positively autocorrelated environment"
                ,r"Negatively autocorrelated environment"
                ][col_i]

    # make each panel
    ax = contour_panel(
            the_fig=the_fig
            ,row=row_i
            ,col=col_i
            ,dataset=data
            ,x_axis="qjuv"
            ,y_axis="sd_vc_noise"
            ,query_str=query_str
            ,trait_selection=[0,1,2,3,4,5,6,7,8]
            ,legend=row_i == 0
            ,title=title)

    col_i += 1

the_fig.fig.text(
        x=0.01
        ,y=0.5
        ,s=r"Maternal vs juvenile cue" + "\n" + r"fidelity $q_{\mathrm{juv}} = 1.5 - q_{\mathrm{mat}}$"
        ,rotation=90
        ,fontsize=18
        ,horizontalalignment="center"
        ,verticalalignment="center"
        ,transform=the_fig.fig.transFigure)


the_fig.fig.text(
        x=0.52
        ,y=0.045
        ,s=r"Probability of environmental change, $1-p$"
        ,fontsize=18
        ,horizontalalignment="center"
        ,verticalalignment="center"
        ,transform=the_fig.fig.transFigure)


the_fig.close(tight=True)

