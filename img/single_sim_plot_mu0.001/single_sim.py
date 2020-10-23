#!/usr/bin/env python3
from pythontools import multipanel
import itertools
import pandas as pd
import sys
import os.path
import argparse
from matplotlib import cm
import matplotlib.lines as lines
import matplotlib.patches as mpatches

import matplotlib as mpl
mpl.use("pgf")

from matplotlib import rc


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
mpl.rcParams['axes.unicode_minus'] = False
mpl.rcParams["axes.labelsize"] = 18
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["ytick.labelsize"] = 16
mpl.rcParams["xtick.labelsize"] = 16
mpl.rcParams["axes.labelpad"] = 16
mpl.rcParams["svg.fonttype"] = "none"

mpl.rcParams["axes.formatter.use_mathtext"]: True

# plot the aggregate of n simulations
def sim_agg_panel(
        the_fig
        ,row
        ,col
        ,dataset
        ,trait_selection
        ,xlim=[0,300000]
        ,ylim=[-7.5,2.5]
        ,ylab=""
        ,legend=False
        ,title=""
        ,y_ticks_major_multiple=2.5
        ):

    id_u = dataset["id"].unique()
    
    the_axis = the_fig.start_block(
        row=row
        ,col=col)

    the_color_map = cm.get_cmap("tab10")

    for trait_idx, trait in enumerate(trait_selection["identifiers"]):

        current_color = the_color_map.colors[trait_idx]

        for id_idx, id_i in enumerate(id_u):

            label = "_nolegend_"

            if id_idx == 0:
                label = trait_selection["labels"][trait_idx]

            the_axis.plot(
                dataset[dataset["id"] == id_i]["generation"]
                ,dataset[dataset["id"] == id_i][trait]
                ,color=current_color
                ,alpha=0.5
                ,label=label)

    # end the figure
    the_fig.end_block(
            the_axis
            ,xlim=xlim
            ,ylim=ylim
            ,y_ticks_minor = 5
            ,x_ticks_minor = 4
            ,x_ticks_major_multiple = 50000
            ,y_ticks_major_multiple = y_ticks_major_multiple
            ,xticks=row == the_fig.rows - 1
            ,yticks=col==0
            ,title=""
            ,xlabel=""
            ,ylabel=ylab
            ,loc_title=True
            ,loc_title_pos=[-0.05,1.05]
            )

    if legend and type(the_axis) != type(None):
        the_axis.legend(loc="best", bbox_to_anchor=(1.0,0.5,0.1,0.5))

    
    return(the_axis)

def read_in_all_sims(
        file_names
        ,max_rows=1501
        ):

    all_files = None

    for i, file_name in enumerate(file_names):
        df = pd.read_csv(filepath_or_buffer=file_name
                ,sep=";"
                ,nrows=max_rows)

        df = df.iloc[::10,:]

        df["id"] = i

        if type(all_files) == type(None):
            all_files = df
        else:
            all_files = all_files.append(
                    other=df
                    ,ignore_index=True)

    return(all_files)



summary_data = pd.read_csv("summary.csv",sep=";")

summary_data = summary_data.loc[summary_data["m"] == 0.1]

all_data = read_in_all_sims(file_names=list(summary_data["file"]),max_rows=6001)

# start the figure
the_fig = multipanel.MultiPanel(
        panel_widths=[1]
        ,panel_heights=[1,1,1,1,1,1]
        ,filename="simulation_example_mu0.001.pdf"
        ,hspace=0.3
        ,width=8
        ,height=15
        )

sim_agg_panel(
        the_fig=the_fig
        ,row=0
        ,col=0
        ,dataset=all_data
        ,ylab=r"Genetic cue" + "\n" + r"sensitivity, $a_{\mathrm{g}}$"
        ,trait_selection={
            "identifiers":["mean_agen"]
            ,"labels":[r"$a_{\mathrm{g}}$"]})

sim_agg_panel(
        the_fig=the_fig
        ,row=1
        ,col=0
        ,dataset=all_data
        ,ylab=r"Individual learning," + "\n" + r"$a_{\mathrm{ind}}$"
        ,trait_selection={
            "identifiers":["mean_ajuv"]
            ,"labels":[r"$a_{ind}$"]})

sim_agg_panel(
        the_fig=the_fig
        ,row=2
        ,col=0
        ,dataset=all_data
        ,ylab=r"Maternal cue" + "\n" + r"sensitivity, $m$"
        ,legend=True
        ,trait_selection={
            "identifiers":["mean_bmat_phen","mean_bmat_envt"]
            ,"labels":[r"$m_{\mathrm{m}}$",r"$m_{\mathrm{e}}$"]})


sim_agg_panel(
        the_fig=the_fig
        ,row=3
        ,col=0
        ,dataset=all_data
        ,legend=True
        ,ylab=r"Vertical social" + "\n" + r"learning, $v$"
        ,ylim=[-7.5,2.8]
        ,trait_selection={
            "identifiers":["mean_vp","mean_vc"]
            ,"labels":[r"$v_{\mathrm{p}}$",r"$v_{\mathrm{c}}$"]})

sim_agg_panel(
        the_fig=the_fig
        ,row=4
        ,col=0
        ,dataset=all_data
        ,legend=True
        ,ylab=r"Horizontal social" + "\n" + r"learning, $h$"
        ,trait_selection={
            "identifiers":["mean_hp","mean_hc"]
            ,"labels":[r"$h_{\mathrm{p}}$",r"$h_{\mathrm{c}}$"]})

sim_agg_panel(
        the_fig=the_fig
        ,row=5
        ,col=0
        ,dataset=all_data
        ,ylim=[0.5,1]
        ,legend=True
        ,ylab=r"Survival, $S$"
        ,y_ticks_major_multiple=0.1
        ,trait_selection={
            "identifiers":["mean_surv0","mean_surv1"]
            ,"labels":[r"$S\left(u_{\mathrm{low}}\right)$",r"$S\left(u_{\mathrm{high}}\right)$"]})

the_fig.fig.text(
        x=0.52
        ,y=0.06
        ,s=r"Generation, $t$"
        ,fontsize=18
        ,horizontalalignment="center"
        ,verticalalignment="center"
        ,transform=the_fig.fig.transFigure)
the_fig.close(tight=True)



