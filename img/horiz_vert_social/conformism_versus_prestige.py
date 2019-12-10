#!/usr/bin/env python3
# plotting relative levels of conformism vs prestige 
# for both horizontal and vertical social learning

import multipanel
import itertools
import pandas as pd
import sys, math
import argparse
from matplotlib import cm
import matplotlib.patches as mpatches

import matplotlib as mpl
mpl.use("pgf")
pgf_with_custom_preamble = {
    "font.family": "serif", # use serif/main font for text elements
    "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": [
        r"\usepackage{units}",         # load additional packages
        r"\usepackage{mathspec}",         # load additional packages
        r"\setmainfont[Path = /usr/share/fonts/personal/ ," +\
            "UprightFont = *-Regular ," +\
            "ItalicFont = *-It ," +\
            "BoldFont = *-Bold ," +\
            "Extension = .otf]{MyriadPro}",
        r"\setmathsfont(Digits,Latin,Greek)[" +\
            "Path = /usr/share/fonts/personal/ ," +\
                        "UprightFont = *-Regular ," +\
                        "ItalicFont = *-It," +\
                        "BoldFont = *-Bold," +\
                        "Extension = .otf]{MyriadPro}",
        r"\setmathrm[" +\
            "Path = /usr/share/fonts/personal/ ," +\
                        "UprightFont = *-Regular ," +\
                        "ItalicFont = *-It," +\
                        "BoldFont = *-Bold," +\
                        "Extension = .otf]{MyriadPro}"
         ]
}
mpl.rcParams.update(pgf_with_custom_preamble)
mpl.rcParams["axes.labelsize"] = 18
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["ytick.labelsize"] = 16
mpl.rcParams["xtick.labelsize"] = 16
mpl.rcParams["axes.labelpad"] = 16
mpl.rcParams["svg.fonttype"] = "none"

# set up argument parsing
parser = argparse.ArgumentParser()

# specify input data file name
parser.add_argument('-i', default="../../data/summary_cue_int_finegrained_p.csv")
# specify output file name
parser.add_argument('-o', default="output_graph_var_components.pdf")

# specify juvenile cue fidelity
parser.add_argument('--qjuv', type=float, default=0.5)

# specify maternal cue fidelity
parser.add_argument('--qmat', type=float, default=1.0)

# specify juvenile survival
parser.add_argument('--juvsurv', type=int, default=0)
args = vars(parser.parse_args())


##### get the data  #####

# the data file is obtained from the command line arguments
data_file_name = args["i"]

# read in the data
data = pd.read_csv(data_file_name
        ,sep=";")

##### data selection #####

# generate the query string dependent on command line args and other things
query_str = "qjuv == " + str(args["qjuv"]) +\
        " & qmat == " + str(args["qmat"]) +\
        " & juvenile_survival == " + str(args["juvsurv"]) +\
        " & sdmat == 0.05 " +\
        " & sdsoc_horiz == 0.05" +\
        " & sdsoc_vert == 0.05"

subset = data.query(query_str).copy(deep=True)

assert(subset.shape[0] > 0)

# get a list of traits
traits_n_labels = {
        "mean_hp":r"$h_{\mathrm{prestige}}$",
        "mean_hc":r"$h_{\mathrm{conformism}}$",
        "mean_vp":r"$v_{\mathrm{prestige}}$",
        "mean_vc":r"$v_{\mathrm{conformism}}$"
        }

missing_keys = []
for key in list(traits_n_labels.keys()):
    if key not in subset.columns.values:
        missing_keys.append(key)

if len(missing_keys) > 0:
    print("following required columns are missing from data frame: " + " ".join(missing_keys))
    sys.exit(1)

# list with traits to select
trait_selection = [ 0, 1, 2, 3 ]
selected_traits = [ list(traits_n_labels.keys())[i] for i in trait_selection]

selected_traits_name = "_".join(selected_traits)

the_color_map = cm.get_cmap("tab10")

# start the figure
the_fig = multipanel.MultiPanel(
        panel_widths=[1]
        ,panel_heights=[1]
        ,filename=args["o"]
        ,width=10
        ,height=5
        )

the_axis = the_fig.start_block(
    row=0
    ,col=0)

subset = subset.sort_values(by="p")

# get unique values for the x axis we are considering here

xval = "p"
x_u = list(subset[xval].unique())
x_u.sort()


# calculate means and sd for each level of p
aggregate = subset.groupby(xval)[selected_traits].agg(
        {'mu':"mean"
        ,'sd':"std"
        })

# df.columns.ravel()
# returns a list of tupels of top columns and sub colums,
# [('top1','sub1'),('top1','sub2'),...,('topn','subm')]

aggregate.columns = ["_".join(x) for x in aggregate.columns.ravel()]

aggregate.reset_index(level=0)
aggregate["pcol"] = aggregate.index

aggregate["1minp"] = 1.0 - aggregate["pcol"]

aggregate.sort_values(by="1minp")

print(aggregate.head())


for i, trait_i in enumerate(selected_traits):

    current_color = the_color_map.colors[trait_selection[i]]

    the_axis.fill_between(
            x=aggregate["1minp"]
            ,y1=aggregate["mu_" + trait_i] - aggregate["sd_" + trait_i]
            ,y2=aggregate["mu_" + trait_i] + aggregate["sd_" + trait_i]
            ,alpha=0.1
            ,color=current_color
            )

    p1 = the_axis.plot(
                    aggregate["1minp"]
                    ,aggregate["mu_" + trait_i]
                    ,linewidth=1
                    ,marker="o"
                    ,color=current_color
                    ,markerfacecolor=current_color
                    ,markeredgecolor=current_color
                    ,label=list(traits_n_labels.values())[i])
the_axis.legend()
xlim = [ -0.05, 1.05 ]
ylim = [ -10, 10.0 ]

# end the figure
the_fig.end_block(
        the_axis
        ,xlim=xlim
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = 0.1
        ,y_ticks_major_multiple = math.floor(abs(ylim[1] - ylim[0]) / 6.0)
        ,xticks=True
        ,yticks=True
        ,title=""
        ,xlabel=r"Rate of environmental change $1 - p$"
        ,ylabel="Social learning coefficients, $v_{i},\ h_{i}$"
        ,loc_title=False
        )


the_fig.close(tight=True)

