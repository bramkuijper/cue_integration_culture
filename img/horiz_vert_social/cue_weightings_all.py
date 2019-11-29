#!/usr/bin/env python3
import multipanel
import itertools
import pandas as pd
import sys
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
         r"\setmainfont{[MyriadPro-Regular.otf]}",         # load additional packages
         r"\setmathsfont(Digits,Greek)[Uppercase=Italic,Lowercase=Italic]{[MyriadPro-Regular.otf]}",
         r"\setmathsfont(Latin)[Uppercase=Italic,Lowercase=Italic]{[MyriadPro-It.otf]}",
         r"\setmathrm{[MyriadPro-Regular.otf]}",
#         r"\setmainfont{DejaVu Serif}", # serif font via preamble
         ]
}
mpl.rcParams.update(pgf_with_custom_preamble)
mpl.rcParams["axes.labelsize"] = 18
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["ytick.labelsize"] = 16
mpl.rcParams["xtick.labelsize"] = 16
mpl.rcParams["axes.labelpad"] = 16
mpl.rcParams["svg.fonttype"] = "none"


##### get the data  #####
data = pd.read_csv("../../data/summary_cue_int_big.csv", sep=";")
data_juv_surv = pd.read_csv("../../data/summary_cue_int_big_surv1.csv", sep=";")
data = data.append(data_juv_surv)

data = pd.read_csv("variance_components.csv",sep=";")


##### data selection #####
subset = data.query(
        "p > 0 & p < 1.0 & qjuv == 1.0 & qmat == 1 & sdmat == 0.05 & sdsoc_horiz == 0.05" +
        " & sdsoc_vert == 0.05 & juvenile_survival == 0").copy(deep=True)

print(subset.shape)
######## make the plot ########

traits_n_labels = {
        "var_component_gen_prop":r"$a_{\text{gen}}$",
        "var_component_ajuv_prop":r"$a_{\mathrm{juv}}$",
        "var_component_amat_prop":r"$a_{\mathrm{mat}}$",
        "var_component_ahoriz_prop":r"$a_{\mathrm{horizontal}}$",
        "var_component_avert_prop":r"$a_{\mathrm{vertical}}$"
        }        

trait_selection = [ 3, 4 ]

selected_traits = [ list(traits_n_labels.keys())[i] for i in trait_selection]

print(selected_traits)

selected_traits_name = "_".join(selected_traits)

the_color_map = cm.get_cmap("tab10")

the_fig = multipanel.MultiPanel(
        panel_widths=[1]
        ,panel_heights=[1]
        ,filename="cue_strength_with_p" + selected_traits_name + ".svg"
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

xlim = [ 0, 1.0 ]
ylim = [ -0.1, 1.0 ]

# end the figure
the_fig.end_block(
        the_axis
        ,xlim=xlim
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = 0.1
        ,y_ticks_major_multiple = 0.2
        ,xticks=True
        ,yticks=True
        ,title=""
        ,ylabel="Proportion of phenotypic variance"
        ,xlabel=r"Rate of environmental change"
        ,loc_title=False
        )


the_fig.close(tight=True)
