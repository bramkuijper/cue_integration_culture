#!/usr/bin/env python3
import multipanel
import itertools
import pandas as pd
import sys
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

# specify output file name
parser.add_argument('-i', default="../../data/summary_cue_int_finegrained_p.csv")
parser.add_argument('-o', default="output_graph_var_components.pdf")
parser.add_argument('--qjuv', type=float, default=0.5)
parser.add_argument('--qmat', type=float, default=1.0)
parser.add_argument('--juvsurv', type=int, default=0)
args = vars(parser.parse_args())


##### get the data  #####

data_file_name = args["i"]
#data_file_name = "../../data/summary_cue_int_finegrained_p_mateffect_only.csv"

data = pd.read_csv(data_file_name
        ,sep=";")

##### data selection #####
query_str = "qjuv == " + str(args["qjuv"]) +\
        " & qmat == " + str(args["qmat"]) +\
        " & juvenile_survival == " + str(args["juvsurv"]) +\
        " & sdmat == 0.05 " +\
        " & sdsoc_horiz == 0.05" +\
        " & sdsoc_vert == 0.05"

print(query_str)

subset = data.query(query_str).copy(deep=True)
assert(subset.shape[0] > 0)

######## make the plot ########

# calculate variance components as proportions of total
def calc_proportional_var_components(row):

    # total genetic var component
    var_component_gen_total = row["var_component_gen"] + \
        row["cov_agen_asoc_vert"] + \
        row["cov_agen_asoc_horiz"] + \
        row["cov_agen_ajuv"]
   
    # total maternal var component
    var_component_amat_total = row["var_component_amat"] + \
        row["cov_amat_asoc_vert"] + \
        row["cov_amat_asoc_horiz"] + \
        row["cov_amat_ajuv"]

    var_component_ajuv_total =\
            row["var_component_ajuv"] +\
            row["cov_ajuv_asoc_vert"] +\
            row["cov_ajuv_asoc_horiz"] +\
            row["cov_agen_ajuv"] +\
            row["cov_amat_ajuv"]
        
    var_component_asoc_vert_total =\
            row["var_component_asoc_vert"] +\
            row["cov_ajuv_asoc_vert"] +\
            row["cov_agen_asoc_vert"] +\
            row["cov_amat_asoc_vert"]

    var_component_asoc_horiz_total =\
            row["var_component_asoc_horiz"] +\
            row["cov_ajuv_asoc_horiz"] +\
            row["cov_agen_asoc_horiz"] +\
            row["cov_amat_asoc_horiz"]

    var_component_asoc_horiz_c_total =\
            row["var_component_asoc_horiz_c"] +\
            row["cov_ajuv_asoc_horiz_c"] +\
            row["cov_agen_asoc_horiz_c"] +\
            row["cov_amat_asoc_horiz_c"]

    var_component_asoc_horiz_p_total =\
            row["var_component_asoc_horiz_p"] +\
            row["cov_ajuv_asoc_horiz_p"] +\
            row["cov_agen_asoc_horiz_p"] +\
            row["cov_amat_asoc_horiz_p"]
        
    var_component_asoc_vert_c_total =\
            row["var_component_asoc_vert_c"] +\
            row["cov_ajuv_asoc_vert_c"] +\
            row["cov_agen_asoc_vert_c"] +\
            row["cov_amat_asoc_vert_c"]
        
    var_component_asoc_vert_p_total =\
            row["var_component_asoc_vert_p"] +\
            row["cov_ajuv_asoc_vert_p"] +\
            row["cov_agen_asoc_vert_p"] +\
            row["cov_amat_asoc_vert_p"]


    var_component_asoc_vert_total =\
            row["var_component_asoc_vert"] +\
            row["cov_ajuv_asoc_vert"] +\
            row["cov_agen_asoc_vert"] +\
            row["cov_amat_asoc_vert"]

    var_component_total =\
            2 * row["cov_agen_asoc_vert"] + \
            2 * row["cov_agen_asoc_horiz"]  + \
            2 * row["cov_agen_asoc_vert"] + \
            2 * row["cov_agen_ajuv"]  + \
            2 * row["cov_ajuv_asoc_vert"]  + \
            2 * row["cov_ajuv_asoc_horiz"] + \
            2 * row["cov_amat_ajuv"] +\
            2 * row["cov_amat_asoc_horiz"] + \
            row["var_component_gen"] + \
            row["var_component_ajuv"] + \
            row["var_component_amat"] + \
            row["var_component_asoc_vert"] + \
            row["var_component_asoc_horiz"]

    var_component_gen_prop = var_component_gen_total / var_component_total
    var_component_ajuv_prop = var_component_ajuv_total / var_component_total
    var_component_ahoriz_prop = var_component_asoc_horiz_total / var_component_total
    var_component_avert_prop = var_component_asoc_vert_total / var_component_total
    var_component_amat_prop = var_component_amat_total / var_component_total

    return(pd.Series(
        {
            "var_component_gen_total":var_component_gen_total,
            "var_component_amat_total":var_component_amat_total,
            "var_component_ajuv_total":var_component_ajuv_total,
            "var_component_asoc_vert_total":var_component_asoc_vert_total,
            "var_component_asoc_horiz_total":var_component_asoc_horiz_total,
            "var_component_asoc_horiz_c_total":var_component_asoc_horiz_c_total,
            "var_component_asoc_horiz_p_total":var_component_asoc_horiz_p_total,
            "var_component_asoc_vert_c_total":var_component_asoc_vert_c_total,
            "var_component_asoc_vert_p_total":var_component_asoc_vert_p_total,
            "var_component_total":var_component_total,
            "var_component_gen_prop":var_component_gen_prop,
            "var_component_ajuv_prop":var_component_ajuv_prop,
            "var_component_ahoriz_prop":var_component_ahoriz_prop,
            "var_component_avert_prop":var_component_avert_prop,
            "var_component_amat_prop":var_component_amat_prop
            }))


subset[["var_component_gen_total",
    "var_component_amat_total",
    "var_component_ajuv_total",
    "var_component_asoc_vert_total",
    "var_component_asoc_horiz_total",
    "var_component_asoc_horiz_c_total",
    "var_component_asoc_horiz_p_total",
    "var_component_asoc_vert_c_total",
    "var_component_asoc_vert_p_total",
    "var_component_total",
    "var_component_gen_prop",
    "var_component_ajuv_prop",
    "var_component_ahoriz_prop",
    "var_component_avert_prop",
    "var_component_amat_prop"
    ]] = subset.apply(calc_proportional_var_components, axis=1)



traits_n_labels = {
        "var_component_gen_prop":r"$a_{\text{gen}}$",
        "var_component_ajuv_prop":r"$a_{\mathrm{juv}}$",
        "var_component_amat_prop":r"$a_{\mathrm{mat}}$",
        "var_component_ahoriz_prop":r"$a_{\mathrm{horizontal}}$",
        "var_component_avert_prop":r"$a_{\mathrm{vertical}}$"
        }

missing_keys = []
for key in list(traits_n_labels.keys()):
    if key not in subset.columns.values:
        missing_keys.append(key)

if len(missing_keys) > 0:
    print("following required columns are missing from data frame: " + " ".join(missing_keys))
    sys.exit(1)

# list with traits to select
trait_selection = [ 0, 1, 2, 3, 4 ]
selected_traits = [ list(traits_n_labels.keys())[i] for i in trait_selection]

selected_traits_name = "_".join(selected_traits)

the_color_map = cm.get_cmap("tab10")

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
        ,xlabel=r"Rate of environmental change $1 - p$"
        ,loc_title=False
        )


the_fig.close(tight=True)
