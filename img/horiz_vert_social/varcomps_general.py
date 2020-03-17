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

## set up argument parsing
#parser = argparse.ArgumentParser()
#
## specify input data file name
#parser.add_argument('-i', default="../../data/summary_cue_int_finegrained_p.csv")
## specify output file name
#parser.add_argument('-o', default="output_graph_var_components.pdf")
#
## specify juvenile cue fidelity
#parser.add_argument('--qjuv', type=float, default=0.5)
#
## specify maternal cue fidelity
#parser.add_argument('--qmat', type=float, default=1.0)
#
## specify juvenile survival
#parser.add_argument('--juvsurv', type=int, default=0)
#args = vars(parser.parse_args())



######## make the plot ########

def clamp(val, min, max):

    if val < min:
        return(min)

    if val > max:
        return(max)

    return(val)
        

# calculate variance components as proportions of total
def calc_proportional_var_components(row):

    trait_names = [
            "aintercept"
            ,"agen_X_g"
            ,"bmat_phen_X_phen_mat_error"
            ,"bmat_envt_X_maternal_envt_cue_erorr"
            ]

    # TODO


    # intercept genetic var component
    var_component_intercept_total = row["var_aintercept.1"] + \
        row["cov_agen_asoc_vert"] + \
        row["cov_agen_asoc_horiz"] + \
        row["cov_agen_ajuv"]

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

    var_component_gen_prop = clamp(var_component_gen_total / var_component_total,0,1)
    var_component_ajuv_prop = clamp(var_component_ajuv_total / var_component_total,0,1)
    var_component_ahoriz_prop = clamp(var_component_asoc_horiz_total / var_component_total,0,1)
    var_component_avert_prop = clamp(var_component_asoc_vert_total / var_component_total,0,1)
    var_component_amat_prop = clamp(var_component_amat_total / var_component_total,0,1)

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


class VarCompPlot:

    def __init__(
        self
        ,data_file
        ,output_file
        ,query_str):

        self.output_file = output_file
        self.data_file = data_file
        self.query_str = query_str

        self.data = pd.read_csv(
            filepath_or_buffer=self.output_file
            ,sep=";")

        # we have no more than 10 traits 
        # hence this color map is safe
        self.the_color_map = cm.get_cmap("tab10")

        self.traits_n_labels = {
                "var_component_gen_prop": r"$\mathrm{var}\left(a_{\text{gen}}\right)$",
                "var_component_ajuv_prop":r"$\mathrm{var}\left(a_{\mathrm{juv}}\right)$",
                "var_component_amat_prop":r"$\mathrm{var}\left(a_{\mathrm{mat}}\right)$",
                "var_component_ahoriz_prop":r"$\mathrm{var}\left(a_{\mathrm{horiz}}\right)$",
                "var_component_avert_prop":r"$\mathrm{var}\left(a_{\mathrm{vert}}\right)$"
                }


    def check_keys(self, subset):
        missing_keys = []
        for key in list(self.traits_n_labels.keys()):
            if key not in subset.columns.values:
                missing_keys.append(key)

        if len(missing_keys) > 0:
            print("following required columns are missing from data frame: " + " ".join(missing_keys))
            sys.exit(1)


    def query(self):
        
        retval = self.data.query(self.query_str) 

        assert(retval.shape[0] > 0)

        return(retval)

    def single_fig(
            self
            ,trait_selection):
       
        # start the figure
        self.fig_obj = multipanel.MultiPanel(
                panel_widths=[1]
                ,panel_heights=[1]
                ,filename=args["o"]
                ,width=10
                ,height=5
                )

        self.single_panel(trait_selection)

        the_fig.close(tight=True)

    def single_panel(
            self
            ,trait_selection
            ,row=0
            ,col=0):

        the_axis = self.fig_obj.start_block(
            row=row
            ,col=col)

        subset = self.query()

        self.check_keys(subset)

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

        # all possible traits that can be plotted
        all_possible_traits = list(self.traits_n_labels.keys())

        for trait_i in selected_traits:

            trait_idx = all_possible_traits.index(trait_i)
            trait_label = self.traits_n_labels[trait_i]

            current_color = self.the_color_map.colors[trait_idx]

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
                            ,label=trait_label
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

vcp = VarCompPlot(
        data_file="../../data/summary_single_logistic.csv"
        ,output_file="one_logistic"
        ,query_str = "mu_agen == 0 &" +\
                "mu_hc == 0 & mu_hp == 0 &" +\ 
                "mu_vp == 0 & mu_vc == 0 &" +\
                "mu_bmat_phen == 0 & mu_bmat_envt == 0 &")

vcp.single_fig(
    trait_selection="")
