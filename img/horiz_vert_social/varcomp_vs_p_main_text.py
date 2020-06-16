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


# plot the mean-a values
def mean_a_panel(
        the_fig
        ,row
        ,col
        ,dataset
        ,query_str
        ,trait_selection
        ,xlim=[-0.05,1.05]
        ,ylim=[-10,10]
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
            "mean_aintercept":r"$a_{\text{0}}$",
            "mean_agen":r"$a_{\mathrm{g}}$",
            "mean_ajuv":r"$a_{\mathrm{ind}}$",
            "mean_bmat_envt":r"$m_{\mathrm{e}}$",
            "mean_bmat_phen":r"$m_{\mathrm{m}}$",
            "mean_vc":r"$v_{\mathrm{c}}$",
            "mean_vp":r"$v_{\mathrm{p}}$",
            "mean_hc":r"$h_{\mathrm{c}}$",
            "mean_hp":r"$h_{\mathrm{p}}$",
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


    the_color_map = cm.get_cmap("tab10")
    
    the_axis = the_fig.start_block(
        row=row
        ,col=col)

    subset = subset.sort_values(by="p")

    xval="p"

    # calculate means and sd for each level of p
    aggregate = subset.groupby(by=xval)[selected_traits].agg(["mean","std"])

    # merge sub and top columns that arise when taking aggregates
    # df.columns.ravel()
    # returns a list of tupels of top columns and sub colums,
    # [('top1','sub1'),('top1','sub2'),...,('topn','subm')]
    # which we then concatenate using "_"
    aggregate.columns = ["_".join(x) for x in aggregate.columns.ravel()]

    aggregate.reset_index(level=0)
    aggregate["pcol"] = aggregate.index

    aggregate["1minp"] = 1.0 - aggregate["pcol"]

    aggregate.sort_values(by="1minp")


    # add a horizontal line indicating randomness
    the_axis.add_artist(lines.Line2D(
        xdata=[0.5,0.5]
        ,ydata=[ylim[0],ylim[1]]
        ,linewidth=1
        ,color="#bcbcbc"))

    if row == 0:
        the_axis.text(x=0.47
                ,y=0.7
                ,s=r"Random"
                ,fontsize=14
                ,rotation=90
                ,transform=the_axis.transAxes)


    # loop through the various traits
    for i, trait_i in enumerate(selected_traits):

        current_color = the_color_map.colors[trait_selection[i]]

        the_axis.fill_between(
                x=aggregate["1minp"]
                ,y1=aggregate[trait_i + "_mean"] - aggregate[trait_i + "_std"]
                ,y2=aggregate[trait_i + "_mean"] + aggregate[trait_i + "_std"]
                ,alpha=0.1
                ,color=current_color
                )

        p1 = the_axis.plot(
                        aggregate["1minp"]
                        ,aggregate[trait_i + "_mean"]
                        ,linewidth=1
                        ,marker="o"
                        ,color=current_color
                        ,markerfacecolor=current_color
                        ,markeredgecolor=current_color
                        ,label=list(traits_n_labels.values())[i])

    if legend:
        the_axis.legend()

    #xlim = [ -0.05, 1.05 ]

    # end the figure
    the_fig.end_block(
            the_axis
            ,xlim=xlim
            ,ylim=ylim
            ,y_ticks_minor = 5
            ,x_ticks_minor = 4
            ,x_ticks_major_multiple = 0.2 
            ,y_ticks_major_multiple = 2.5
            ,xticks=row == the_fig.rows - 1
            ,yticks=col==0
            ,title=title
            ,xlabel=""
            ,ylabel=""
            ,loc_title=True
            ,loc_title_pos=[-0.05,1.05]
            )
    
    return(the_axis)
# plot the mean-a values
def eta_panel(
        the_fig
        ,row
        ,col
        ,dataset
        ,query_str
        ,trait_selection
        ,xlim=[-0.05,1.05]
        ,ylim=[-0.05,1.05]
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
            "eta2_ajuv_X_cue_juv_envt_high":r"$a_{\mathrm{ind}}$",
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

    the_color_map = cm.get_cmap("tab10")
    
    the_axis = the_fig.start_block(
        row=row
        ,col=col)

    subset = subset.sort_values(by="p")

    xval="p"

    # calculate means and sd for each level of p
    aggregate = subset.groupby(by=xval)[selected_traits].agg(["mean","std"])

    # merge sub and top columns that arise when taking aggregates
    # df.columns.ravel()
    # returns a list of tupels of top columns and sub colums,
    # [('top1','sub1'),('top1','sub2'),...,('topn','subm')]
    # which we then concatenate using "_"
    aggregate.columns = ["_".join(x) for x in aggregate.columns.ravel()]

    aggregate.reset_index(level=0)
    aggregate["pcol"] = aggregate.index

    aggregate["1minp"] = 1.0 - aggregate["pcol"]

    aggregate.sort_values(by="1minp")


    # add a horizontal line indicating randomness
    the_axis.add_artist(lines.Line2D(
        xdata=[0.5,0.5]
        ,ydata=[ylim[0],ylim[1]]
        ,linewidth=1
        ,color="#bcbcbc"))

    if row == 0:
        the_axis.text(x=0.47
                ,y=0.7
                ,s=r"Random"
                ,fontsize=14
                ,rotation=90
                ,transform=the_axis.transAxes)

        the_axis.text(x=0.25
                ,y=1.0
                ,s=r"Autocorrelation $+$ve"
                ,fontsize=14
                ,horizontalalignment="center"
                ,transform=the_axis.transAxes)

        the_axis.text(x=0.75
                ,y=1.0
                ,s=r"Autocorrelation $-$ve"
                ,fontsize=14
                ,horizontalalignment="center"
                ,transform=the_axis.transAxes)

    # loop through the various traits
    for i, trait_i in enumerate(selected_traits):

        current_color = the_color_map.colors[trait_selection[i]]

        the_axis.fill_between(
                x=aggregate["1minp"]
                ,y1=aggregate[trait_i + "_mean"] - aggregate[trait_i + "_std"]
                ,y2=aggregate[trait_i + "_mean"] + aggregate[trait_i + "_std"]
                ,alpha=0.1
                ,color=current_color
                )

        p1 = the_axis.plot(
                        aggregate["1minp"]
                        ,aggregate[trait_i + "_mean"]
                        ,linewidth=1
                        ,marker="o"
                        ,color=current_color
                        ,markerfacecolor=current_color
                        ,markeredgecolor=current_color
                        ,label=traits_n_labels[trait_i])

    # end the figure
    the_fig.end_block(
            the_axis
            ,xlim=xlim
            ,ylim=ylim
            ,y_ticks_minor = 5
            ,x_ticks_minor = 4
            ,x_ticks_major_multiple = 0.2 
            ,y_ticks_major_multiple = 0.2
            ,xticks=row == the_fig.rows - 1
            ,yticks=col==0
            ,title=title
            ,xlabel=""
            ,ylabel=""
            ,loc_title=True
            ,loc_title_pos=[-0.05,1.05]
            )
    
    return(the_axis)

##### get the data  #####

data_dir = os.path.join(os.path.expanduser("~"),"Projects/cue_integration_culture/data/")

# the data file is obtained from the command line arguments
data_file_name = "summary_single_logistic.csv"

# read in the data
data = pd.read_csv(os.path.join(data_dir,data_file_name)
        ,sep=";")

data_before_migration = "summary_learning_moment.csv"

data_learning_moment = pd.read_csv(os.path.join(data_dir,data_before_migration)
        ,sep=";")

learning_moment2 = "summary_learning_moment2.csv"

data_learning_moment2 = pd.read_csv(os.path.join(data_dir,learning_moment2)
        ,sep=";")


file_full_factorial = "summary_vary_p_social.csv"
data_full_factorial = pd.read_csv(os.path.join(data_dir,file_full_factorial)
        ,sep=";")

##### data selection #####


######## make the eta plot ########

# start the figure
the_fig = multipanel.MultiPanel(
        panel_widths=[1]
        ,panel_heights=[1,1,1,1]
        ,filename="plot_varcomp_vary_horiz.pdf"
        ,hspace=0.3
        ,width=8
        ,height=15
        )

title = "No social learning" 

#query_str = "qmat == 0.5 & qjuv == 1.0 & mu_hp == 0.0 & adult_survival == 1 & juvenile_survival == 0"

query_str = "qmat == 0.5 & qjuv == 1.0 & mu_hp > 0.0 & mu_vp > 0.0 & sd_hc_noise == 1.0 & sd_vc_noise == 1.0 & juvenile_learns_remote_envt == 0 & adult_survival == 1 & juvenile_survival == 0 & envt_change_at_birth == 0"

# make each panel
ax = eta_panel(
        the_fig=the_fig
        ,row=0
        ,col=0
        ,dataset=data_learning_moment
        ,query_str=query_str
        ,trait_selection=[1,2,3,4]
        ,legend=False
        ,title=title)

ax.legend(loc="best", bbox_to_anchor=(1.0,0.5,0.1,0.5))

title = "With social learning" 

query_str = "qmat == 0.5 & qjuv == 1.0 & mu_hp > 0.0 & sd_hc_noise == 0.0 & juvenile_learns_remote_envt == 1 & adult_survival == 1 & juvenile_survival == 0 & envt_change_at_birth == 0"

ax = eta_panel(
        the_fig=the_fig
        ,row=1
        ,col=0
        ,dataset=data_full_factorial
        ,query_str=query_str
        ,trait_selection=[1,2,3,4,5,6,7,8]
        ,legend=False
        ,title=title)

ax.legend(loc="best", bbox_to_anchor=(1.0,0.5,0.1,0.5))

title = "With social learning; noise in individual learning" 
query_str = "qmat == 0.5 & qjuv == 0.5 & mu_hp > 0.0 & sd_hc_noise == 0.0 & juvenile_learns_remote_envt == 0 & adult_survival == 1 & juvenile_survival == 0"

ax = eta_panel(
        the_fig=the_fig
        ,row=2
        ,col=0
        ,dataset=data_learning_moment2
        ,query_str=query_str
        ,trait_selection=[1,2,3,4,5,6,7,8]
        ,legend=False
        ,title=title)

title = "With social learning; noise in individual \& horizontal social learning" 
query_str = "qmat == 0.5 & qjuv == 0.5 & mu_hp > 0.0 & sd_hc_noise == 0.5 & juvenile_learns_remote_envt == 0 & adult_survival == 1 & juvenile_survival == 0"

ax = eta_panel(
        the_fig=the_fig
        ,row=3
        ,col=0
        ,dataset=data_learning_moment2
        ,query_str=query_str
        ,trait_selection=[1,2,3,4,5,6,7,8]
        ,legend=False
        ,title=title)


the_fig.fig.text(
        x=0.01
        ,y=0.5
        ,s=r"Proportion of variance explained in adult phenotype, $\eta^{2}_{\bar{x}_{\mathrm{ad}}}$"
        ,rotation=90
        ,fontsize=18
        ,horizontalalignment="center"
        ,verticalalignment="center"
        ,transform=the_fig.fig.transFigure)

the_fig.fig.text(
        x=0.52
        ,y=0.06
        ,s=r"Probability of environmental change, $1-p$"
        ,fontsize=18
        ,horizontalalignment="center"
        ,verticalalignment="center"
        ,transform=the_fig.fig.transFigure)

the_fig.close(tight=True)









print("Varcomp plot done.")
######## make the slope plot ########

# start the figure
the_fig = multipanel.MultiPanel(
        panel_widths=[1]
        ,panel_heights=[1,1,1,1]
        ,filename="plot_meana_vary_horiz.pdf"
        ,hspace=0.3
        ,width=8
        ,height=15
        )

title = "No social learning" 

#query_str = "qmat == 0.5 & qjuv == 1.0 & mu_hp == 0.0 & adult_survival == 1 & juvenile_survival == 0"

query_str = "qmat == 0.5 & qjuv == 1.0 & mu_hp > 0.0 & mu_vp > 0.0 & sd_hc_noise == 1.0 & sd_vc_noise == 1.0 & juvenile_learns_remote_envt == 0 & adult_survival == 1 & juvenile_survival == 0 & envt_change_at_birth == 0"

# make each panel
ax = mean_a_panel(
        the_fig=the_fig
        ,row=0
        ,col=0
        ,dataset=data_learning_moment
        ,query_str=query_str
        ,trait_selection=[1,2,3,4,5,6,7,8]
        ,legend=False
        ,title=title)

ax.legend(loc="best", bbox_to_anchor=(1.0,0.5,0.1,0.5))

title = "With social learning" 

query_str = "qmat == 0.5 & qjuv == 1.0 & mu_hp > 0.0 & sd_hc_noise == 0.0 & juvenile_learns_remote_envt == 0 & adult_survival == 1 & juvenile_survival == 0"

ax = mean_a_panel(
        the_fig=the_fig
        ,row=1
        ,col=0
        ,dataset=data_learning_moment
        ,query_str=query_str
        ,trait_selection=[1,2,3,4,5,6,7,8]
        ,legend=False
        ,title=title)

#ax.legend(loc="best", bbox_to_anchor=(1.0,0.5,0.1,0.5))

title = "With social learning; noise in individual learning" 
query_str = "qmat == 0.5 & qjuv == 0.5 & mu_hp > 0.0 & sd_hc_noise == 0.0 & juvenile_learns_remote_envt == 0 & adult_survival == 1 & juvenile_survival == 0"

ax = mean_a_panel(
        the_fig=the_fig
        ,row=2
        ,col=0
        ,dataset=data_learning_moment2
        ,query_str=query_str
        ,trait_selection=[1,2,3,4,5,6,7,8]
        ,legend=False
        ,title=title)

title = "With social learning; noise in individual \& horizontal social learning" 
query_str = "qmat == 0.5 & qjuv == 0.5 & mu_hp > 0.0 & sd_hc_noise == 0.5 & juvenile_learns_remote_envt == 0 & adult_survival == 1 & juvenile_survival == 0"

ax = mean_a_panel(
        the_fig=the_fig
        ,row=3
        ,col=0
        ,dataset=data_learning_moment2
        ,query_str=query_str
        ,trait_selection=[1,2,3,4,5,6,7,8]
        ,legend=False
        ,title=title)


the_fig.fig.text(
        x=0.01
        ,y=0.5
        ,s=r"Average sensitivity to each cue"
        ,rotation=90
        ,fontsize=18
        ,horizontalalignment="center"
        ,verticalalignment="center"
        ,transform=the_fig.fig.transFigure)

the_fig.fig.text(
        x=0.52
        ,y=0.06
        ,s=r"Probability of environmental change, $1-p$"
        ,fontsize=18
        ,horizontalalignment="center"
        ,verticalalignment="center"
        ,transform=the_fig.fig.transFigure)

the_fig.close(tight=True)

print("Slopes plot done too.")






