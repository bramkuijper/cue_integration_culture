#!/usr/bin/env python3

# image to recapitulate result from Leimer & McNamara

import multipanel
import pandas as pd
import numpy as np
import sys
import string
from matplotlib import cm
import matplotlib.patches as mpatches

# stats for kernel density estimatation
from scipy import stats

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

# declare functions


# calculate cue according to eqns 6, 7 in mcnamara & leimar
def calculate_x_off(bmat=0.0, dmat=0.0, hi_envt=False):

    # make data frame with x axis: maternal phenotype, varying from 0 to 1
    d = {"mat_phen":list(np.linspace(0,1,100)),
            "cue":np.NaN}
    df = pd.DataFrame(data=d)

    # define eqns 6,7
    def xoff(row,hi_envt,bmat,dmat):

        mat_phen = row["mat_phen"]

        if hi_envt:
            return(1.0 / (1.0 + np.exp(-bmat * (mat_phen - 0.5) - dmat))) 
        else:
            return(1.0 / (1.0 + np.exp(-bmat * (mat_phen - 0.5) + dmat))) 

    # calculate eqns 6,7
    df["cue"] = df.apply(xoff, axis=1, hi_envt=hi_envt, bmat=bmat, dmat=dmat)

    return(df)


# calculate phenotype according to eqn (5)
def calculate_u_vs_mom(parvals):

    # make dataframe for various maternal cues
    d = {"xmat":list(np.linspace(0,1,100))
            ,"offspring_phenotype":np.NaN}
    df = pd.DataFrame(data=d)

    def calc_u(row, parvals):
        return(1.0 / (1.0 + np.exp(
            -parvals["amat"] * row["xmat"]
            -parvals["agen"] * parvals["xgen"]
            -parvals["ajuv"] * parvals["xjuv"]
            )))

    df["offspring_phenotype"] = df.apply(
            func=calc_u
            ,axis=1
            ,parvals=parvals)

    return(df)

##### get the data of a single simulation which works #####
data_file_name = sys.argv[1]

data = pd.read_csv(data_file_name
        ,nrows=1500
        ,sep=";")

# distribution of datapoints in last generation
dist_file_name = data_file_name + "_dist"
dist = pd.read_csv(dist_file_name
        ,sep=";")


dist = dist.iloc[list(range(0,dist.shape[0],10))]


#### start the figgo
the_fig = multipanel.MultiPanel(
        panel_widths=[1,1]
        ,panel_heights=[1,1]
        ,filename=sys.argv[2]
        ,width=12
        ,height=12
        )

# panel A: maternal phenotype vs cue to offspring

panel_A = the_fig.start_block(
    row=0
    ,col=0)


# calculate logistic eqns based on Leimar & McNamara eqns 6,7
# see function declared above
mean_bmat = data["mean_bmat_phen"].iloc[-1]
mean_dmat = data["mean_bmat_envt"].iloc[-1]

xoff_line_hi = calculate_x_off(
        bmat = mean_bmat
        ,dmat = mean_dmat
        ,hi_envt = True)

xoff_line_low = calculate_x_off(
        bmat = mean_bmat
        ,dmat = mean_dmat
        ,hi_envt = False)

# get the high and the low environmental strategies
dist_sub0 = dist[(dist["envt"] == 0)]
dist_sub1 = dist[(dist["envt"] == 1)]

colors = ["#63c14f","#f9a538"]
line_colors = ["green","#ff8b00"] 
panel_A.plot(
        dist_sub0["phen_mat"]
        ,dist_sub0["xmat"]
        ,linestyle=""
        ,markersize=1.0
        ,marker="."
        ,markerfacecolor=colors[0]
        ,markeredgecolor=colors[0]
        ,label="Low")

panel_A.plot(
        dist_sub1["phen_mat"]
        ,dist_sub1["xmat"]
        ,linestyle=""
        ,markersize=1.0
        ,marker="."
        ,markerfacecolor=colors[1]
        ,markeredgecolor=colors[1]
        ,label="High")


panel_A.plot(
        xoff_line_low["mat_phen"]
        ,xoff_line_low["cue"]
        ,color=line_colors[0]
        ,linewidth=1)

panel_A.plot(
        xoff_line_hi["mat_phen"]
        ,xoff_line_hi["cue"]
        ,color=line_colors[1]
        ,linewidth=1)


panel_A.legend()



xlim = [ -0.05,1.05]
ylim = [ -0.05,1.05]

# end the figure
the_fig.end_block(
        panel_A
        ,xlim=xlim
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = 0.2
        ,y_ticks_major_multiple = 0.2
        ,xticks=True
        ,yticks=True
        ,title=""
        ,xlabel=r"Maternal phenotype, $u_{\mathrm{mat}}$"
        ,ylabel=r"Cue to offspring, $x_{\mathrm{off}}$"
        ,loc_title=True
        )


# panel B: Maternal cue vs offspring phenotype
panel_B = the_fig.start_block(
    row=0
    ,col=1)

print(data.columns.values)

parameter_values = {
        "amat": data["mean_amat"].iloc[-1],
        "agen": data["mean_agen"].iloc[-1],
        "xgen": data["mean_g"].iloc[-1],
        "ajuv": data["mean_ajuv"].iloc[-1]}
parameter_values_lo = parameter_values

parameter_values_lo["xjuv"] = 0

xmat_vs_phen_mean_low = calculate_u_vs_mom(parameter_values_lo)

color="#a5afe4"
panel_B.plot(
        dist_sub0["xmat"]
        ,dist_sub0["phen_ad"]
        ,linestyle=""
        ,markersize=1.0
        ,marker="."
        ,markerfacecolor=color
        ,markeredgecolor=color
        ,label="Low")

color="#6675b3"
panel_B.plot(
        xmat_vs_phen_mean_low["xmat"]
        ,xmat_vs_phen_mean_low["offspring_phenotype"]
        ,color=color)

the_fig.end_block(
        panel_B
        ,xlim=xlim
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = 0.2
        ,y_ticks_major_multiple = 0.2
        ,xticks=True
        ,yticks=True
        ,title=""
        ,xlabel=r"Maternal cue, $x_{\mathrm{mat}}$"
        ,ylabel=r"Offspring phenotype, $u_{\mathrm{off}}$"
        ,loc_title=True
        )

# panel C: Maternal cue vs offspring phenotype, hi envt
panel_C = the_fig.start_block(
    row=1
    ,col=0)

parameter_values_hi = parameter_values
parameter_values_hi["xjuv"] = 1
xmat_vs_phen_mean_hi= calculate_u_vs_mom(parameter_values_hi)

color="#f565a8"
panel_C.plot(
        dist_sub1["xmat"]
        ,dist_sub1["phen_ad"]
        ,linestyle=""
        ,markersize=1.0
        ,marker="."
        ,markerfacecolor=color
        ,markeredgecolor=color
        ,label="High")

color="#d2417d"
panel_C.plot(
        xmat_vs_phen_mean_hi["xmat"]
        ,xmat_vs_phen_mean_hi["offspring_phenotype"]
        ,color=color)

the_fig.end_block(
        panel_C
        ,xlim=xlim
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = 0.2
        ,y_ticks_major_multiple = 0.2
        ,xticks=True
        ,yticks=True
        ,title=""
        ,xlabel=r"Maternal cue, $x_{\mathrm{mat}}$"
        ,ylabel=r"Offspring phenotype, $u_{\mathrm{off}}$"
        ,loc_title=True
        )

# panel D: Maternal cue vs offspring phenotype, hi envt
panel_D = the_fig.start_block(
    row=1
    ,col=1)

color = {"low": "#707bc4", "high":"#da4c7f"}

kernel_lo = stats.gaussian_kde(dist_sub0["phen_ad"])
kernel_hi = stats.gaussian_kde(dist_sub1["phen_ad"])

xvals = np.linspace(0,1,200)

panel_D.plot(
        xvals
        ,kernel_lo.evaluate(xvals)
        ,color=color["low"]
        ,label="Low")

panel_D.plot(
        xvals
        ,kernel_hi.evaluate(xvals)
        ,color=color["high"]
        ,label="High")

ylim = [ 0, 10 ]

the_fig.end_block(
        panel_D
        ,xlim=xlim
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = 0.2
        ,y_ticks_major_multiple = 2
        ,xticks=True
        ,yticks=True
        ,title=""
        ,xlabel=r"Maternal cue, $x_{\mathrm{mat}}$"
        ,ylabel=r"Offspring phenotype, $u_{\mathrm{off}}$"
        ,loc_title=True
        )

the_fig.close(tight=True)

