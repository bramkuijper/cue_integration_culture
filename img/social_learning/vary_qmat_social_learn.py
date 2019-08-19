#!/usr/bin/env python3
import multipanel
import pandas as pd
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
         r"\setmainfont{FreeSans}",         # load additional packages
         r"\setmathsfont(Digits,Latin,Greek)[Uppercase=Italic,Lowercase=Italic]{FreeSans}",
#         r"\setmathsfont{[STIXMath-Regular.otf]}",
#         r"\setmainfont{DejaVu Serif}", # serif font via preamble
         ]
}
mpl.rcParams.update(pgf_with_custom_preamble)
mpl.rcParams["axes.labelsize"] = 16
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["axes.labelpad"] = 16
mpl.rcParams["svg.fonttype"] = "none"


##### get the data  #####
data = pd.read_csv("../../data/summary_social_learn_qmat.csv", sep=";")

##### data selection #####

subset = data.query(
        "np == 2 & nc == 0 & qmat > 0.5 & m == 0.1 & qjuv == 0.5 & p == 0.9").copy(deep=True)

######## make the plot ########

the_fig = multipanel.MultiPanel(
        panel_widths=[1]
        ,panel_heights=[1]
        ,filename="social_learn_vary_qmat_np2.svg"
        ,width=5
        ,height=5
        )


the_axis = the_fig.start_block(
    row=0
    ,col=0)

subset = subset.sort_values(by="qmat")

# get unique values for qmat
qmat_u =list(subset["qmat"].unique())
qmat_u.sort()

traits = ["mean_agen","mean_asoc","mean_amat","mean_ajuv"]
labels = [
        r"$a_{\mathrm{gen}}$"
        ,r"$a_{\mathrm{soc}}$"
        ,r"$a_{\mathrm{mat}}$"
        ,r"$a_{\mathrm{juv}}$"
        ]

for i, trait_i in enumerate(traits):
    mean = []
    mean_plus_sd = []
    mean_minus_sd = []
    for qmat_i in qmat_u:
        mean_trait = subset[(subset["qmat"] == qmat_i)][trait_i].mean()
        sd_trait = subset[(subset["qmat"] == qmat_i)][trait_i].std()
        mean += [mean_trait]
        mean_plus_sd += [mean_trait + sd_trait]
        mean_minus_sd += [mean_trait - sd_trait]

    p1 = the_axis.plot(
                    qmat_u
                    ,mean
                    ,linewidth=1
                    ,marker="o"
                    ,label=labels[i])

    the_axis.fill_between(
            x=qmat_u
            ,y1=mean_minus_sd
            ,y2=mean_plus_sd
            ,alpha=0.1
            )

xlim = [ 0.48, 1.02 ]
ylim = [ 0, 2.6 ]

# end the figure
the_fig.end_block(
        the_axis
        ,xlim=xlim
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = 0.1
        ,y_ticks_major_multiple = 0.5
        ,xticks=True
        ,yticks=True
        ,title=""
        ,ylabel="Weight given to cue"
        ,xlabel=r"Accuracy of maternal cue, $q_{\text{mat}}$"
        ,loc_title=False
        )


the_fig.close(tight=True)
