# makefile to plot figures containing relative variance components
#
# python script for variance components plotting
VARCOMPPY=./varcomps_general.py
# python script for plotting the a coefficients
ACOEFFPY=./mean_a_plot.py

# python script for plotting the prestige versus conformism coefficients
HVCOEFFPY=./conformism_versus_prestige.py

# data file with the original analysis of variance components
DATAVARCOMP=../../data/summary_cue_int_finegrained_p.csv

# data file where only maternal effects are evolvingwhere only maternal effects are evolvingwhere only maternal effects are evolvingwhere only maternal effects are evolving
DATAVARCOMPMAT=../../data/summary_cue_int_finegrained_p_mateffect_only.csv

EXT=pdf

VARCOMPONENT_GRAPHS=variance_components_full.$(EXT) variance_components_only_mateffect_qjuv1_qmat1.$(EXT) variance_components_only_mateffect_qjuv0.5_qmat1.$(EXT) 

MEANA_GRAPHS=mean_a_plot_full.$(EXT) mean_a_plot_mateffect_qjuv1_qmat1.$(EXT) mean_a_plot_mateffect_qjuv0.5_qmat1.$(EXT) mean_a_plot_full_juvsurv.$(EXT) mean_a_plot_mateffect_qjuv0.5_qmat1_juvsurv.$(EXT) mean_a_plot_mateffect_qjuv1_qmat1_juvsurv.$(EXT)


all : $(VARCOMPONENT_GRAPHS) $(MEANA_GRAPHS) conformism_vs_prestige_plot_full.$(EXT)

# make the conformism versus prestige bias plots
conformism_vs_prestige_plot_full.$(EXT) : $(HVCOEFFPY) $(DATAVARCOMP)
	$(HVCOEFFPY) -o conformism_vs_prestige_plot_full.$(EXT) -i $(DATAVARCOMP) --qjuv 1.0 --qmat 1.0 --juvsurv 0

# make the a coefficient plots without and with juvenile survival
mean_a_plot_full.$(EXT) : $(ACOEFFPY) $(DATAVARCOMP)
	$(ACOEFFPY) -o mean_a_plot_full.$(EXT) -i $(DATAVARCOMP) --qjuv 1.0 --qmat 1.0 --juvsurv 0

mean_a_plot_full_juvsurv.$(EXT) : $(ACOEFFPY) $(DATAVARCOMP)
	$(ACOEFFPY) -o mean_a_plot_full_juvsurv.$(EXT) -i $(DATAVARCOMP) --qjuv 1.0 --qmat 1.0 --juvsurv 1


mean_a_plot_mateffect_qjuv1_qmat1.$(EXT) : $(ACOEFFPY) $(DATAVARCOMPMAT)
	$(ACOEFFPY) -o mean_a_plot_mateffect_qjuv1_qmat1.$(EXT) -i $(DATAVARCOMPMAT) --qjuv 1.0 --qmat 1.0 --juvsurv 0

mean_a_plot_mateffect_qjuv1_qmat1_juvsurv.$(EXT) : $(ACOEFFPY) $(DATAVARCOMPMAT)
	$(ACOEFFPY) -o mean_a_plot_mateffect_qjuv1_qmat1_juvsurv.$(EXT) -i $(DATAVARCOMPMAT) --qjuv 1.0 --qmat 1.0 --juvsurv 1


mean_a_plot_mateffect_qjuv0.5_qmat1.$(EXT) : $(ACOEFFPY) $(DATAVARCOMPMAT)
	$(ACOEFFPY) -o mean_a_plot_mateffect_qjuv0.5_qmat1.$(EXT) -i $(DATAVARCOMPMAT) --qjuv 0.5 --qmat 1.0 --juvsurv 0

mean_a_plot_mateffect_qjuv0.5_qmat1_juvsurv.$(EXT) : $(ACOEFFPY) $(DATAVARCOMPMAT)
	$(ACOEFFPY) -o mean_a_plot_mateffect_qjuv0.5_qmat1_juvsurv.$(EXT) -i $(DATAVARCOMPMAT) --qjuv 0.5 --qmat 1.0 --juvsurv 1


# make the various variance components
variance_components_full.$(EXT) : $(VARCOMPPY) $(DATAVARCOMP)
	$(VARCOMPPY) -o variance_components_full.$(EXT) -i $(DATAVARCOMP) --qjuv 1.0 --qmat 1.0 --juvsurv 0

variance_components_only_mateffect_qjuv1_qmat1.$(EXT) : $(VARCOMPPY) $(DATAVARCOMPMAT)
	$(VARCOMPPY) -o variance_components_only_mateffect_qjuv1_qmat1.$(EXT) -i $(DATAVARCOMPMAT) --qjuv 1.0 --qmat 1.0 --juvsurv 0

variance_components_only_mateffect_qjuv0.5_qmat1.$(EXT) : $(VARCOMPPY) $(DATAVARCOMPMAT)
	$(VARCOMPPY) -o variance_components_only_mateffect_qjuv0.5_qmat1.$(EXT) -i $(DATAVARCOMPMAT) --qjuv 0.5 --qmat 1.0 --juvsurv 0

.PHONY: clean

clean :
	rm -rf *.$(EXT)
