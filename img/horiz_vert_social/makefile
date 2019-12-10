# makefile to plot figures containing relative variance components

# data file with the original analysis of variance components
DATAVARCOMP=../../data/summary_cue_int_finegrained_p.csv

# data file where only maternal effects are evolvingwhere only maternal effects are evolvingwhere only maternal effects are evolvingwhere only maternal effects are evolving
DATAVARCOMPMAT=../../data/summary_cue_int_finegrained_p_mateffect_only.csv

all : variance_components_full.pdf variance_components_only_mateffect_qjuv1_qmat1.pdf variance_components_only_mateffect_qjuv0.5_qmat1.pdf


variance_components_full.pdf : cue_weightings_all.py $(DATAVARCOMP)
	./cue_weightings_all.py -o variance_components_full.pdf -i $(DATAVARCOMP) --qjuv 1.0 --qmat 1.0 --juvsurv 0

variance_components_only_mateffect_qjuv1_qmat1.pdf : cue_weightings_all.py $(DATAVARCOMPMAT)
	./cue_weightings_all.py -o variance_components_only_mateffect_qjuv1_qmat1.pdf -i $(DATAVARCOMPMAT) --qjuv 1.0 --qmat 1.0 --juvsurv 0

variance_components_only_mateffect_qjuv0.5_qmat1.pdf : cue_weightings_all.py $(DATAVARCOMPMAT)
	./cue_weightings_all.py -o variance_components_only_mateffect_qjuv0.5_qmat1.pdf -i $(DATAVARCOMPMAT) --qjuv 0.5 --qmat 1.0 --juvsurv 0

.PHONY: clean


