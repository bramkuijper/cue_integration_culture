# makefile to generate python plots
# recapitulating everything done by Leimar & McNamara

SIMFILE_FIG2 = sim_cue_integration_9_12_2019_181703_1298713285


fig2.pdf : recapitulate_fig2.py $(SIMFILE_FIG2)
	./recapitulate_fig2.py $(SIMFILE_FIG2) fig2.pdf

clean :
	rm -rf fig2.pdf

.PHONY: clean

