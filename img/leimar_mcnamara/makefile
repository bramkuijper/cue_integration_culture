# makefile to generate python plots
# recapitulating everything done by Leimar & McNamara

SIMFILE_FIG2 = sim_cue_integration_9_12_2019_181703_1298713285
SIMFILE_FIG3 = sim_cue_integration_9_12_2019_201701_846428936


all : fig2.pdf fig3.pdf

fig2.pdf : recapitulate_fig2.py $(SIMFILE_FIG2)
	./recapitulate_fig2.py $(SIMFILE_FIG2) fig2.pdf

fig3.pdf : recapitulate_fig2.py $(SIMFILE_FIG3)
	./recapitulate_fig2.py $(SIMFILE_FIG3) fig3.pdf

clean :
	rm -rf fig2.pdf

.PHONY: clean

