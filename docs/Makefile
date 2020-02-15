
EXECS = doc.pdf

all: $(EXECS)

doc.pdf: doc.tex running.tex hamiltonian.tex preamble.tex lipid_composition.tex npt.tex introduction.tex srd.tex dynamics.tex timescales.tex normal_modes.tex sans.tex options.tex
	pdflatex doc
	pdflatex doc
	bibtex doc
	pdflatex doc

