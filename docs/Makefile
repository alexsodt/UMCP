
EXECS = hamiltonian.pdf doc.pdf

all: $(EXECS)

hamiltonian.pdf: hamiltonian_chapter.tex hamiltonian.tex
	pdflatex hamiltonian_chapter
	pdflatex hamiltonian_chapter
	bibtex hamiltonian_chapter
	pdflatex hamiltonian_chapter
	cp hamiltonian_chapter.pdf hamiltonian.pdf

doc.pdf: doc.tex running.tex hamiltonian.tex preamble.tex lipid_composition.tex npt.tex introduction.tex srd.tex dynamics.tex
	pdflatex doc
	pdflatex doc
	bibtex doc
	pdflatex doc

