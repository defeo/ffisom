ffisom.pdf: ffisom.tex defeo.bib #arith_prog.pdf
	latexmk -pdf ffisom.tex

arith_prog.pdf: plots.py
	./plots.py    # WARNING: this takes ~2mins!

all: arith_prog.pdf ffisom.pdf 

.PHONY: all

rains.pdf: rains.tex
	latexmk -pdf rains.tex
