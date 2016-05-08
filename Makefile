BENCHMARKS := bench-rains.pdf bench-allombert.pdf bench-allombert-vs-ell.pdf bench-all.pdf

ffisom.pdf: ffisom.tex defeo.bib $(BENCHMARKS) #arith_prog.pdf
	latexmk -pdf ffisom.tex

arith_prog.pdf: plots.py
	./plots.py    # WARNING: this takes ~2mins!

benchdata/index.html: 
	cd benchdata && wget -r -nH --cut-dirs=2 http://perso.telecom-paristech.fr/~flori/ffisom/

$(BENCHMARKS): benchdata/index.html benchmarks.py
	./benchmarks.py

all: arith_prog.pdf ffisom.pdf 

.PHONY: all

rains.pdf: rains.tex
	latexmk -pdf rains.tex
