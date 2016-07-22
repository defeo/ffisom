BENCHMARKS := bench-rains.pdf bench-allombert.pdf bench-allombert-vs-ell.pdf bench-all.pdf
BENCHMARKS := $(patsubst %, plots/%, $(BENCHMARKS))

ffisom.pdf: ffisom.tex refs.bib $(BENCHMARKS) #plots/arith_prog.pdf
	latexmk -pdf ffisom.tex

plots/arith_prog.pdf: plots/plots.py
	cd plots && ./plots.py    # WARNING: this takes ~2mins!

$(BENCHMARKS): plots/benchdata/index.html plots/benchmarks.py
	cd plots && ./benchmarks.py

plots/benchdata/index.html: 
	cd plots/benchdata && wget -r -nH --cut-dirs=2 http://perso.telecom-paristech.fr/~flori/ffisom/
