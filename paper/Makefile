BENCHMARKS := bench-magma.pdf bench-rains.pdf bench-allombert-lowaux.pdf bench-allombert-anyaux.pdf bench-all.pdf
BENCHMARKS := $(patsubst %, plots/%, $(BENCHMARKS))

ffisom.pdf: ffisom.tex refs.bib $(BENCHMARKS) #plots/arith_prog.pdf
	latexmk -pdf ffisom.tex

plots/arith_prog.pdf: plots/plots.py
	cd plots && ./plots.py    # WARNING: this takes ~2mins!

$(BENCHMARKS): plots/benchmarks.py
	cd plots && ./benchmarks.py
