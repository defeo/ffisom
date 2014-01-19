embed.pdf: embed.tex creat.pdf magma.pdf defeo.bib
	pdflatex embed
	bibtex embed
	pdflatex embed
	pdflatex embed

magma.eps: magma.gp test_irred.dat test_ext.dat test_relative.dat test_compositum.dat
	gnuplot magma.gp

magma.pdf:	magma.eps
	epstopdf magma.eps

creat.eps: creat.gp bench.dat
	gnuplot creat.gp

creat.pdf: creat.eps
	epstopdf creat.eps
