embed.pdf: embed.tex creat.pdf
	pdflatex embed
	bibtex embed
	pdflatex embed
	pdflatex embed

creat.eps: creat.gp bench.dat
	gnuplot creat.gp

creat.pdf: creat.eps
	epstopdf creat.eps
