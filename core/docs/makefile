# simple latex compilation makefile - aslak stubsgaard

# Source files
sources = $(wildcard *.tex)
out = $(sources:.tex=.pdf)

# PdfLaTeX compilation options
latexopt = -halt-on-error -file-line-error

.PHONY: all
all: $(out)

%.pdf: %.tex %.bib
	pdflatex $(latexopt) $<
	biber $*
	pdflatex $(latexopt) $<
	pdflatex $(latexopt) $<
	make clean

.PHONY: clean
clean:
	rm --force *.log *.aux *.out *.bbl *.bcf *.blg *.xml *.fdb_latexmk *.fls *.toc
