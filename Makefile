# Build the manuscript + supplement PDFs (pdflatex + bibtex).
.PHONY: all main supp clean

all: main supp

main: main.pdf
supp: supplement.pdf

main.pdf: main.tex references.bib
	pdflatex -interaction=nonstopmode -halt-on-error main.tex
	bibtex main
	pdflatex -interaction=nonstopmode -halt-on-error main.tex
	pdflatex -interaction=nonstopmode -halt-on-error main.tex

supplement.pdf: supplement.tex references.bib software_table.tex
	pdflatex -interaction=nonstopmode -halt-on-error supplement.tex
	bibtex supplement
	pdflatex -interaction=nonstopmode -halt-on-error supplement.tex
	pdflatex -interaction=nonstopmode -halt-on-error supplement.tex

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.lof *.lot *.fls \
	      *.fdb_latexmk *.synctex.gz

distclean: clean
	rm -f main.pdf supplement.pdf
