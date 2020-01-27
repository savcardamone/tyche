### Generic LaTeX Makefile
### Does all the bibtex bullshit as well.
###
### AUTHOR: Salvatore Cardamone
###         University of Cambridge, 2017

# File name where \begin{document} ... \end{document} is
NAME = manual

# Make the lot
all:
	@rm -f *.aux *.bbl *.blg 
	@pdflatex $(NAME)
	@if test -f $(NAME).aux &&  test `grep citation $(NAME).aux | wc -l` -ge 1; then bibtex $(NAME); fi
	@pdflatex $(NAME)
	@pdflatex $(NAME)

# Remove all the temporary crap LaTeX dumps out
clean:
	rm -f *.out *.aux *.bbl *.blg *.dvi *.log *.toc *Notes.bib *~ *.backup $(NAME).ps $(NAME).pdf $(NAME).synctex.gz

