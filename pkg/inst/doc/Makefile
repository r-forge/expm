TEXI2PDF = texi2dvi --pdf
RM = rm -rf
R_PROG = R

all: expm.pdf

%.tex: %.Rnw
	${R_PROG} CMD Sweave $<

%.pdf : %.tex %.bib
	$(TEXI2PDF) $<

clean:
	$(RM) *-[0-9][0-9][0-9].pdf *.log *.blg *.out *.rel Rplots.ps

realclean: clean
	$(RM) *.tex *.aux *.bbl auto *~
