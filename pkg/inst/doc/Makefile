TEXI2DVI = texi2dvi --pdf
RM = rm -rf
R_PROG = R

all: expm.pdf

%.tex: %.Rnw
	${R_PROG} CMD Sweave $<

%.pdf : %.tex %.bib
	$(TEXI2DVI) $<

clean:
	$(RM) *-[0-9][0-9][0-9].eps *-[0-9][0-9][0-9].pdf *.aux *.log  *.blg \
	*.bbl *.out *.tex *.rel *~ auto Rplots.ps
