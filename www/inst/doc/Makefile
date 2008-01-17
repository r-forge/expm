TEXI2DVI = texi2dvi --pdf
RM = rm -rf

all: expm.pdf clean

expm.pdf : expm.tex
	$(TEXI2DVI) expm.tex

clean:
	$(RM) *-[0-9][0-9][0-9].eps *-[0-9][0-9][0-9].pdf *.aux *.log  *.blg \
	*.bbl *.out *.tex *.rel *~ auto Rplots.ps
