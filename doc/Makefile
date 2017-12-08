
.PHONY: output presentation autorebuild

all: presentation

presentation:
	pdflatex -shell-escape -interaction=batchmode -halt-on-error -output-directory=build presentation/presentation.tex

verbose:
	pdflatex -shell-escape -output-directory=build presentation/presentation.tex

autorebuild:
	fswatch -0 presentation | xargs -0 -n 1 make presentation

