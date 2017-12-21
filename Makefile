.PHONY: run build show slides paper paper-verbose paper-autorebuild

alles: build run show

jupyterlatex:
	jupyter nbconvert Paper.ipynb --to latex --output paper/jupyter-paper.tex

jupytermarkdown:
	jupyter nbconvert Paper.ipynb --to markdown --output paper/paper.md
	pandoc -o paper/jupyter-paper.tex --filter pandoc-minted paper/paper.md

jupyterverbose:
	pdflatex -shell-escape -output-directory=build jupyter-paper.tex

paper:
	cd paper; pdflatex -shell-escape -interaction=batchmode -halt-on-error -output-directory=build paper.tex

paper-verbose:
	cd paper; pdflatex -shell-escape -output-directory=build paper.tex
	cd paper/build; bibtex paper
	cd paper; pdflatex -shell-escape -output-directory=build paper.tex

paper-autorebuild:
	fswatch -0 paper | xargs -0 -n 1 make paper

slides:
	jupyter-nbconvert --to slides Presentation.ipynb --reveal-prefix=reveal.js --output presentation/presentation

build:
	scons

run:
	cd output; ../build/ADV1D -t 800 -s 1000

show:
	#paraview --state=plot.pvsm
	paraview --data=output/adv1d_..vtr
