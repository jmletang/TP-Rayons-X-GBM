pdflatex exo1.tex
pdfcrop exo1.pdf exo1.pdf
convert -density 200 exo1.pdf -resize 80% -quality 100 exo1.png
#pdf2svg exo1.pdf exo1.svg
rm exo1.aux exo1.log exo1.pdf
