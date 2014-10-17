all: 	powderday_manual.pdf

powderday_manual.pdf: powderday_manual.tex classnotes.sty 
	pdflatex powderday_manual.tex
	sleep 0.5
	pdflatex powderday_manual.tex

clean:
	rm *.pdf
	rm *.out *.log *.aux
