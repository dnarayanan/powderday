all: 	powderday_manual.pdf

powderday_manual.pdf: powderday_manual.tex classnotes.sty 
	pdflatex powderday_manual.tex powderday_manual.pdf
	sleep 0.5
	pdflatex powderday_manual.tex powderday_manual.pdf

clean:
	rm *.pdf
	rm *.out *.log *.aux
