MANNAME =  PVMOS_manual
MANSRC =  $(MANNAME).tex Octave_PVMOS_Mesh_Topology.pdf
if BUILDDOC
MANPDF=$(MANNAME).pdf
dist_doc_DATA=$(MANPDF) $(MANSRC)  version.tex
PVMOS_manual.pdf: PVMOS_manual.tex Octave_PVMOS_Mesh_Topology.pdf  version.tex
	pdflatex  PVMOS_manual.tex
	pdflatex  PVMOS_manual.tex	
clean:
	rm -f $(MANPDF) *.aux *.log version.tex
else
dist_doc_DATA=$(MANSRC)  version.tex
endif

version.tex:
	echo "\newcommand{\version}{$(VERSION)}" > version.tex
