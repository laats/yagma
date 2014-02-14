################################################################################
#
# File:         Makefile
# RCS:          $Header: $
# Description:  Gnu make makefile for yagma
# Author:       Staal Vinterbo
# Created:      Sun Jul 31 14:27:49 2011
# Modified:     Fri May 10 03:54:41 2013 (Staal Vinterbo) staal@mats
# Language:     BSDmakefile
# Package:      N/A
# Status:       Experimental
#
# Makefile is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Makefile is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Makefile; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# (c) Copyright 2011, 2013, Staal Vinterbo, all rights reserved.
#
################################################################################

SRCDIR = ./pysource
BINSTALLDIR = /usr/local/bin


PSRC = $(SRCDIR)/*.py yagm.py setup.py
CSS = html4css1.css voidspace.css
uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')

ifeq ($(uname_S),Darwin)
 RST2MAN = rst2man.py
 RST2HTML = rst2html.py
 RST2LATEX = rst2latex.py
else
 RST2MAN = rst2man
 RST2HTML = rst2html
 RST2LATEX = rst2latex
endif

INSTALL = install
PINSTALL = easy_install

OSRC = utils.ml utils.mli tree.ml tree.mli state.ml state.mli ga.ml ga.mli yagmproto.ml
YAGMCMX = utils.cmx tree.cmx ga.cmx state.cmx ga.cmx
YAGMCMO = utils.cmo tree.cmo ga.cmo state.cmo ga.cmo
OCAMLOPT = ocamlopt
OCAMLC = ocamlc
OCAMLDEP = ocamldep

PACKFLAGS = -package batteries,batteries.pa_comprehension.syntax -syntax camlp4o
OFLAGS = -thread $(PACKFLAGS) 
OPTFLAGS = -unsafe -inline 5 -ccopt -O3 
DEBFLAGS = -g
LDFLAGS = -linkpkg

.SUFFIXES: .mli .cmi .cmx .ml .cmo

options:
	@echo 'targets are: builds sdist egg html pdf paper readme yagm install dist love'
	@echo 'to build and install do: make builds; sudo make install'
	@echo 'to see what each target does, do: make -n target'

.ml.cmo:
	ocamlfind $(OCAMLC) $(OFLAGS) $(DEBFLAGS) -c $<

.mli.cmi:
	ocamlfind $(OCAMLC) $(OFLAGS) -c $<

.ml.cmx:
	ocamlfind $(OCAMLOPT) $(OFLAGS) $(OPTFLAGS)  -c $<

yagm: $(YAGMCMX) yagmproto.ml
	ocamlfind $(OCAMLOPT) $(OFLAGS) $(OPTFLAGS) $(LDFLAGS) -o $@ $^

yagm.deb: $(YAGMCMO) yagmproto.ml
	ocamlfind $(OCAMLC) $(OFLAGS) $(DEBFLAGS) $(LDFLAGS) -o $@ $^


builds: egg html yagm

clean: 
	rm -f *.o *.cm[oix] *.aux *.bbl *.out *.log 

superclean: clean
	rm -f yagm yagma-paper.tex yagma.tex yagma.pdf yagma-paper.pdf yagma.html yagmproto

README.txt: $(SRCDIR)/docstring.py
	python $< | $(RST2MAN) | nroff -man -c  | sed -f sedman.sed > $@

yagma.html: $(SRCDIR)/docstring.py html4css1.css voidspace.css
	python $< | $(RST2HTML) --time --stylesheet=`echo $(CSS) | tr ' ' ',' ` > $@

MANIFEST.in: MANIFEST.template Makefile
	cat MANIFEST.template > MANIFEST.in
	for t in $(OSRC); do \
           echo "include $$t" >> MANIFEST.in ; \
        done

sdist: $(PSRC) $(CSS) readme MANIFEST.in
	python setup.py sdist

egg: $(PSRC) $(CSS) README.txt 
	python setup.py bdist_egg

readme: README.txt

html: yagma.html

pdf: yagma.tex
	pdflatex -interaction nonstopmode "\input" yagma.tex	

latex: yagma.tex


yagma.tex:  $(SRCDIR)/docstring.py
	@python $< | \
          $(RST2LATEX) --use-verbatim-when-possible --stylesheet=amsmath.sty,amssymb.sty,geometry.sty,iwona.sty,parskip.sty --documentoptions=letterpaper --use-latex-citations  | sed -e 's/\\\$$/\$$/g;s/\\{/{/g;s/\\}/}/g;s/\\textbackslash{}/\\/g;s/\\_/_/g;s/\\textasciicircum{}/\^/g;s/\\\&/\&/g' > yagma.tex


paper: yagma.tex
	@(echo '\\documentclass[10pt]{article}';\
        awk '/^% generated/,/\\begin\{document\}/' $< ;\
        echo '\\title{Functional color coding applied to graph matching\\\\\\small{DBMI Technical report Oct 2011}}' ; \
        echo '\\author{Staal A. Vinterbo and Jialan Que\\\\Division of Biomedical informatics, UCSD\\\\\\href{mailto:sav@ucsd.edu}{sav@ucsd.edu}, \\href{mailto:jque@ucsd.edu}{jque@ucsd.edu}}';\
        echo '\\maketitle';\
        awk '/^\\textbf{Abstract}/,/\\end\{document\}/' $<) > yagma-paper.tex
	@echo "Created LaTeX file yagma-paper.tex. Compiling..."
	@pdflatex --interaction nonstopmode "\input" yagma-paper.tex
	@pdflatex --interaction nonstopmode "\input" yagma-paper.tex
	@echo "Created files yagma-paper.tex and yagma-paper.pdf."

install.egg:
	$(PINSTALL) `ls -t dist/yagma-*.egg | head -1`

install.ocaml: yagm
	$(INSTALL) yagm $(BINSTALLDIR) 

install: install.egg install.ocaml

dist: html sdist egg pdf paper
	rm -f dist/index.html
	echo "<html><head></head><body><h3>Files:</h3>" > dist/index.html
	for t in `ls -r dist/yagma-*.tar.gz | head -1`; do \
		p=`basename $$t`; \
		echo "<li> <a href=\"$$p\">$$p</a></li>" >> dist/index.html;\
	done
	echo "</body></html>" >> dist/index.html
	tar cvzf yagma.tgz yagma.html yagm.ml yagma.png yagma.pdf yagma-paper.pdf yagma.ogg dist/* 
	rm -f dist/index.html
	echo "Created yagma.tgz!"

love:
	@echo 'your place or mine?'

ptg: dist
	scp $(pname).tgz ptg.ucsd.edu:www/$(pname)
	echo ssh ptg.ucsd.edu "(cd www/$(pname); ls)"

deps:
	ocamlfind $(OCAMLDEP) $(PACKFLAGS) $(OSRC) > yagma.deps

include yagma.deps

