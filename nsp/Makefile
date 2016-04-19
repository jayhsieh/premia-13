do: compile

compile:
	(cd libpremia; \
	nsp -nw -e 'exec builder.sce; quit'; )


manual: nsp-premia-manual.tex
	latex $<
	makeindex -s index.isty ${<:%.tex=%}
	latex $<
	dvipdf ${<:%.tex=%.dvi}

clean:
	@for i in aux idx ilg ind log out dvi ;\
	do \
		$(RM) nsp-premia-manual.$$i ;\
	done

cleanall: clean
	(cd libpremia; if test -f Makelib; then make -f Makelib clean; fi)
