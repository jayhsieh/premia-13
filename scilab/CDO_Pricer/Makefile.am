AUTOMAKE_OPTIONS = foreign
SUBDIRS = src demo test

doc: $(DISTFILES) webdoc/Doxyfile
	doxygen webdoc/Doxyfile

clean-tmp::
	@for i in nsp2 scilab ;\
	do \
		echo "Cleaning .losses files in $$i" ;\
		(cd $$i; rm -f .losses*) ;\
	done

