all: sage_local libjavad_nmod.so libellmul.so cythonized

download:
	cp -r ../ffisom_impl ./javad
	git clone https://github.com/jpflori/ellmul.git ./ellmul

sage_local:
	@if [ -z $(SAGE_LOCAL) ]; then echo "make should be run from within a Sage shell" >&2; exit 1; fi

libjavad_nmod.so:
	$(MAKE) -C javad
	cp javad/nmod_poly_isom/libjavad_nmod.so .

libellmul.so:
	cd ellmul && ./configure --with-gmp=$(SAGE_LOCAL) --with-mpfr=$(SAGE_LOCAL) --with-flint=$(SAGE_LOCAL)
	$(MAKE) -C ellmul
	cp ellmul/libellmul.so .

CYTHON_SOURCES = $(wildcard *.pyx) $(wildcard *.pxd)

cythonized: $(CYTHON_SOURCES)
	python setup.py build_ext
	cp cythonized/lib.*/*.so ./

clean:
	$(MAKE) -C javad clean || true
	$(MAKE) -C ellmul clean
	rm -rf cythonized
	rm -rf *.so
	rm -rf *.pyc

distclean: clean
	rm -rf javad ellmul

.PHONY: clean distclean download sage_local 
