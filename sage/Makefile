all: libjavad_nmod.so libellmul.so cythonized

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
	$(MAKE) -C javad clean
	$(MAKE) -C ellmul clean
	rm -rf cythonized
	rm -rf *.so

.PHONY: clean 
