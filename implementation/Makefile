all: cythonized benchmark

benchmark: sage_local libkummer.so
	g++ -I$(SAGE_LOCAL)/include -L$(SAGE_LOCAL)/lib -L. benchmark.cpp -lpari -lkummer -lflint -o $@

sage_local:
	@if [ -z $(SAGE_LOCAL) ]; then echo "make should be run from within a Sage shell" >&2; exit 1; fi

libkummer.so: sage_local
	$(MAKE) -C kummer_c++_flint
	cp kummer_c++_flint/$@ .

libellmul.so: sage_local
	cd ellmul && ./configure --with-gmp=$(SAGE_LOCAL) --with-mpfr=$(SAGE_LOCAL) --with-flint=$(SAGE_LOCAL)
	$(MAKE) -C ellmul
	cp ellmul/$@ .

CYTHON_SOURCES = $(wildcard *.pyx) $(wildcard *.pxd)

cythonized: $(CYTHON_SOURCES) libkummer_nmod.so libellmul.so
	python setup.py build_ext
	cp cythonized/lib.*/*.so ./

clean:
	$(MAKE) -C kummer_c++_flint clean || true
	$(MAKE) -C ellmul clean
	rm -rf cythonized
	rm -rf *.so
	rm -rf *.pyc
	rm benchmark

.PHONY: clean sage_local 
