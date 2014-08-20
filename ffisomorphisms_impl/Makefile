SUBDIRS=fmpz_mod_poly_isom nmod_poly_isom
GCC = g++
GCCFLAGS="-O3 -I$(SAGE_LOCAL)/include"
LIBS = -lflint -lgmp -lmpfr
export

all: libjavad_nmod.so libjavad_fmpz_mod.so

clean: 
	$(foreach dir, $(SUBDIRS), $(MAKE) -C $(dir) clean;)
	rm libjavad_nmod.so libjavad_fmpz_mod.so

libjavad_nmod.so:
	$(MAKE) -C nmod_poly_isom GCCFLAGS=$(GCCFLAGS)
	$(GCC) -shared -L$(SAGE_LOCAL)/lib $(LIBS) $(wildcard nmod_poly_isom/build/*.o) -o $@

libjavad_fmpz_mod.so:
	$(MAKE) -C fmpz_mod_poly_isom GCCFLAGS=$(GCCFLAGS)
	$(GCC) -shared -L$(SAGE_LOCAL)/lib $(LIBS) $(wildcard fmpz_mod_poly_isom/build/*.o) -o $@
