SUBDIRS=fmpz_mod_poly_isom nmod_poly_isom
GCC = g++
GCCFLAGS="-g -O3 -fPIC -I$(INC_DIR)/include"
export

all: libkummer_nmod.so libkummer_fmpz_mod.so

clean: 
	$(foreach dir, $(SUBDIRS), $(MAKE) -C $(dir) clean;)

libkummer_nmod.so:
	$(MAKE) -C nmod_poly_isom GCCFLAGS=$(GCCFLAGS) LIBRARY=$@

libkummer_fmpz_mod.so:
	$(MAKE) -C fmpz_mod_poly_isom GCCFLAGS=$(GCCFLAGS) LIBRARY=$@
