SUBDIRS=fmpz_mod_poly_isom nmod_poly_isom
GCC = g++
GCCFLAGS="-O3 -I$(INC_DIR)/include"
export

all: libjavad_nmod.so libjavad_fmpz_mod.so

clean: 
	$(foreach dir, $(SUBDIRS), $(MAKE) -C $(dir) clean;)

libjavad_nmod.so:
	$(MAKE) -C nmod_poly_isom GCCFLAGS=$(GCCFLAGS) LIBRARY=$@

libjavad_fmpz_mod.so:
	$(MAKE) -C fmpz_mod_poly_isom GCCFLAGS=$(GCCFLAGS) LIBRARY=$@
