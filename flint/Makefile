

DIRS = util sage_output nmod_vec_extra nmod_poly_extra embeddings pq_nmod

all:
	$(foreach dir, $(DIRS), $(MAKE) -C $(dir) ;)

clean:
	rm -f build/*.o
	$(foreach dir, $(DIRS), $(MAKE) -C $(dir)/test ;)
