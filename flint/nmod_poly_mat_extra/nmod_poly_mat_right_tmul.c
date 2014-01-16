#include <flint/nmod_poly_mat.h>

void nmod_poly_mat_right_tmul(nmod_poly_mat_t C, const nmod_poly_mat_t L, const nmod_poly_mat_t Bin, long m, long k){
  long qq = nmod_poly_mat_nrows(Bin);
  long rr = nmod_poly_mat_ncols(Bin);
  long pp = nmod_poly_mat_nrows(L);

  nmod_poly_mat_t B;
  nmod_poly_mat_init(B, rr, qq, nmod_poly_mat_modulus(Bin));

  long i, j;

  for (i = 0; i < rr; i++)
    for (j = 0; j < qq; j++)
      nmod_poly_reverse(nmod_poly_mat_entry(B, i, j), nmod_poly_mat_entry(Bin, j, i), m+1);

  nmod_poly_mat_mul(C, L, B);

  for (i = 0; i < pp; i++)
    for (j = 0; j < qq; j++){
      nmod_poly_truncate(nmod_poly_mat_entry(C, i, j), k+m);
      nmod_poly_shift_right(nmod_poly_mat_entry(C, i, j), nmod_poly_mat_entry(C, i, j), m);
    }
	
  nmod_poly_mat_clear(B);
}
