#include <stdio.h>
#include <flint/nmod_mat.h>

/*---------------------------------------------------------------*/
/* prints s = M                                                  */
/*---------------------------------------------------------------*/
void sage_output_assign_mat(nmod_mat_t M, char *s){
  long n = nmod_mat_nrows(M);
  long m = nmod_mat_ncols(M);
  printf("%s = Matrix([", s);
  long i, j;
  for (i = 0; i < n-1; i++){
    printf("[");
    for (j = 0; j < m-1; j++)
      printf("k(%lu), ", nmod_mat_entry(M, i, j));
    printf("k(%lu)],", nmod_mat_entry(M, i, j));
  }
  printf("[");
  for (j = 0; j < m-1; j++)
    printf("k(%lu), ", nmod_mat_entry(M, i, j));
  printf("k(%lu)]])\n", nmod_mat_entry(M, i, j));
}
