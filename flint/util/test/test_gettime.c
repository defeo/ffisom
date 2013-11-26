#include <stdio.h>
#include <stdlib.h>
#include <flint/fmpz.h>

#include "util.h"

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
  long i;

  double t = util_gettime();

  // a dummy function call, so it is not optimized away
  for (i = 0; i < 1<<20L; i++)
    util_gettime();

  t = util_gettime() - t;

  printf("%f\n", t);
}


/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if not argument is given, runs timings                     */
/* if the argument 1 is given, runs check                     */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
