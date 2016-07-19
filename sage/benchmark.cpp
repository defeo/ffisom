/*
 * Compile doing stg like:
 *   g++ -g -I"$(SAGE_LOCAL)/include" -L"$(SAGE_LOCAL)/lib" -L"." benchmark.cpp -lpari -ljavad_nmod -lflint -o benchmark
 * and run doing stg like:
 *   LD_LIBRARY_PATH="$(SAGE_LOCAL)/local/lib:." ./benchmark 
 */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <pari/pari.h>

#include <flint/nmod_poly.h>
#include "javad/nmod_poly_isom/ff_embedding.h"

int main (int argc, char **argv)
{
 clock_t t;
 unsigned long l;
 long n, v;
 GEN p, f, g, S1, S2;

 nmod_poly_t g_flint, gen1, gen2;
 int i;
 FFEmbedding *ffemb;
 unsigned long lathr, mpthr;

 pari_init(1<<30, 1UL<<20);

 v = 0L;

 //l = unextprime(1UL<<20);
 l = 2053UL;
 n = 79L;
 if (argc == 3)
 {
  l = atoi(argv[1]);
  n = atoi(argv[2]);
 }
 p = stoi(l);
 pari_printf("%lu %ld %Ps\n", l, n, p);

 f = init_Fq(p, n, v);
 pari_printf("%Ps\n", f);

 g = ZX_to_Flx(f, l);
 t = clock();
 Flx_ffintersect(g, g, n, l, &S1, &S2, NULL, NULL);
 t = clock()-t;
 pari_printf("%lf\n", (double) t/CLOCKS_PER_SEC);
 pari_printf("%Ps\n", S1);

 lathr = 0;
 mpthr = 1<<20;
 nmod_poly_init(g_flint, l);
 nmod_poly_init(gen1, l);
 nmod_poly_init(gen2, l);
 for (i = 0; i <= n; i++)
  nmod_poly_set_coeff_ui(g_flint, i, ((long*)g)[i+2]);
 nmod_poly_print_pretty(g_flint, "x");
 printf("\n");
 ffemb = new FFEmbedding(g_flint, g_flint);
 t = clock();
 ffemb->compute_generators(gen1, gen2, 0, 0);
 t = clock()-t;
 printf("%lf\n", (double) t/CLOCKS_PER_SEC);
 t = clock();
 ffemb->compute_generators(gen1, gen2, 0, 1<<20);
 t = clock()-t;
 printf("%lf\n", (double) t/CLOCKS_PER_SEC);
 t = clock();
 ffemb->compute_generators(gen1, gen2, 1<<20, 1<<20);
 t = clock()-t;
 printf("%lf\n", (double) t/CLOCKS_PER_SEC);
 t = clock();
 ffemb->compute_generators(gen1, gen2, 1<<20, 0);
 t = clock()-t;
 printf("%lf\n", (double) t/CLOCKS_PER_SEC);
 nmod_poly_print_pretty(gen1, "x");
 printf("\n");
 pari_close();

 return 0;
}
