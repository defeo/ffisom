#ifndef NMOD_FAST_IRRED_H_
#define NMOD_FAST_IRRED_H_

#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/fq_nmod_poly.h>

class NmodFastIrred {
public:
	void compute_qpower(nmod_poly_t result, const nmod_poly_t g, slong k, slong r, const fq_nmod_ctx_t cyclo_ctx);
	void compute_qmtpower(fq_nmod_t alphaAk, const fq_nmod_t alpha, slong m, slong t, slong r, const fq_nmod_ctx_t cyclo_ctx);
	slong test_residue(const fq_nmod_t alpha, slong r, const fq_nmod_ctx_t cyclo_ctx);
	void compute_qpower_ext(fq_nmod_poly_t result, const fq_nmod_poly_t psi, slong j,
	const fq_nmod_t xiAj, const fq_nmod_t xi, slong n, slong r, const fq_nmod_ctx_t cyclo_ctx);
	void compute_trace_xi(fq_nmod_poly_t trace, const fq_nmod_poly_t psi, slong i, fq_nmod_t xiAi, const fq_nmod_t xiA1, const fq_nmod_t xi, slong n, slong r, const fq_nmod_ctx_t cyclo_ctx);
	void compute_trace(fq_nmod_poly_t trace, const fq_nmod_poly_t psi, const fq_nmod_ctx_t cyclo_ctx, slong r, slong e, const fq_nmod_t xi);
	void irred_prime_power(nmod_poly_t irred, slong r, slong e, mp_limb_t p);
}
;

#endif /* NMOD_FAST_IRRED_H_ */
