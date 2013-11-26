#ifndef EMBEDDINGS_H
#define EMBEDDINGS_H

#include <flint/nmod_poly.h>

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/* Data structures and functions for finite field embeddings:             */
/*   an embeddings_t stores an extension Fq of Fp, using 4 polynomials:   */
/*   * P, such that Fq=Fp[x]/P(x)                                         */
/*   * SP = 1/rev(P) mod x^{deg(P)-1}                                     */
/*   * iP = 1/P' mod P                                                    */
/*   * TP = \sum_{i < deg(P)} trace(z^i mod P) x^i                        */
/*                                                                        */
/* Standard trick: we define an embeddings_t as an array of size 1        */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

typedef struct{
  nmod_poly_t P;
  nmod_poly_t SP;
  nmod_poly_t iP;
  nmod_poly_t TP;
} embeddings_struct;

typedef embeddings_struct embeddings_t[1];

/*------------------------------------------------------------------------*/
/* initializes all polynomials in F                                       */
/*------------------------------------------------------------------------*/
void embeddings_init(embeddings_t F, mp_limb_t n);

/*------------------------------------------------------------------------*/
/* clears all polynomials in F                                            */
/*------------------------------------------------------------------------*/
void embeddings_clear(embeddings_t F);

/*------------------------------------------------------------------------*/
/* assigns all polynomials in F from the defining polynomial P            */
/*------------------------------------------------------------------------*/
void embeddings_set_parameters(embeddings_t F, const nmod_poly_t P);

/*------------------------------------------------------------------------*/
/* computes the composed product R of P and Q                             */
/* assumes that char(Fp) > deg(P) deg(Q) and that Fp is a field           */
/*------------------------------------------------------------------------*/
void embeddings_composed_product_large_char(nmod_poly_t R, const nmod_poly_t P, const nmod_poly_t Q);

/*------------------------------------------------------------------------*/
/* computes the composed product R of P and Q                             */
/* assumes that Fp is a field and that R is irreducible                   */
/*------------------------------------------------------------------------*/
void embeddings_composed_product_small_char(nmod_poly_t R, const nmod_poly_t P, const nmod_poly_t Q);

/*------------------------------------------------------------------------*/
/* computes G = phi(F), where phi is the embedding FP -> FR               */
/*------------------------------------------------------------------------*/
void embeddings_embed(nmod_poly_t G, const nmod_poly_t F, 
		      const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR);

/*------------------------------------------------------------------------*/
/* transpose version of the following:                                    */
/* computes G = phi(F), where phi is the embedding FP -> FR               */
/* Equivalently, computes ell(S^i mod R), for i < m*n, where S=phi(X)     */
/*------------------------------------------------------------------------*/
void embeddings_tembed(mp_ptr ellstar, mp_srcptr ell, 
		       const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR);

/*------------------------------------------------------------------------*/
/* given G = phi(F), where phi is the embedding FP -> FR, computes F      */
/* Equivalently, computes ell(S^i mod R), for i < m*n, where S=phi(X)     */
/*------------------------------------------------------------------------*/
void embeddings_embed_inverse(nmod_poly_t F, const nmod_poly_t G, 
			      const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR);

/*------------------------------------------------------------------------*/
/* computes G = Phi(F), where Phi is the isomorphism FPxFQ -> FR          */
/* F is an array of size m*n (m=deg(P), n=deg(Q))                         */
/* of the form F0..F_{m-1}, with Fi = coeff(F,X^i) \in Fp[y]/Q(y)         */
/*------------------------------------------------------------------------*/
void embeddings_isomorphism(nmod_poly_t G, mp_srcptr F, 
			    const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR);


/*------------------------------------------------------------------------*/
/* input G is a linear form defined modulo R                              */
/* output F is an array of size m*n (m=deg(P), n=deg(Q))                  */
/* of the form F0..F_{m-1}, with Fi = G(S^i T^j), j=0..n-1y]/Q(y)         */
/* with S=Phi(X), T=Phi(Y), Phi = isomorphism FPxFQ -> FR                 */
/*------------------------------------------------------------------------*/
void embeddings_tisomorphism(mp_ptr F, mp_srcptr G,
			     const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR);

/*------------------------------------------------------------------------*/
/* given G, computes F such that G = Phi(F),                              */
/*  where Phi is the isomorphism FPxFQ -> FR                              */
/* F is an array of size m*n (m=deg(P), n=deg(Q))                         */
/* of the form F0..F_{m-1}, with Fi = coeff(F,X^i) \in Fp[y]/Q(y)         */
/*------------------------------------------------------------------------*/
void embeddings_isomorphism_inverse(mp_ptr F, const nmod_poly_t G, 
				    const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR);

#endif
