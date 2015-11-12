/* 
 * File:   ff_embedding.h
 * Author: javad
 *
 * Created on November 11, 2015, 3:23 PM
 */

#ifndef FF_EMBEDDING_H
#define FF_EMBEDDING_H

#include <flint/nmod_poly.h>

class FFEmbedding {
    
    nmod_poly_t modulus1;
    nmod_poly_t modulus2;
    nmod_poly_t x_image;
    
    void compute_trace(nmod_poly_t alpha, nmod_poly_t xi, const nmod_poly_t alpha_init, 
        const nmod_poly_t xi_init, const nmod_poly_t modulus, slong i);
    void compute_xi_init(nmod_poly_t xi_init, const nmod_poly_t modulus, slong r);
    void find_subfield(nmod_poly_t subfield_modulus, nmod_poly_t embedding_image, 
        const nmod_poly_t modulus, slong degree);
    
public:
    
    /**
     * @param f1 Defining modulus of the first extension
     * @param f2 Defining modulus of the second extension
     */
    FFEmbedding(const nmod_poly_t f1, const nmod_poly_t f2);
    ~FFEmbedding();

    /**
     * Computes generators g1 of ctx_1, and g2 of ctx_2 such that
     * h: ctx_1 --> ctx_2
     *       g1 --> g2
     * is an embedding  
     * @param g1
     * @param g2
     */
    void compute_generators(nmod_poly_t g1, nmod_poly_t g2);
    
    /**
     * Given generators g1 of ctx_1, and g2 of ctx_2 such that
     * h: ctx_1 --> ctx_2
     *	 g1 --> g2
     * is an embedding, this method builds and embedding
     * h: ctx_1 --> ctx_2
     *        x --> g
     * for some g.
     * @param f1
     * @param f2
     */
    void build_embedding(const nmod_poly_t g1, const nmod_poly_t g2);

    /**
     * @param x_image The image of x under the embedding
     */
    void get_x_image(nmod_poly_t x_image);

    /**
     * Computes the image of {@code f} under the embedding using modular composition.
     * @param image
     * @param f
     */
    void compute_image(nmod_poly_t image, const nmod_poly_t f);
};


#endif /* FF_EMBEDDING_H */

