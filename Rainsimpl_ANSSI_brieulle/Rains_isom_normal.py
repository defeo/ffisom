# Every review, remarks and corrections are welcomed  
# and to be sent to l.brieulle(at)gmail(dot)com

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

def convert(z, v, k, normal_basis = None): 
    '''
    INPUT : an element z of K, an normal element of K, a base field k for 
            which K is an extension of degree n.
    OUTPUT : the coefficient of z in the normal basis generated by v.
    
    Algorithm : 

    The goal of this function is to find the coefficients of an element z 
    of a finite field in the normal basis generated by the normal element v
    in the same finite field.
    
    First we have to compute the normal basis generated by v, then we compute 
    the first line of the matrice B = Tr(v^(p^i).v^(p^j)) which a circulant 
    matrix of size n, i.e. :
    
    b_ij = b_kl if and only if i-j = k-l mod n
    
    We now need to compute the vector Z = Tr(v.z^(p^(n-j))); by some 
    calculations detailed in the proposition ? of some paper, we have :
    
    Z = BC,
    
    where C is the vector containing the coefficient of z in v^(p^i). So,
    we need to invert B. The core of the algorithm remains in the fact 
    that the ring of circulant matrix of size n in GF(q) is isomorphic 
    to the ring R = GF(q)[U]/(U^n - 1), via :
    
    (b_ij) -> sum_{i = 0 to n-1}{b_0j.U^i}
    
    Thus, we compute the inverse in R and we "construct" B^(-1) from the 
    coefficient of that element. Each c_i is then equal to :
    
    c_i = d_ij*Z_j,
    
    if we wrote d_ij the coefficient of B^(-1).
     
    '''

    q = k.cardinality()
    n = v.parent().degree()
    R = PolynomialRing(GF(q), 'U')
    U = R.gen()

    if normal_basis is None:
        normal_basis = [v]
        for i in range(n-1):
            normal_basis.append(normal_basis[-1]**p)

    # We need only to compute Tr(v*v^(p^(i))) since the matrix is circulant.
    B = []
    for i in range(n):              
        B.append((v*normal_basis[i]).trace())


    inv = R(B).inverse_mod(U**n - 1)    # We compute the inverse of the image of
                                        # the circulant matrix B in the cyclotomic ring 
                                        # GF(p^n)[U]/(U^n - 1) 
    
    val_trz = []   
    for i in range(n):
        val_trz.append((v*(z**q**(n-i))).trace())
    
    # We will now compute the coefficients c_i while keeping in mind that they 
    # are computed from each rows of the matrix B. So we need to take that into 
    # account and push the coefficient of inv "to the right". Practically, we 
    # just need to index them with (j-i)%n.
    c = []
    for i in range(n):                
        c.append(sum(inv[(j-i)%n]*val_trz[j] for j in range(n)))

    
    tuple(c)
    return (c, normal_basis, B, inv)



def isom_normal(v, w, k1, k2, k, normal_basis_w = None, normal_basis_v = None):
    '''
    INPUT : v a normal element of k1, w a normal element of k2, k1 and k2 two 
            extensions of degree n over a base field k, k a finite field of 
            cardinality q
    OUTPUT : The image of the generator of k1 by the isomorphsim phi defined by
            phi(v) = w
            
    Algorithm :
    
    Function that explicitly computes the isomorphism phi such as phi(v) = w.

    So we need that v and w are normal elements of k1 and k2 respectively and that
    there exists an isomorphism such as phi(v) = w.

    Concretely, we have the following :

    x = sum_i c_i*v^{p^i} => phi(x) = sum_i c_i*w^{p^i}

    we are computing phi(x).
    '''
    q = k.cardinality()
    n = k1.degree()

    # We need to compute the normal basis generated by w.
    if normal_basis_w is None:
        normal_basis_w = [w]
        for i in range(n):
            normal_basis_w.append(normal_basis_w[-1]**q)

    # We begin by computing the coefficients of x in the normal basis
    # generated by v.
    normal_coefficients = convert(F.gen(), v, k, normal_basis_v)[0] 

    # Then we simply return its image in respect of the normal basis 
    # generated by w.
    return sum(normal_coefficients[i]*normal_basis_w[i]
                for i in range(n))


def calcul_isom_normal(elem, k2, img_x):
    '''
    Function that given an element of k1 computes its image in k2 by the 
    isomorphism defined by the image of k1.gen() named img_x.
    '''

    n = k2.degree()
    elem_vector = elem.vector()
    puis_img = [1]

    for i in range(n):
            puis_img.append(puis_img[-1]*img_x)

    return sum([elem_vector[i]*puis_img[i] for i in range(n)])

