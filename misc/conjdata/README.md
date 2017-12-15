The sagemath file conj.sage defines functions to check the conjecture
using the three strategies described in the long version of the paper:
* strategy one is "brute-force": for a given range of p and ell check all
  curves defined over GF(p) by finding torsion points over the extension
  of degree r and computing the period (trace-like or variation of it
  using other symmetric functions).
* strategy two is similar but instead of computing a torsion point over
  the extension of degree r, polynomials gving the multiplication map
  over the base field are used.
* strategy three is using modular curves to generate curves over number
  fields with torsion over extensions of a given degree r and trying to
  find a prime p where reduction would yield a counter example over
  the finite field of characteristic p.

The text file conj.dat gives data on the ranges where the conjecture
and its variations were checked.
To summarize:
* No counter example was found for the actual conjecture (using the trace)
  though some millions of curves where tested.
* For variations using other symmetric functions counter examples were
  quickly found.
