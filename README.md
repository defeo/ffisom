# Computing isomorphisms and embeddings of finite fields

*a research project in Computer Algebra*

by [Ludovic Brieulle](https://github.com/brieulle/), [Luca De
Feo](http://defeo.lu/), [Javad
Doliskani](https://github.com/javad-doliskani), [Jean-Pierre
Flori](https://github.com/jpf-anssi/) and [Éric
Schost](http://www.csd.uwo.ca/~eschost/).

#### Abstract

Let *𝔽<sub>q</sub>* be a finite field.  Given two irreducible
polynomials *f*, *g* over *𝔽<sub>q</sub>*, with deg *f* dividing deg
*g*, the finite field embedding problem asks to compute an explicit
description of a field embedding of *𝔽<sub>q</sub>[X]/f(X)* into
*𝔽<sub>q</sub>[Y]/g(Y)*.  When deg *f* = deg *g*, this is also known
as the isomorphism problem.

This problem, a special instance of polynomial factorization, plays a
central role in computer algebra software.  We review previous
algorithms, due to Lenstra, Allombert, Rains, and Narayanan, and
propose improvements and generalizations.  Our detailed complexity
analysis shows that our newly proposed variants are at least as
efficient as previously known algorithms, and in many cases
significantly better.

We also implement most of the presented algorithms, compare them with
the state of the art computer algebra software, and make the code
available as open source.  Our experiments show that our new variants
consistently outperform available software.


## Research paper

The research paper has been submitted for publication in [AMS
Mathematics of
Computation](http://www.ams.org/publications/journals/journalsframework/mcom).

A preprint version is available at <https://arxiv.org/abs/1705.01221>.

The LaTeX sources to the paper can be found in the [`paper`](paper)
folder.  They are subject to the copying restriction enforced by the
editor.


## Software [![launch binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/defeo/ffisom/master?filepath=notebooks%2Fexample.ipynb)

The [`implementation`](implementation) folder contains source code for
the algorithms implemented by the project. See the folder `README.md`
for instructions on compiling and running the software on your own
machine.

The [`notebooks`](notebooks) folder contains Jupyter notebooks showing
sample usage of the software and benchmarks. Notebooks can be
statically viewed on GitHub, or **executed** in
[Binder](https://mybinder.org/v2/gh/defeo/ffisom/master?filepath=notebooks%2Fexample.ipynb),
thus saving you the time to compile and install the software.

All source code is distributed under the MIT license.


## Related material

The [`misc`](misc) folder contains material related to the project,
such as

- conference posters,
- presentation slides,
- experimental data gathered on the conjecture formulated in the
  paper.


## Contributing

We do not accept pull requests, but if you see any problem with this
repository, plese open an issue.
