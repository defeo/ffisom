#from sage.env import sage_include_directories
#include_dirs = sage_include_directories(use_sources=True)
#import pyximport
#pyximport.install(setup_args={'include_dirs': include_dirs})

from finite_field_flint_fq_nmod import FiniteField_flint_fq_nmod as GF_flint

from rains import find_root_order
from rains import find_gens as find_gens_cyclorains
from ellrains import find_gens as find_gens_ellrains
rains_functions = [lambda k, l: find_gens_cyclorains(k, l, use_lucas=False), lambda k, l: find_gens_cyclorains(k, l, use_lucas=True), lambda k, l: find_gens_ellrains(k, l)]

from pari_ffisom import find_gens as find_gens_pari

from kummer_nmod import algolist as kummer_algolist
from kummer_nmod import namelist as kummer_namelist
from kummer_nmod import find_gens as find_gens_kummer
