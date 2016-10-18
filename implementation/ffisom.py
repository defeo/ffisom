#from sage.env import sage_include_directories
#include_dirs = sage_include_directories(use_sources=True)
#import pyximport
#pyximport.install(setup_args={'include_dirs': include_dirs})

from rains import find_gens as find_gens_rains
from ellrains import find_gens as find_gens_ellrains
from pari_ffisom import find_gens as find_gens_pari
from kummer_nmod import find_gens as find_gens_kummer
from finite_field_flint_fq_nmod import FiniteField_flint_fq_nmod as GF_flint
