from distutils.core import setup
from Cython.Build import cythonize

#python setupXZ.py build_ext --inplace
setup(
    ext_modules = cythonize("XZ.pyx")
)
