from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from sage.env import SAGE_ENV, sage_include_directories

extensions = [
              Extension("*", sources = ["*.pyx"],
                        include_dirs = sage_include_directories() + [SAGE_ENV['SAGE_INC'] + '/flint'] + ['ellmul'])
             ]
setup(
    ext_modules = cythonize(extensions, include_path=SAGE_ENV['SITE_PACKAGES'] + ['.'], build_dir='cythonized')
)
