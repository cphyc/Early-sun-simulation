from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "Menou paper cython implementation",
    ext_modules = cythonize('crunner.pyx'), # accepts a glob pattern
)
