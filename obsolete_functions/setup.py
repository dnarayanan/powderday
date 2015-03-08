from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

import numpy

ext = Extension("particle_smooth_cython",
                ["particle_smooth_cython.pyx"],
                include_dirs = [numpy.get_include()]
      #          extra_compile_args=['-fopenmp'],
      #          extra_link_args=['-fopenmp']
                )

setup(ext_modules=[ext],
      cmdclass = {'build_ext' : build_ext})
