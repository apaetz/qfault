'''
Use this file to generate C-code bindings for countExRecForConfigs_c.
'''

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = 'Golay Counting',
  cmdclass = {'build_ext': build_ext},
  ext_modules = [Extension("exrecCython", ["exrecCython.pyx"], 
			 extra_compile_args=['-std=c99'],
	                 libraries=['gmp']
		)]

)

