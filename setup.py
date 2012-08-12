'''
Use this file to generate C-code bindings for braidview.pyx.
'''

#from distutils.core import setup
#from distutils.extension import Extension
from setuptools import setup, find_packages, Extension

setup(
  name = 'qfault',
  version='0.2a',
  description='Threshold and overhead calculation for fault-tolerant quantum circuits.',
  author='Adam Paetznick',
  author_email='adampaetznick@gmail.com',
  
  # Python packages
  package_dir = {'': 'qfault'},
  packages = find_packages('qfault'),
  
  install_requires = ['sympy >= 0.7.1rc1', 'gmpy >= 1.11', 'numpy >= 1.3.0'],
)

