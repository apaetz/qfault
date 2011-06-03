'''
Created on 2010-10-19

This file contains classes for creating and manipulating polynomials.

@author: adam
'''
import logging
logger = logging.getLogger('polynomial')
import operator

def disableSympyCache():
	'''
	Disabled the sympy global cache.  This function must be called before
	importing sympy. 
	By default sympy uses a global cache.  This cache seems to be buggy
	when using multiple threads.  The workaround is to disable the cache
	entirely.
	'''
		
	import os
	os.putenv('SYMPY_USE_CACHE', 'no')
	import sympy
	if sympy.cache.usecache != 'no':
		raise Exception('Unable to turn off sympy cache.  Sympy may have been imported already.')
	logger.info('Disabled sympy cache.')

disableSympyCache()

import sympy
from sympy import Symbol
from sympy.simplify.simplify import powsimp, collect
from sympy.functions.combinatorial.factorials import Factorial




def chebyshevT(n, symbol='x'):
	'''
	Returns the degree-n Chebyshev polynomial of the first kind T_n(x).
	Note: This implemenation is optimized for taking derivatives, and only 
	works for n > 0.
	
	>>> T3 = chebyshevT(3)
	>>> T3(2)
	26
	'''
	x = Symbol(symbol)
	term = lambda k: (-2)**k * Factorial(n+k-1) / (Factorial(n-k) * Factorial(2*k)) * (1-x)**(k)
	Tn = n * sum(term(k) for k in range(n + 1))
	return SymPolyWrapper(Tn) 

def sympoly1d(coeffs, symbol='x'):
	'''
	Returns a Sympy expression for the univariate one-dimensional
	polynomal represented by the given coefficients.
	
	coeffs - A list of coefficients. coeffs[0] is the highest order coefficient,
	         coeffs[-1] is the zero-order coefficient.
	'''
	x = Symbol(symbol)
	poly = 0
	for k, c in enumerate(reversed(coeffs)):
		poly += c * (x ** k)
	
	return poly
		

class SymPolyWrapper(object):
	'''
	Provides is nicer interface to sympy expressions.
	In particular, allows for substitution without having specify the
	symbol.  For example poly(1.5) rather than poly.subs(x,1.5).
		
	Arguments
	----------
	coeffs : array_like
		polynomial coefficients, in decreasing powers.  For example,
		``[1, 2, 3)]`` implies :math:`x^2 + 2x + 3`.

	'''
	
	__slots__ = ('_poly')
	
	def __init__(self, sympoly):
		if isinstance(sympoly, SymPolyWrapper):
			sympoly = sympoly._poly
		
		self._poly = sympoly  
		
	def simplify(self):
		syms = list(self._poly.atoms(Symbol))
		simp = collect(self._poly, syms)
		simp = powsimp(simp)
		return self.__class__(simp)

	def _operateBinary(self, other, op):
		other = self.__class__(other)
		return self.__class__(op(self._poly, other._poly))
	
	def __add__(self, other):
		return self._operateBinary(other, operator.add)
	
	def __radd__(self, other):
		return self + other
	
	def __sub__(self, other):
		return self._operateBinary(other, operator.sub)
	
	def __rsub__(self, other):
		return -(self - other)
	
	def __mul__(self, other):
		return self._operateBinary(other, operator.mul)
	
	def __rmul__(self, other):
		return self * other
	
	def __pow__(self, exp):
		return self._operateBinary(exp, operator.pow)
	
	def __div__(self, other):
		return self.__truediv__(other)
	
	def __rdiv__(self, other):
		other = self.__class__(other)
		return other / self

	def __truediv__(self, other):
		return self._operateBinary(other, operator.truediv)
	
	def __divmod__(self, other):
		return NotImplemented
	
	def __neg__(self):
		return self*(-1)
	
	def __call__(self, val, symbols=None):
		if None == symbols:
			symbols = self._poly.atoms(Symbol)
			
		# Substitute val, for all symbols in the
		# expression.
		result = self._poly
		for x in symbols:
			result = result.subs(x, val)
			
		return result
	
	def diff(self, symbol=Symbol('x'), order=1):
		return self.__class__(sympy.diff(self._poly, symbol, order))
	
	def asNumDen(self):
		n, d = self._poly.as_numer_denom()
		return self.__class__(n), self.__class__(d)
	
	def getsympy(self):
		return self._poly
	
	def __eq__(self, other):
		return self._poly == other
		
	def __repr__(self):
		return 'SymPolyWrapper(%s)' % self._poly
	
	def __str__(self):
		return str(self._poly)

	
class Composite(object):
	'''
	Represents a composite function f(g(x)).
	'''
	
	def __init__(self, f, g):
		self._f = f
		self._g = g
	
	def __call__(self, *args, **kwargs):
		return self._f(self._g(*args, **kwargs))

if __name__ == "__main__":
	import doctest
	doctest.testmod()