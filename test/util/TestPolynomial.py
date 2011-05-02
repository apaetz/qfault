'''
Created on Nov 3, 2010

@author: adam
'''
import unittest
from util.polynomial import SymPolyWrapper, sympoly1d


#def assertEqualTypes(list1, list2):
#	assert len(list1) == len(list2)
#	typesEqual = [type(list1[i]) == type(list2[i]) for i in xrange(len(list1))]
#	assert all(typesEqual)
#
#class TestPoly1d(unittest.TestCase):
#
#	def testInit(self):
#		coeffs = [int(1), long(2), Fraction(3,4), 0]
#		p = Poly1d(coeffs)
#		for i, c in enumerate(p.coeffs):
#			assert coeffs[i] == c
#			assert type(coeffs[i]) == type(c) 
#			
#		zeros = [int(0), long(0), 0., Fraction(0,1)]
#		for z in zeros:
#			p = Poly1d(z)
#			assert 0 == p.degree
#			assert z == p[0]
#			assert type(z) == type(p[0])
#
#	def testAdd(self):
#		c11 = 2<<500
#		c10 = 5<<1000
#		c21 = 3<<500
#		c20 = 6<<1000
#		p1 = Poly1d([c11, c10])
#		p2 = Poly1d([c21, c20])
#		
#		p = p1 + p2
#		
#		assert p[0] == c10 + c20
#		assert p[1] == c11 + c21
#		assertEqualTypes(p, p1)
#		
#		typeC0 = type(p1.coeffs[-1])
#		assert typeC0 == type((p1 + 0).coeffs[-1])
#		
#	def testMul(self):
#		c11 = 2<<500
#		c10 = 5<<1000
#		c21 = 3<<500
#		c20 = 6<<1000
#		p1 = Poly1d([c11, c10])
#		p2 = Poly1d([c21, c20])
#		
#		p = p1 * p2
#		
#		assert p[2] == c11 * c21
#		assert p[1] == (c11 * c20) + (c10 * c21)
#		assert p[0] == c10 * c20
#		
#		typeC0 = type(p1.coeffs[-1])
#		p0 = (p1 * 0)
#		assert typeC0 == type((p1 * 0).coeffs[-1])
#		
##	def testGcd(self):
##		p = Poly1d([1,1])
##		g = p.gcd(p)
##		print g
##		assert g == p
#		
#		
#class TestUpperBoundedPoly(unittest.TestCase):
#	
#	def testBound(self):
#		xMin = 1.
#		xMax = 100.
#		maxDegree = 2
#		coeffs = [500, 400, 300, 200, 100, 1]
#		coeffs2 =  [250, 200, 300, 200, 100, 1]
#		p1 = Poly1d(coeffs)
#		p2 = Poly1d(coeffs2)
#		p = p1/p2
#		p = UpperBoundedPoly(p, maxDegree, xMin, xMax)
#		print p
#		print p1(xMin) / p2(xMin), p1(xMax)/p2(xMax)
#		print p(xMin), p(xMax), abs(p(xMin)-p(xMax))
#	
#	def testAdd(self):
#		xMin = 0.001/15
#		xMax = 0.002/15
#		exp = 10
#		n = 20
#		maxDegree = 2
#		p1 = Poly1d([-1,1])
#		p2 = UpperBoundedPoly(p1, maxDegree, xMin, xMax) ** exp
#		p2L = LowerBoundedPoly(p1, maxDegree, xMin, xMax) ** exp
#		p3 = p2L/p2
#		p4 = sum(p3 for _ in range(n))
#		
##		print n * (p1(xMax) ** exp) / (p1(xMax) ** exp)
##		print n* p2(xMax) /p2(xMax)
##		print n * p3(xMax)
##		print p4(xMax)
##		
##		print n * (p1(xMin) ** exp) / (p1(xMin) ** exp)
##		print n * p2(xMin) / p2(xMin)
##		print n * p3(xMin)
##		print p4(xMin)
##		
##		print 'n*p3=', n * p3
##		print 'p4=', p4


class TestSymPolyWrapper(unittest.TestCase):

	def testSum(self):
		a = SymPolyWrapper(sympoly1d([12,0]))
		b = 1-a
		c = b**200 * (a**10)
		d = sum(c for _ in range(12))
		print c
		print d
		
if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testName']
	unittest.main()