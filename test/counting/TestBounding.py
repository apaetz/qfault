'''
Created on 2011-04-17

@author: adam
'''
from counting.bounding import computeWeights
from util.polynomial import SymPolyWrapper, sympoly1d
import logging
import unittest


class TestBounding(unittest.TestCase):

	def testComputeWeights(self):
		polys = [SymPolyWrapper(sympoly1d([i, 0, i])) for i in range(1,5)]
		refPoly, weights = computeWeights(polys, 0, 1)
		print refPoly, weights


if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testComputeWeights']
	logging.basicConfig(level=logging.DEBUG)
	unittest.main()