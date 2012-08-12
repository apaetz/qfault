'''
Created on Nov 3, 2010

@author: adam
'''
from counting.probability import prMinFailures
from counting import count_parallel
from counting.location import LocationCount
from settings.noise import NoiseModelXSympy
import unittest


class Test(unittest.TestCase):

#
#    def testPrecision(self):
#        locTotals = LocationCount(10, 10, 10, 10, 10, 10)
#        noise = NoiseModelXSympy(gMin=0.001/15, gMax=0.002/15)
#        p = (1, locTotals, noise, None)
#        print repr(p)
#        print repr(p.simplify())
#        print p(0.001/15)
		
	def testPrMinFailures(self):
		#locTotals = LocationCount(10, 10, 10, 10, 10, 10)
		locTotals = LocationCount(25, 0, 0, 0, 0, 0)
		noise = NoiseModelXSympy(gMin=0.001/15, gMax=0.002/15) 
		pr = prMinFailures(5, locTotals, noise)
		print repr(pr)
		print pr(0.001/15)


if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testName']
	count_parallel.configureMultiProcess(0)
	unittest.main()