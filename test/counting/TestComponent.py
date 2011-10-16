'''
Created on May 3, 2011

@author: Adam
'''
import unittest
import golay
from counting.component import Prep


class TestPrepZero(unittest.TestCase):


	def testGolay(self):
		locations = golay.ancillaPrep.getSteaneRandomPreps()[0]
		kGood = 1
		kBest = 1
		code = golay.golayCode.
		prep = Prep(kGood, kBest, locations, code) 


if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testName']
	unittest.main()