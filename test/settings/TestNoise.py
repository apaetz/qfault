'''
Created on 2011-04-29

@author: adam
'''
import unittest
from settings.noise import NoiseModelXSympy
from sympy.mpmath.mptypes import mpf

locCnot = {'type': 'cnot'}
locRest = {'type': 'rest'}
locPrepX = {'type': 'prepX'}
locPrepZ = {'type': 'prepZ'}
locMeasX = {'type': 'measX'}
locMeasZ = {'type': 'measZ'}

locations = [locCnot, locRest, locPrepX, locPrepZ, locMeasX, locMeasZ]

class TestNoiseModelXSympy(unittest.TestCase):

	def setUp(self):
		self.noiseModel = NoiseModelXSympy(0,1)

	def testWeights(self):
		expectedWeights = [4,8,4,4,4,4]
		
		for l in range(len(locations)):
			loc = locations[l]
			expWt = expectedWeights[l]
			for error in range(self.noiseModel.numErrors(loc)):
				assert self.noiseModel.getWeight(loc, error) == expWt
		
	def testNumErrors(self):
		expectedErrors = [3,1,1,1,1,1]
	
		for l in range(len(locations)):
			loc = locations[l]
			expNum = expectedErrors[l]
			assert self.noiseModel.numErrors(loc) == expNum
			
	def testLikelyhood(self):
		# check interval gamma \in [0,1]
		gammas = [mpf(.0001*i) for i in range(1001)]
		likelyhood = self.noiseModel.likelyhood()
		for gamma in gammas:
			# check accuracy up to ten digits.
			assert int(likelyhood(gamma) * 10**10) == int(gamma/(1-12*gamma) * 10**10)

if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testName']
	unittest.main()
	