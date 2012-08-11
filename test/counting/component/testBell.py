'''
Created on May 3, 2011

@author: Adam
'''

from counting.component.bell import BellPair, BellMeas
from qec.error import Pauli
from settings.noise import CountingNoiseModelX
import logging
import testBase
import testComponent
import unittest
logging.basicConfig(level=logging.DEBUG)

		
class TestBellPair(testComponent.ComponentTestCase):
	
	@staticmethod
	def BellPair(kGood, code):
		plus = testBase.TestPrep.Plus(kGood, code)
		zero = testBase.TestPrep.Zero(kGood, code)
		return BellPair(kGood, plus, zero, kGood)				
	
	def testX(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = {Pauli.X: CountingNoiseModelX()}
		code = self.trivialCode
		
		bp = self.BellPair(kGood, code)
		result = bp.count(noise, Pauli.X)
		expected = [{(0, 0): 1}, {(0, 1): 2, (1, 0): 1, (0, 0): 1}]
		assert result.counts == expected
		
	def _getComponent(self, kGood, code):
		return self.BellPair(kGood, code)

class TestBellMeas(testComponent.ComponentTestCase):
	
	@staticmethod
	def BellMeas(kGood, code):
		return BellMeas(kGood, code)
	
	def testCount(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = self.countingNoiseModels
		code = self.trivialCode
		
		bm = self.BellMeas(kGood, code)
		
		eKey = {Pauli.X: (0,1), Pauli.Z: (2,0)}
		for pauli in (Pauli.X, Pauli.Z):
			result = bm.count(noise, pauli)
			expected = [{(0, 0): 1}, {(0, 0): 1, eKey[pauli]: 3}]
			assert result.counts == expected
		
	def _getComponent(self, kGood, code):
		return self.BellMeas(kGood, code)
		

if __name__ == "__main__":

	import util.cache
	util.cache.enableFetch(False)
	util.cache.enableMemo(False)

	complog = logging.getLogger('counting.component')
	complog.setLevel(logging.DEBUG)
	
	unittest.main()