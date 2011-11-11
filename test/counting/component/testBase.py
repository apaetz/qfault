'''
Created on May 3, 2011

@author: Adam
'''

from counting.component.base import Prep
import testComponent
from counting.location import Locations
from qec.error import Pauli, xType, zType
from qec.qecc import TrivialStablizerCode, StabilizerState
from settings.noise import NoiseModelXSympy, CountingNoiseModel, \
	CountingNoiseModelX, CountingNoiseModelZ
from util import counterUtils
import logging
import unittest
logging.basicConfig(level=logging.DEBUG)


class TestPrep(testComponent.ComponentTestCase):
	
	noise = {Pauli.X: CountingNoiseModelX(),
			 Pauli.Z: CountingNoiseModelZ()}
	code = TrivialStablizerCode()
	
	@staticmethod
	def Zero(kGood, code):
		return Prep(kGood, Locations([counterUtils.locZprep('|0>', 0)], '|0>'), StabilizerState(code, [zType]))
	
	@staticmethod
	def Plus(kGood, code):
		return Prep(kGood, Locations([counterUtils.locXprep('|+>', 0)], '|+>'), StabilizerState(code, [xType]))
	
	def _getComponent(self, kGood):
		return self.Zero(kGood, self.trivialCode)
	
	def testCount(self):
		expectedZero = [{(0,): 1}, {(1,): 1}]
		expectedPlus = [{(0,): 1}, {}]
		
		self.doCountTest(Pauli.X, expectedPlus, expectedZero)
		
		expectedPlus = [{(0,): 1}, {(2,): 1}]
		expectedZero = [{(0,): 1}, {}]
		self.doCountTest(Pauli.Z, expectedPlus, expectedZero)
		
	def doCountTest(self, pauli, expectedPlus, expectedZero):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		
		plus = self.Plus(kGood, self.code)
		result = plus.count(self.noise, pauli)
		assert result.counts == expectedPlus
		
		zero = self.Zero(kGood, self.code)
		result = zero.count(self.noise, pauli)
		assert result.counts == expectedZero
		
	def testPrBad(self):
		logger = logging.getLogger('counting.probability')
		logger.setLevel(logging.DEBUG)
		
		kGood = {Pauli.X: 0, Pauli.Z: 1}
		pauli = Pauli.Z
		
		plus = self.Plus(kGood, self.code)
		
		# In the Z-error case we are counting all (1) locations
		# so there are no bad cases.
		prBad = plus.prBad(self.noise[Pauli.Z], pauli, 1)
		assert 0 == prBad
		
		# In the X-error case, |+> does not cause any
		# generate X-errors, so there are again no bad cases.
		prBad = plus.prBad(self.noise[Pauli.X], pauli, 1)
		assert 0 == prBad
		
		# Now we count zero of the (1) locations.  So
		# everything is bad.
		kGood[Pauli.Z] = 0
		plus = self.Plus(kGood, self.code)
		prBad = plus.prBad(self.noise[Pauli.Z], pauli, 1)
		assert 1 == prBad
		
	def testPrAccept(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		plus = self.Plus(kGood, self.code)
		zero = self.Zero(kGood, self.code)
		
		paulis = (Pauli.X, Pauli.Z)
		assert all((1 == plus.prAccept(self.noise, 1) for pauli in paulis))
		assert all((1 == zero.prAccept(self.noise, 1) for pauli in paulis))		

if __name__ == "__main__":

	import util.cache
	util.cache.enableFetch(False)
	util.cache.enableMemo(False)

	complog = logging.getLogger('counting.component')
	complog.setLevel(logging.DEBUG)
	
	#import sys;sys.argv = ['', 'Test.testName']
	unittest.main()