'''
Created on May 3, 2011

@author: Adam
'''

from counting.block import Block
from counting.component.transversal import TransCnot, TransMeas, TransRest
from counting.key import SyndromeKeyGenerator, MultiBlockSyndromeKeyGenerator
from qec.error import Pauli, PauliError, xType, zType
import logging
import testComponent
import unittest
from qec.ed422 import ED412Code

	
	
class TestCnot(testComponent.ComponentTestCase):
	
	def testCount(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = self.countingNoiseModels
		code = self.trivialCode
		
		cnot = self._getComponent(kGood, code)
		
		eKey = {Pauli.X: 1, Pauli.Z: 2}
		for pauli in (Pauli.X, Pauli.Z):
			result = cnot.count(noise, pauli)
		
			e = eKey[pauli]
			expected = [{(0, 0): 1}, {(0, e): 1, (e, 0): 1, (e, e): 1}]
			print result.counts
			assert result.counts == expected
		
	def testKeyPropagator(self):
		cnot = self._getComponent({}, self.trivialCode)
		blocks = [Block('0', self.trivialCode), Block('1', self.trivialCode)]
		generator = MultiBlockSyndromeKeyGenerator(blocks)
		propagator = cnot.keyPropagator()
		
		for e0 in (Pauli.I, Pauli.X, Pauli.Z, Pauli.Y):
			for e1 in (Pauli.I, Pauli.X, Pauli.Z, Pauli.Y):
				e = {'0': e0, '1': e1}
				key = generator.getKey(e)
				propagated = propagator(key)
				expected = {'0': e0 * PauliError(zbits=e1[zType]), '1': e1 * PauliError(xbits=e0[xType])}
				expectedKey = generator.getKey(expected)
#				print 'e=', e, 'key=', key, 'propagated=', propagated, 'expected=', expected
				assert propagated == expectedKey
				

	def _getComponent(self, kGood, code):
		return TransCnot(kGood, code, code)
	
	
class TestMeas(testComponent.ComponentTestCase):
	
	@staticmethod
	def Filter(e, basis):
		if Pauli.X == basis:
			e = PauliError(zbits=e[zType])
		elif Pauli.Z == basis:
			e = PauliError(xbits=e[xType])
			
		return e
	
	def testCount(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1, Pauli.Y: 1}
		noise = self.countingNoiseModels
		code = self.trivialCode
		generator = SyndromeKeyGenerator(code, None)
		
		for basis in (Pauli.Z, Pauli.X):
			meas = self._getComponent(kGood, code, basis)
			
			if basis == Pauli.X:
				expected = {Pauli.X: [{(0,): 1}, {}],
					    Pauli.Z: [{(0,): 1}, {(generator.getKey(Pauli.Z),): 1}],
					    Pauli.Y: [{(0,): 1}, {(generator.getKey(Pauli.Z),): 1}],}
			else:
				expected = {Pauli.X: [{(0,): 1}, {(generator.getKey(Pauli.X),): 1}],
					    Pauli.Z: [{(0,): 1}, {}],
					    Pauli.Y: [{(0,): 1}, {(generator.getKey(Pauli.X),): 1}],}

			
			for pauli in (Pauli.X, Pauli.Z, Pauli.Y):
				result = meas.count(noise, pauli)
#				print result.counts, expected[pauli]
				assert result.counts == expected[pauli]
		
	def testKeyPropagator(self):
		
		for basis in (Pauli.Z, Pauli.X):
			meas = self._getComponent({}, self.trivialCode, basis)
			blocks = [Block('0', self.trivialCode)]
			generator = MultiBlockSyndromeKeyGenerator(blocks)
			propagator = meas.keyPropagator()
			
			for e0 in (Pauli.I, Pauli.X, Pauli.Z, Pauli.Y):
				e = {'0': e0}
				key = generator.getKey(e)
				propagated = propagator(key)
				# Filter the part of the error that can't be detected by measurement.
				expected = {'0': self.Filter(e0, basis)}
				expectedKey = generator.getKey(expected)
#				print 'e=', e, 'key=', key, 'propagated=', propagated, 'expected=', expected
				assert propagated == expectedKey

	def _getComponent(self, kGood, code, basis=Pauli.Z):
		return TransMeas(kGood, code, basis)
		
		
class TestRest(testComponent.ComponentTestCase):
	
	def testED412(self):
		kGood = {Pauli.Y: 1}
		noiseModels = self.countingNoiseModels
		code = ED412Code()
		rest = self._getComponent(kGood, code)
		
		result = rest.count(noiseModels, Pauli.Y)
		print 'Rest [[4,2,1]]'
		print result.counts
	
	def _getComponent(self, kGood, code):
		return TransRest(kGood, code)
		
if __name__ == "__main__":

	import util.cache
	util.cache.enableFetch(False)
	util.cache.enableMemo(False)

	complog = logging.getLogger('counting.component')
	complog.setLevel(logging.DEBUG)
	
	unittest.main()