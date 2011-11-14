'''
Created on May 3, 2011

@author: Adam
'''

from counting.block import Block
from counting.component.transversal import TransCnot
from counting.key import SyndromeKeyGenerator, MultiBlockSyndromeKeyGenerator
from qec.error import Pauli
from qec.qecc import TrivialStablizerCode
from settings.noise import CountingNoiseModelX
import logging
import testComponent
import unittest

	
	
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
		propagator = cnot.keyPropagator(generator.keyMeta())
		
		for e0 in (Pauli.I, Pauli.X, Pauli.Z, Pauli.Y):
			for e1 in (Pauli.I, Pauli.X, Pauli.Z, Pauli.Y):
				e = {'0': e0, '1': e1}
				key = generator.getKey(e)
				propagated = propagator(key)
				print 'e=', e, 'key=', key, 'propagated=', propagated

	def _getComponent(self, kGood, code):
		return TransCnot(kGood, code, code)
		
		
if __name__ == "__main__":

	import util.cache
	util.cache.enableFetch(False)
	util.cache.enableMemo(False)

	complog = logging.getLogger('counting.component')
	complog.setLevel(logging.DEBUG)
	
	unittest.main()