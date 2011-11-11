'''
Created on May 3, 2011

@author: Adam
'''

from counting.component.transversal import TransCnot
from qec.error import Pauli
from qec.qecc import TrivialStablizerCode
from settings.noise import CountingNoiseModelX
import logging
import unittest
import testComponent

	
	
class TestCnot(testComponent.ComponentTestCase):
	
	def testX(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = {Pauli.X: CountingNoiseModelX()}
		code = self.trivialCode
		
		cnot = self._getComponent(kGood, code)
		result = cnot.count(noise, Pauli.X)
		
		expected = [{(0, 0): 1}, {(0, 1): 1, (1, 0): 1, (1, 1): 1}]
		assert result.counts == expected
		
	def _getComponent(self, kGood, code):
		return TransCnot(kGood, code, code)
		
		
if __name__ == "__main__":

	import util.cache
	util.cache.enableFetch(False)
	util.cache.enableMemo(False)

	complog = logging.getLogger('counting.component')
	complog.setLevel(logging.DEBUG)
	
	unittest.main()