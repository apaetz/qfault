'''
Created on May 3, 2011

@author: Adam
'''

from counting.component.adapter import DecodeAdapter
from counting.component.base import Empty
from counting.component.exrec import ExRec
from counting.key import SyndromeKeyGenerator
from qec.error import Pauli
import logging
import testComponent
import testTeleport
import unittest

	
class TestExRecTeleportED(testComponent.ComponentTestCase):

	def testCount(self):
		kGood = {Pauli.X: 2, Pauli.Z: 2}
		noise = self.countingNoiseModels
		
		def checkCount(code, expected):
			exRec = self._getComponent(kGood, code)
			for pauli in (Pauli.X, Pauli.Z):
				result = exRec.count(noise, pauli)
				assert result.counts == expected[pauli]
				
		# There are 7 ways to produce an X/Z-error from the output of the
		# an ED with one fault, and 2 ways to produce an I-error with one fault.  
		# There are 10 ways to get an X with two faults in an ED, and 20 ways to
		# get an I-error.
		key = {pauli: (SyndromeKeyGenerator(self.trivialCode, None).get_key(pauli),) for pauli in kGood.keys()}
		expected = {pauli: [{(0,): 1}, {(0,): 4, key[pauli]: 2*7}, {(0,): 7*7 + 2*2 + 2*20, key[pauli]: 2*7*2 + 2*10}] for pauli in (Pauli.X, Pauli.Z)}
		checkCount(self.trivialCode, expected)
	
		
	def _getComponent(self, kGood, code):
		empty = Empty(code) 
		
		ec = testTeleport.TestTeleportED.TeleportED(kGood, code)
		lec = ec
		tec = DecodeAdapter(ec)
		
		exRec = ExRec(kGood, lec, empty, tec)
		return exRec

if __name__ == "__main__":

	import util.cache
	util.cache.enableFetch(False)
	util.cache.enableMemo(False)

	complog = logging.getLogger('counting.component')
	complog.setLevel(logging.DEBUG)
	
	unittest.main()