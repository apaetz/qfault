'''
Created on May 3, 2011

@author: Adam
'''

from copy import copy
from counting.component.adapter import InputAdapter
from counting.component.base import Empty
from counting.component.ec import DecodeAdapter, LECSyndromeAdapter
from counting.component.exrec import ExRec
from qec.error import Pauli
import logging
import testComponent
import testTeleport
import unittest
from counting.key import SyndromeKeyGenerator

	
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
		key = {pauli: (SyndromeKeyGenerator(self.trivialCode, None).getKey(pauli),) for pauli in kGood.keys()}
		expected = {pauli: [{(0,): 1}, {(0,): 4, key[pauli]: 2*7}, {(0,): 7*7 + 2*2 + 2*20, key[pauli]: 2*7*2 + 2*10}] for pauli in (Pauli.X, Pauli.Z)}
		checkCount(self.trivialCode, expected)
			

#	def testXCnot(self):
#		kGood = {Pauli.X: 1, Pauli.Z: 1}
#		noise = {Pauli.X: CountingNoiseModelX()}
#		code = self.trivialCode
#		cnot = TransCnot(kGood, code, code) 
#		
#		ec = TestTeleportED.TeleportED(kGood, code)
#		lec = InputAdapter(ec, (0,))
#		lec = ConcatenatedComponent(kGood, lec, lec)
#		tec = DecodeAdapter(ec)
#		tec = ConcatenatedTEC(kGood, tec, tec)
#					
#		exRec = ExRec(kGood, lec, cnot, tec)
#		result = exRec.count(noise, Pauli.X)
#		expected = [{Pauli.I: 1}, {Pauli.X: 2, Pauli.I: 4}]
#		#print 'cnot exRec:' + str(result.counts)# == expected
		
#
#		
#		cnot = TransCnot(kGood, code, code) 
#		code = self.trivialCode
#		ec = TestTeleportED.TeleportED(kGood, code)
#		exRec = ExRec(kGood, lec, cnot, tec)
#		kGood = {Pauli.X: 1, Pauli.Z: 1}
#		lec = ConcatenatedComponent(kGood, lec, lec)
#		lec = InputAdapter(ec, (0,))
#		noise = self.countingNoiseModels
#		pr = exRec.prAccept(noise)
#		print 'pr(0.001)', pr(0.001)		
#		print pr
#		print pr(0.000)
#		tec = ConcatenatedTEC(kGood, tec, tec)
#		tec = DecodeAdapter(ec)
#	@SkipTest
#	def testPrAccept(self):

	def testProbabilities(self):
		# The normal probability test doesn't work for the Forward ExRec because
		# some of the counts are filtered out after the Ga.
		pass
	
		
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