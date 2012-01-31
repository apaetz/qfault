'''
Created on May 3, 2011

@author: Adam
'''

from counting.component.teleport import TeleportED,\
	TeleportEDFilter, Teleport
from qec.error import Pauli, PauliError, xType, zType
from unittest.case import SkipTest
import logging
import testComponent
import testBell
import unittest
from counting.key import MultiBlockSyndromeKeyGenerator
from counting.block import Block
from qec import ed422
from counting.result import CountResult, TrivialResult
import counting
import qec


class TestTeleport(testComponent.ComponentTestCase):
	
	def _getComponent(self, kGood, code):
		bellPair = testBell.TestBellPair.BellPair(kGood, code)
		bellMeas = testBell.TestBellMeas.BellMeas(kGood, code)
		return Teleport(kGood, bellPair, bellMeas)
	
	def testCount(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		teleport = self._getComponent(kGood, self.trivialCode)
		
		expected = {Pauli.X: [{(0, 0, 0): 1}, {(0, 1, 0): 1, (0, 1, 1): 4, (0, 0, 0): 1, (0, 0, 1): 3}],
				    Pauli.Z: [{(0, 0, 0): 1}, {(2, 0, 0): 1, (2, 0, 2): 5, (0, 0, 0): 1, (0, 0, 2): 2}]}
		
		for pauli in (Pauli.X, Pauli.Z):
			inputResult = testComponent.trivialInput(teleport)
			result = teleport.count(self.countingNoiseModels, pauli, inputResult)
			print result.counts
			assert expected[pauli] == result.counts
		
	
class TestTeleportED(testComponent.ComponentTestCase):
	
	@staticmethod
	def TeleportED(kGood, code):
		bellPair = testBell.TestBellPair.BellPair(kGood, code)
		bellMeas = testBell.TestBellMeas.BellMeas(kGood, code)
		return TeleportED(kGood, bellPair, bellMeas)
	
	def testCount(self):
		kGood = {Pauli.X: 2, Pauli.Z: 2}
		teleportED = self._getComponent(kGood, self.trivialCode)
		
		eKey = {Pauli.X: (1,), Pauli.Z: (2,)}
		for pauli in (Pauli.X, Pauli.Z):
			result = teleportED.count(self.countingNoiseModels, pauli)
			expected = [{(0,): 1}, {(0,): 2, eKey[pauli]: 7}, {(0,): 20, eKey[pauli]: 10}]
			print result.counts
			assert expected == result.counts
		
	def testPrAccept(self):
		kGood = {Pauli.Y: 6}
		teleportED = self._getComponent(kGood, self.trivialCode)
		prAccept = teleportED.prAccept(self.depolarizingNoiseModels, TrivialResult(teleportED.inBlocks()), None)
		
		# With a trivial code, i.e., no code, the acceptance
		# probability should approach 1 as k increases.
		print prAccept
		print prAccept(0), prAccept(0.001)
		
	def testPrBad(self):
		kGood = {Pauli.X: 0}
		teleportED = self.TeleportED(kGood, self.trivialCode)
		prBad = teleportED.prBad(self.depolarizingNoiseModels[Pauli.X], Pauli.X, kMax=1)
		print prBad
		print prBad(0), prBad(0.001)
				
	def _getComponent(self, kGood, code):
		return self.TeleportED(kGood, code)		
#	
class TestTeleportEDFilter(testComponent.ComponentTestCase):
	
	def testCount(self):
		kGood = {Pauli.X: 2, Pauli.Z: 2}
		
		code = ed422.ED412Code()
		Xl = code.logicalOperators()[0][xType]
		Zl = code.logicalOperators()[0][zType]
		
		teleportED = self._getComponent(kGood, code)
		
		# These errors should be undetected, and pass through
		III = {0: Pauli.I, 1: Pauli.I, 2: Pauli.I}
		IIX = {0: Pauli.I, 1: Pauli.I, 2: Pauli.X}
		IXlX = {0: Pauli.I, 1: Xl, 2: Pauli.X}
		ZlIX = {0: Zl, 1: Pauli.I, 2: Pauli.X}
		
		# These errors should be detected, and not pass through.
		IXI = {0: Pauli.I, 1: Pauli.X, 2: Pauli.I}
		ZII = {0: Pauli.Z, 1: Pauli.I, 2: Pauli.I}
		
		blocks = tuple(Block(i, code, None) for i in range(3))
		
		generator = MultiBlockSyndromeKeyGenerator(blocks)
		counts = {}
		expected = {}
		notDetected = set([Pauli.I, Xl, Zl])
		for i, err in enumerate([III, IIX, IXlX, ZlIX, IXI, ZII]):
			key = generator.getKey(err)
			counts[key] = i+1
			
			if err[0] in notDetected and err[1] in notDetected:
				k = (key[2],)
				expected[k] = expected.get(k, 0) + i+1
		
		expected = [expected]
		
		inputResult = CountResult([counts], teleportED.inBlocks())
		result = teleportED.count(self.countingNoiseModels, Pauli.Y, inputResult=inputResult)
#		print result.counts, expected
		assert expected == result.counts
	
	def _getComponent(self, kGood, code):
		teleport = TeleportEDFilter(code)
		return teleport
	
if __name__ == "__main__":

	import util.cache
	util.cache.enableFetch(False)
	util.cache.enableMemo(False)

	complog = logging.getLogger('counting.component')
	complog.setLevel(logging.DEBUG)
	
	unittest.main()