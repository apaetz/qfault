'''
Created on May 3, 2011

@author: Adam
'''

import testComponent, testBell
from counting.component.base import InputAdapter
from counting.component.teleport import TeleportED, UncorrectedTeleport
from qec.error import Pauli
import logging
import unittest
from unittest.case import SkipTest


class TestUncorrectedTeleport(testComponent.ComponentTestCase):
	
	def _getComponent(self, kGood, code):
		bellPair = testBell.TestBellPair.BellPair(kGood, code)
		bellMeas = testBell.TestBellMeas.BellMeas(kGood, code)
		return UncorrectedTeleport(kGood, bellPair, bellMeas)
	
	def testCount(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		teleport = self._getComponent(kGood, self.trivialCode)
		
		expected = {Pauli.X: [{(0, 0, 0): 1}, {(0, 1, 1): 1, (0, 1, 0): 4, (0, 0, 0): 1, (0, 0, 1): 3}],
				    Pauli.Z: [{(0, 0, 0): 1}, {(2, 0, 2): 1, (2, 0, 0): 5, (0, 0, 0): 1, (0, 0, 2): 2}]}
		
		for pauli in (Pauli.X, Pauli.Z):
			result = teleport.count(self.countingNoiseModels, pauli)
			print result.counts
			assert expected[pauli] == result.counts
		
	
class TestTeleportED(testComponent.ComponentTestCase):
	
	@staticmethod
	def TeleportED(kGood, code):
		bellPair = testBell.TestBellPair.BellPair(kGood, code)
		bellMeas = testBell.TestBellMeas.BellMeas(kGood, code)
		return TeleportED(kGood, bellPair, bellMeas)
	
	def testCount(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		teleportED = self._getComponent(kGood, self.trivialCode)
		
		eKey = {Pauli.X: (1,), Pauli.Z: (2,)}
		for pauli in (Pauli.X, Pauli.Z):
			result = teleportED.count(self.countingNoiseModels, pauli)
			expected = [{(0,): 1}, {(0,): 2, eKey[pauli]: 7}]
			print result.counts
			assert expected == result.counts
	
	@SkipTest
	def testPrAccept(self):
		kGood = {Pauli.Y: 6}
		teleportED = self._getComponent(kGood, self.trivialCode)
		prAccept = teleportED.prAccept(self.depolarizingNoiseModels)
		
		# With a trivial code, i.e., no code, the acceptance
		# probability should approach 1 as k increases.
		print prAccept
		print prAccept(0), prAccept(0.001)
		
	@SkipTest
	def testPrBad(self):
		kGood = {Pauli.X: 0}
		teleportED = self.TeleportED(kGood, self.trivialCode)
		prBad = teleportED.prBad(self.depolarizingNoiseModels[Pauli.X], Pauli.X, kMax=1)
		print prBad
		print prBad(0), prBad(0.001)
				
	def _getComponent(self, kGood, code):
		return InputAdapter(self.TeleportED(kGood, code), (0,))		
	
if __name__ == "__main__":

	import util.cache
	util.cache.enableFetch(False)
	util.cache.enableMemo(False)

	complog = logging.getLogger('counting.component')
	complog.setLevel(logging.DEBUG)
	
	unittest.main()