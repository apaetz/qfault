'''
Created on May 3, 2011

@author: Adam
'''

import testComponent, testBell
from counting.component.base import InputAdapter
from counting.component.teleport import TeleportED
from qec.error import Pauli
import logging
import unittest

	
	
class TestTeleportED(testComponent.ComponentTestCase):
	
	@staticmethod
	def TeleportED(kGood, code):
		bellPair = testBell.TestBellPair.BellPair(kGood, code)
		bellMeas = testBell.TestBellMeas.BellMeas(kGood, code)
		return TeleportED(kGood, bellPair, bellMeas)
	
	def testPrAccept(self):
		kGood = {Pauli.Y: 6}
		teleportED = self.TeleportED(kGood, self.trivialCode)
		teleportED = InputAdapter(teleportED, (0,))
		prAccept = teleportED.prAccept(self.depolarizingNoiseModels)
		
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
		return InputAdapter(self.TeleportED(kGood, code), (0,))			
	
if __name__ == "__main__":

	import util.cache
	util.cache.enableFetch(False)
	util.cache.enableMemo(False)

	complog = logging.getLogger('counting.component')
	complog.setLevel(logging.DEBUG)
	
	unittest.main()