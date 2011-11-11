'''
Created on May 3, 2011

@author: Adam
'''

import logging
from counting.component.ec import TECDecodeAdapter, ConcatenatedTEC
from unittest.case import SkipTest
from counting import probability
import testBase, testComponent
logging.basicConfig(level=logging.DEBUG)

from counting.component.transversal import TransCnot
from counting.component.base import Prep, Empty, InputAdapter,\
	ConcatenatedComponent
from counting.component.bell import BellPair, BellMeas
from counting.component.exrec import ExRec
from counting.component.rectangle import CnotRectangle
from counting.component.teleport import UncorrectedTeleport, TeleportED
from counting.countErrors import mapCounts
from counting.location import Locations
from qec.error import Pauli, xType, zType
from qec.qecc import TrivialStablizerCode, StabilizerState
from settings.noise import NoiseModelXSympy, CountingNoiseModel, \
	CountingNoiseModelX, CountingNoiseModelZ, NoiseModelZSympy, NoiseModelXZLowerSympy
from util import counterUtils
import golay
import unittest

	
	


		
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
		expected = [{(0, 0): 1}, {(0, 1): 2, (1, 0): 1, (1, 1): 1}]
		assert result.counts == expected
		
	def _getComponent(self, kGood, code):
		return self.BellPair(kGood, code)

class TestBellMeas(testComponent.ComponentTestCase):
	
	@staticmethod
	def BellMeas(kGood, code):
		return BellMeas(kGood, code)
	
	def testX(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = {Pauli.X: CountingNoiseModelX()}
		code = self.trivialCode
		
		bm = self.BellMeas(kGood, code)
		result = bm.count(noise, Pauli.X)
		expected = [{(0, 0): 1}, {(0, 0): 1, (0, 1): 3}]
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