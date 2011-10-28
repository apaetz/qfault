'''
Created on May 3, 2011

@author: Adam
'''

import logging
from counting.component.transversal import TransCnot
logging.basicConfig(level=logging.DEBUG)

from counting.component.base import Prep, Empty, InputAdapter
from counting.component.bell import BellPair, BellMeas
from counting.component.exrec import ExRec
from counting.component.rectangle import CnotRectangle
from counting.component.teleport import UncorrectedTeleport, TeleportED
from counting.countErrors import mapCounts
from counting.location import Locations
from qec.error import Pauli, xType, zType
from qec.qecc import TrivialStablizerCode, StabilizerState
from settings.noise import NoiseModelXSympy, CountingNoiseModel, \
	CountingNoiseModelX, CountingNoiseModelZ
from util import counterUtils
import golay
import unittest






#class TestTeleport(unittest.TestCase):
#
#
#	def testX(self):
#		kGood = {Pauli.X: 1, Pauli.Z: 1}
#		noise = {Pauli.X: CountingNoiseModelX()}
#		
#		code = TrivialStablizerCode()
#		plus = Prep(kGood, Locations([counterUtils.locXprep('|+>', 0)], '|+>'), StabilizerState(code, [xType]))
#		zero = Prep(kGood, Locations([counterUtils.locZprep('|0>', 0)], '|0>'), StabilizerState(code, [zType]))
#		bellPair = BellPair(kGood, plus, zero, kGood)
#		bellMeas = BellMeas(kGood, code)
#		data = Empty(code, 'data').count(noise, Pauli.X)
#		
#		teleport = UncorrectedTeleport(kGood, bellPair, bellMeas)
#		result = teleport.count(noise, Pauli.X)
#		
#		expectedCounts = [{(0, 0, 0): 1}, {(1, 0, 0): 1, (1, 1, 0): 1, (0, 1, 0): 3, (0, 0, 1): 2, (0, 1, 1): 1}]
#		assert result.counts == expectedCounts
#		
#		#print result.counts
#		
#	def testZ(self):
#		kGood = {Pauli.Z: 1}
#		noise = {Pauli.Z: CountingNoiseModelZ()}
#		
#		code = TrivialStablizerCode()
#		plus = Prep(kGood, Locations([counterUtils.locXprep('|+>', 0)], '|+>'), StabilizerState(code, [xType]))
#		zero = Prep(kGood, Locations([counterUtils.locZprep('|0>', 0)], '|0>'), StabilizerState(code, [zType]))
#		bellPair = BellPair(kGood, plus, zero, kGood)
#		bellMeas = BellMeas(kGood, code)
#		data = Empty(code, 'data').count(noise, Pauli.Z)
#		
#		teleport = UncorrectedTeleport(kGood, bellPair, bellMeas)
#		result = teleport.count(noise, Pauli.Z)
#		
#		expectedCounts = [{(0, 0, 0): 1}, {(2, 0, 0): 2, (0, 2, 0): 1, (2, 2, 0): 3, (0, 0, 2): 1, (2, 2, 2): 1}]
#		assert result.counts == expectedCounts
#		
#		#print result.counts
		
		
class TestPrep(unittest.TestCase):
	
	@staticmethod
	def Zero(kGood, code):
		return Prep(kGood, Locations([counterUtils.locZprep('|0>', 0)], '|0>'), StabilizerState(code, [zType]))
	
	@staticmethod
	def Plus(kGood, code):
		return Prep(kGood, Locations([counterUtils.locXprep('|+>', 0)], '|+>'), StabilizerState(code, [xType]))
	
	def testX(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = {Pauli.X: CountingNoiseModelX()}
		code = TrivialStablizerCode()
		
		plus = self.Plus(kGood, code)
		result = plus.count(noise, Pauli.X)
		expected = [{(0,): 1}, {}]
		assert result.counts == expected
		
		zero = self.Zero(kGood, code)
		result = zero.count(noise, Pauli.X)
		expected = [{(0,): 1}, {(1,): 1}]
		assert result.counts == expected
	
class TestCnot(unittest.TestCase):
	
	def testX(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = {Pauli.X: CountingNoiseModelX()}
		code = TrivialStablizerCode()
		
		cnot = TransCnot(kGood, code, code)
		result = cnot.count(noise, Pauli.X)
		
		expected = [{(0, 0): 1}, {(0, 1): 1, (1, 0): 1, (1, 1): 1}]
		assert result.counts == expected
		

class TestBellPair(unittest.TestCase):
	
	@staticmethod
	def BellPair(kGood, code):
		plus = TestPrep.Plus(kGood, code)
		zero = TestPrep.Zero(kGood, code)
		return BellPair(kGood, plus, zero, kGood)				
	
	def testX(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = {Pauli.X: CountingNoiseModelX()}
		code = TrivialStablizerCode()
		
		bp = self.BellPair(kGood, code)
		result = bp.count(noise, Pauli.X)
		expected = [{(0, 0): 1}, {(0, 1): 2, (1, 0): 1, (1, 1): 1}]
		assert result.counts == expected

class TestBellMeas(unittest.TestCase):
	
	@staticmethod
	def BellMeas(kGood, code):
		return BellMeas(kGood, code)
	
	def testX(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = {Pauli.X: CountingNoiseModelX()}
		code = TrivialStablizerCode()
		
		bm = self.BellMeas(kGood, code)
		result = bm.count(noise, Pauli.X)
		expected = [{(0, 0): 1}, {(0, 0): 1, (0, 1): 3}]
		assert result.counts == expected
	
class TestTeleportED(unittest.TestCase):
	
	@staticmethod
	def TeleportED(kGood, code):
		bellPair = TestBellPair.BellPair(kGood, code)
		bellMeas = TestBellMeas.BellMeas(kGood, code)
		return TeleportED(kGood, bellPair, bellMeas)
		
class TestExRec(unittest.TestCase):

	@unittest.skip('')
	def testX(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = {Pauli.X: CountingNoiseModelX()}
		code = TrivialStablizerCode()
		empty = Empty(code) 
		
		ec = TestTeleportED.TeleportED(kGood, code)
		lec = InputAdapter(ec, (0,))
		
		exRec = ExRec(kGood, lec, empty, ec)
		result = exRec.count(noise, Pauli.X)
		print result.counts
		
		#print result.counts

#class TestBellPair(unittest.TestCase):
#
#
#	def testGolay(self):
#		kGood = {Pauli.X: 1}
#		noise = { Pauli.X: CountingNoiseModel(),
#		  #Pauli.Z: NoiseModelZSympy(),
#		  #Pauli.Y: NoiseModelXZSympy() 
#		 }
#		
#		code = TrivialStablizerCode()
#		plus = Prep(kGood, Locations([counterUtils.locXprep('|+>', 0)], '|+>'), StabilizerState(code, [xType]))
#		zero = Prep(kGood, Locations([counterUtils.locZprep('|0>', 0)], '|0>'), StabilizerState(code, [zType]))
#		bellPair = BellPair(kGood, plus, zero, kGood)
#
#		result = bellPair.count(noise)
#		print result.counts
		
#class TestTransCnot(unittest.TestCase):
#
#
#	def testGolay(self):
#		kGood = {Pauli.X: 1}
#		noise = { Pauli.X: CountingNoiseModel(),
#		  #Pauli.Z: NoiseModelZSympy(),
#		  #Pauli.Y: NoiseModelXZSympy() 
#		 }
#		
#		code = TrivialStablizerCode()
#		cnot = TransCnot(kGood, code, code)
#
#		result = cnot.count(noise)
#		print result.counts
		

if __name__ == "__main__":
	

	
	import util.cache
	util.cache.enableFetch(False)
	util.cache.enableMemo(False)
	
	#import sys;sys.argv = ['', 'Test.testName']
	unittest.main()