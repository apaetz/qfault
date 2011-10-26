'''
Created on May 3, 2011

@author: Adam
'''

import logging
logging.basicConfig(level=logging.DEBUG)

from counting.component.base import Prep, Empty
from counting.component.bell import BellPair, BellMeas
from counting.component.teleport import UncorrectedTeleport, TeleportED
from counting.location import Locations
from qec.error import Pauli, xType, zType
from qec.qecc import TrivialStablizerCode, StabilizerState
from settings.noise import NoiseModelXSympy, CountingNoiseModel, \
	CountingNoiseModelX, CountingNoiseModelZ
from util import counterUtils
import golay
import unittest
from counting.component.rectangle import CnotRectangle
from counting.key import keyConcatenator
from counting.countErrors import mapCounts





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
		
		
class TestCnotRectangle(unittest.TestCase):


	def testX(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = {Pauli.X: CountingNoiseModelX()}
		
		code = TrivialStablizerCode()
		plus = Prep(kGood, Locations([counterUtils.locXprep('|+>', 0)], '|+>'), StabilizerState(code, [xType]))
		zero = Prep(kGood, Locations([counterUtils.locZprep('|0>', 0)], '|0>'), StabilizerState(code, [zType]))
		bellPair = BellPair(kGood, plus, zero, kGood)
		bellMeas = BellMeas(kGood, code)
		data = Empty(code).count(noise, Pauli.X)
		concatenator, data.keyMeta = keyConcatenator(data.keyMeta, data.keyMeta)
		concat = lambda key: concatenator(key, key)
		data.counts = mapCounts(data.counts, concat)
		
		tec = TeleportED(kGood, bellPair, bellMeas)
		rect = CnotRectangle(kGood, kGood, tec)
		result = rect.count(noise, Pauli.X, data)
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
	
	#import sys;sys.argv = ['', 'Test.testName']
	unittest.main()