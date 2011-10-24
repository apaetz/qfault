'''
Created on May 3, 2011

@author: Adam
'''
from counting.component import Prep, BellPair, BellMeas, Empty, Teleport,\
	TransCnot
from qec.error import Pauli, xType, zType
from qec.qecc import TrivialStablizerCode, StabilizerState
from settings.noise import NoiseModelXSympy, CountingNoiseModel
import golay
import unittest
from util import counterUtils
from counting.location import Locations


class TestTeleport(unittest.TestCase):


	def testGolay(self):
		kGood = {Pauli.X: 1}
		noise = { Pauli.X: CountingNoiseModel(),
		  #Pauli.Z: NoiseModelZSympy(),
		  #Pauli.Y: NoiseModelXZSympy() 
		 }
		
		code = TrivialStablizerCode()
		plus = Prep(kGood, Locations([counterUtils.locXprep('|+>', 0)], '|+>'), StabilizerState(code, [xType]))
		zero = Prep(kGood, Locations([counterUtils.locZprep('|0>', 0)], '|0>'), StabilizerState(code, [zType]))
		bellPair = BellPair(kGood, plus, zero, kGood)
		bellMeas = BellMeas(kGood, code)
		data = Empty(code, 'data').count(noise)
		
		teleport = Teleport(kGood, data, bellPair, bellMeas)
		
		
		result = teleport.count(noise)
		print result.counts

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
	#import sys;sys.argv = ['', 'Test.testName']
	unittest.main()