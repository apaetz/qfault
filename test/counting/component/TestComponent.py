'''
Created on May 3, 2011

@author: Adam
'''

import logging
from counting.component.ec import TECDecodeAdapter, ConcatenatedTEC
from unittest.case import SkipTest
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
	CountingNoiseModelX, CountingNoiseModelZ, NoiseModelZSympy, NoiseModelXZSympy
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
	
	noise = {Pauli.X: CountingNoiseModelX(),
			 Pauli.Z: CountingNoiseModelZ()}
	code = TrivialStablizerCode()
	
	@staticmethod
	def Zero(kGood, code):
		return Prep(kGood, Locations([counterUtils.locZprep('|0>', 0)], '|0>'), StabilizerState(code, [zType]))
	
	@staticmethod
	def Plus(kGood, code):
		return Prep(kGood, Locations([counterUtils.locXprep('|+>', 0)], '|+>'), StabilizerState(code, [xType]))
	
	def testCount(self):
		expectedZero = [{(0,): 1}, {(1,): 1}]
		expectedPlus = [{(0,): 1}, {}]
		
		self.doCountTest(Pauli.X, expectedPlus, expectedZero)
		
		expectedPlus = [{(0,): 1}, {(2,): 1}]
		expectedZero = [{(0,): 1}, {}]
		self.doCountTest(Pauli.Z, expectedPlus, expectedZero)
		
	def doCountTest(self, pauli, expectedPlus, expectedZero):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		
		plus = self.Plus(kGood, self.code)
		result = plus.count(self.noise, pauli)
		assert result.counts == expectedPlus
		
		zero = self.Zero(kGood, self.code)
		result = zero.count(self.noise, pauli)
		assert result.counts == expectedZero
		
	def testPrBad(self):
		logger = logging.getLogger('counting.probability')
		logger.setLevel(logging.DEBUG)
		
		kGood = {Pauli.X: 0, Pauli.Z: 1}
		pauli = Pauli.Z
		
		plus = self.Plus(kGood, self.code)
		
		# In the Z-error case we are counting all (1) locations
		# so there are no bad cases.
		prBad = plus.prBad(self.noise[Pauli.Z], pauli, 1)
		assert 0 == prBad
		
		# In the X-error case, |+> does not cause any
		# generate X-errors, so there are again no bad cases.
		prBad = plus.prBad(self.noise[Pauli.X], pauli, 1)
		assert 0 == prBad
		
		# Now we count zero of the (1) locations.  So
		# everything is bad.
		kGood[Pauli.Z] = 0
		plus = self.Plus(kGood, self.code)
		prBad = plus.prBad(self.noise[Pauli.Z], pauli, 1)
		assert 1 == prBad
		
	def testPrAccept(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		plus = self.Plus(kGood, self.code)
		zero = self.Zero(kGood, self.code)
		
		paulis = (Pauli.X, Pauli.Z)
		assert all((1 == plus.prAccept(self.noise, 1) for pauli in paulis))
		assert all((1 == zero.prAccept(self.noise, 1) for pauli in paulis))
		

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
	
	noise = {Pauli.X: NoiseModelXSympy(),
			 Pauli.Z: NoiseModelZSympy(),
			 Pauli.Y: NoiseModelXZSympy()}
	code = TrivialStablizerCode()
	
	@staticmethod
	def TeleportED(kGood, code):
		bellPair = TestBellPair.BellPair(kGood, code)
		bellMeas = TestBellMeas.BellMeas(kGood, code)
		return TeleportED(kGood, bellPair, bellMeas)
	
	def testPrAccept(self):
		kGood = {Pauli.Y: 6}
		teleportED = self.TeleportED(kGood, self.code)
		teleportED = InputAdapter(teleportED, (0,))
		prAccept = teleportED.prAccept(self.noise)
		
		# With a trivial code, i.e., no code, the acceptance
		# probability should approach 1 as k increases.
		print prAccept
		print prAccept(0), prAccept(0.001)
		
	def testPrBad(self):
		kGood = {Pauli.X: 0}
		teleportED = self.TeleportED(kGood, self.code)
		prBad = teleportED.prBad(self.noise[Pauli.X], Pauli.X, kMax=1)
		print prBad
		print prBad(0), prBad(0.001)
				

class TestExRec(unittest.TestCase):

	def testX(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = {Pauli.X: CountingNoiseModelX()}
		code = TrivialStablizerCode()
		empty = Empty(code) 
		
		ec = TestTeleportED.TeleportED(kGood, code)
		lec = InputAdapter(ec, (0,))
		tec = TECDecodeAdapter(ec)
		
		exRec = ExRec(kGood, lec, empty, tec)
		result = exRec.count(noise, Pauli.X)
		expected = [{Pauli.I: 1}, {Pauli.X: 3, Pauli.I: 5}]
		print result.counts
		assert result.counts == expected
		

	def testXCnot(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = {Pauli.X: CountingNoiseModelX()}
		code = TrivialStablizerCode()
		cnot = TransCnot(kGood, code, code) 
		
		ec = TestTeleportED.TeleportED(kGood, code)
		lec = InputAdapter(ec, (0,))
		lec = ConcatenatedComponent(kGood, lec, lec)
		tec = TECDecodeAdapter(ec)
		tec = ConcatenatedTEC(kGood, tec, tec)
					
		exRec = ExRec(kGood, lec, cnot, tec)
		result = exRec.count(noise, Pauli.X)
		expected = [{Pauli.I: 1}, {Pauli.X: 2, Pauli.I: 4}]
		#print 'cnot exRec:' + str(result.counts)# == expected
		
	@SkipTest
	def testPrAccept(self):
		kGood = {Pauli.X: 1, Pauli.Z: 1}
		noise = {Pauli.X: CountingNoiseModelX()}
		code = TrivialStablizerCode()
		cnot = TransCnot(kGood, code, code) 
		
		ec = TestTeleportED.TeleportED(kGood, code)
		lec = InputAdapter(ec, (0,))
		lec = ConcatenatedComponent(kGood, lec, lec)
		tec = TECDecodeAdapter(ec)
		tec = ConcatenatedTEC(kGood, tec, tec)

		exRec = ExRec(kGood, lec, cnot, tec)
		pr = exRec.prAccept(noise)
		print pr
		print pr(0.000)
		print 'pr(0.001)', pr(0.001)

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

	complog = logging.getLogger('counting.component')
	complog.setLevel(logging.DEBUG)
	
	
	
	#import sys;sys.argv = ['', 'Test.testName']
	suite = unittest.TestSuite()
	suite.addTest(TestTeleportED('testPrAccept'))
	suite.addTest(TestTeleportED('testPrBad'))
	unittest.TextTestRunner().run(suite)
	unittest.main()