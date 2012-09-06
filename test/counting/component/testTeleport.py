'''
Created on May 3, 2011

@author: Adam
'''

from counting.block import Block
from counting.component.base import Prep
from counting.component.bell import BellPair, BellMeas
from counting.component.teleport import TeleportED, TeleportEDFilter, TeleportWithMeas, \
	EDInputFilter, Teleport
from counting.key import MultiBlockSyndromeKeyGenerator, SyndromeKeyDecoder
from counting.result import CountResult, TrivialResult
from qec import ed422, error
from qec.error import Pauli, PauliError, xType, zType
from qec.qecc import StabilizerState
from settings import noise
from unittest.case import SkipTest
import counting
import logging
import qec
import testBell
import testComponent
import unittest
from counting.component.adapter import SyndromeAdapter

class TestTeleportWithMeas(testComponent.ComponentTestCase):
	
	@staticmethod
	def TeleportWithMeas(kGood, code):
		bellPair = testBell.TestBellPair.BellPair(kGood, code)
		bellMeas = testBell.TestBellMeas.BellMeas(kGood, code)
		return TeleportWithMeas(kGood, bellPair, bellMeas)
		
	def _getComponent(self, kGood, code):
		return self.TeleportWithMeas(kGood, code)
	
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
	def TeleportED(kGood, code, enableRest=True):
		bellPair = testBell.TestBellPair.BellPair(kGood, code)
		bellMeas = testBell.TestBellMeas.BellMeas(kGood, code)
		return TeleportED(kGood, bellPair, bellMeas, enableRest)
	
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
				
	def testED412(self):
		kGood = {Pauli.Y: 2}
		# Ignore rests, since in Knill's scheme we condition on global acceptance.
		self.code = ed422.ED412Code(gaugeType=error.xType)
		
		prepZ = Prep(kGood, ed422.prepare(Pauli.Z, Pauli.X), StabilizerState(self.code, [error.zType]))
		prepX = Prep(kGood, ed422.prepare(Pauli.X, Pauli.X), StabilizerState(self.code, [error.xType]))
		
		bellPair = BellPair(kGood, prepX, prepZ, kGood)
		bellMeas = BellMeas(kGood, self.code, kGoodCnot=kGood, kGoodMeasX=kGood, kGoodMeasZ=kGood)
		
		bp0 = BellPair({}, prepX, prepZ, {})
		bm0 = BellMeas({}, self.code)
		
		teleportED = TeleportED(kGood, bellPair, bellMeas, enableRest=False)
		ed0 = TeleportED(kGood, bp0, bm0, enableRest=False)
		edFilter = EDInputFilter(teleportED)
		
		noiseModels = {Pauli.Y: noise.NoiseModelXZSympy()}
		noiseModels = self.countingNoiseModels
		
		result = teleportED.count(noiseModels, Pauli.Y)
		print '[[4,1,2]] counts:'
		print 'A=', len(teleportED.locations()), [l['type'] for l in teleportED.locations()]
		print result.counts
		print [sum(count.values()) for count in result.counts]
		
		decoder412 = SyndromeKeyDecoder(ed422.ED412Code(gaugeType=error.xType))
		decoderTrivial = SyndromeKeyDecoder(self.trivialCode)
		decodedCounts = []
		for count in result.counts:
			decoded = {}
			for key, val in count.iteritems():
				d = decoderTrivial.asPauli(decoder412.decode(key[0]))
				decoded[d] = decoded.get(d, 0) + val
			decodedCounts.append(decoded)
		print decodedCounts
		
		result = edFilter.count(noiseModels, Pauli.Y, result, kMax=2)
		print '[[4,1,2]] counts after TED filter:'
		print result.counts
		print [sum(count.values()) for count in result.counts]
		
		decoder412 = SyndromeKeyDecoder(ed422.ED412Code(gaugeType=error.xType))
		decoderTrivial = SyndromeKeyDecoder(self.trivialCode)
		decodedCounts = []
		for count in result.counts:
			decoded = {}
			for key, val in count.iteritems():
				d = decoderTrivial.asPauli(decoder412.decode(key[0]))
				decoded[d] = decoded.get(d, 0) + val
			decodedCounts.append(decoded)
		print decodedCounts
				
	def _getComponent(self, kGood, code, enableRest=True):
		return self.TeleportED(kGood, code, enableRest)		

class TestTeleport(testComponent.ComponentTestCase):
	
	@staticmethod
	def Teleport(kGood, code, enableRest=True):
		bellPair = testBell.TestBellPair.BellPair(kGood, code)
		bellMeas = testBell.TestBellMeas.BellMeas(kGood, code)
		return Teleport(kGood, bellPair, bellMeas, enableRest)
	
					
	def testED412(self):
		kGood = {Pauli.Y: 1}
		# Ignore rests, since in Knill's scheme we condition on global acceptance.
		self.code = ed422.ED412Code(gaugeType=error.xType)
		
		prepZ = Prep(kGood, ed422.prepare(Pauli.Z, Pauli.X), StabilizerState(self.code, [error.zType]))
		prepX = Prep(kGood, ed422.prepare(Pauli.X, Pauli.X), StabilizerState(self.code, [error.xType]))
		
		bellPair = BellPair(kGood, prepX, prepZ, kGood)
		bellMeas = BellMeas(kGood, self.code, kGoodCnot=kGood, kGoodMeasX=kGood, kGoodMeasZ=kGood)
		
		bp0 = BellPair({}, prepX, prepZ, {})
		bm0 = BellMeas({}, self.code)
		
		teleport = SyndromeAdapter(Teleport(kGood, bellPair, bellMeas, enableRest=False))
		
		noiseModels = {Pauli.Y: noise.NoiseModelXZSympy()}
#		noiseModels = self.countingNoiseModels
		
		result = teleport.count(noiseModels, Pauli.Y)
		print '[[4,1,2]] Teleport counts:'
		print 'A=', len(teleport.locations()), [l['type'] for l in teleport.locations()]
		print result.counts
		print [sum(count.values()) for count in result.counts]
		
		decoder412 = SyndromeKeyDecoder(ed422.ED412Code(gaugeType=error.xType))
		decoderTrivial = SyndromeKeyDecoder(self.trivialCode)
		decodedCounts = []
		for count in result.counts:
			decoded = {}
			for key, val in count.iteritems():
				d = decoderTrivial.asPauli(decoder412.decode(key[0]))
				decoded[d] = decoded.get(d, 0) + val
			decodedCounts.append(decoded)
		print decodedCounts
		
		result = prepZ.count(noiseModels, Pauli.Y)
		print '[[4,1,2]] Prep |0> counts:'
		print 'A=', len(prepZ.locations()), [l['type'] for l in prepZ.locations()]
		print result.counts
		print [sum(count.values()) for count in result.counts]
		
		result = prepX.count(noiseModels, Pauli.Y)
		print '[[4,1,2]] Prep |+> counts:'
		print 'A=', len(prepX.locations()), [l['type'] for l in prepX.locations()]
		print result.counts
		print [sum(count.values()) for count in result.counts]
				
	def _getComponent(self, kGood, code, enableRest=True):
		return self.Teleport(kGood, code, enableRest)


#class TestTeleportEDFilter(testComponent.ComponentTestCase):
#	
#	def testCount(self):
#		kGood = {Pauli.X: 2, Pauli.Z: 2}
#		
#		code = ed422.ED412Code()
#		Xl = code.logicalOperators()[0][xType]
#		Zl = code.logicalOperators()[0][zType]
#		
#		teleportED = self._getComponent(kGood, code)
#		
#		# These errors should be undetected, and pass through
#		III = {0: Pauli.I, 1: Pauli.I, 2: Pauli.I}
#		IIX = {0: Pauli.I, 1: Pauli.I, 2: Pauli.X}
#		IXlX = {0: Pauli.I, 1: Xl, 2: Pauli.X}
#		ZlIX = {0: Zl, 1: Pauli.I, 2: Pauli.X}
#		
#		# These errors should be detected, and not pass through.
#		IXI = {0: Pauli.I, 1: Pauli.X, 2: Pauli.I}
#		ZII = {0: Pauli.Z, 1: Pauli.I, 2: Pauli.I}
#		
#		blocks = tuple(Block(i, code) for i in range(3))
#		
#		generator = MultiBlockSyndromeKeyGenerator(blocks)
#		counts = {}
#		expected = {}
#		notDetected = set([Pauli.I, Xl, Zl])
#		for i, err in enumerate([III, IIX, IXlX, ZlIX, IXI, ZII]):
#			key = generator.get_key(err)
#			counts[key] = i+1
#			
#			if err[0] in notDetected and err[1] in notDetected:
#				k = (key[2],)
#				expected[k] = expected.get(k, 0) + i+1
#		
#		expected = [expected]
#		
#		inputResult = CountResult([counts], teleportED.inBlocks())
#		result = teleportED.count(self.countingNoiseModels, Pauli.Y, inputResult=inputResult)
##		print result.counts, expected
#		assert expected == result.counts
#	
#	def _getComponent(self, kGood, code):
#		teleport = TeleportEDFilter(TestTeleportWithMeas.TeleportWithMeas(kGood, code))
#		return teleport
#	
#class TestEDInputFilter(testComponent.ComponentTestCase):
#	
#	
#	def testCount(self):
#		kGood = {Pauli.X: 2, Pauli.Z: 2}
#		teleportED = self._getComponent(kGood, self.trivialCode)
#		
#
#		counts = {(2,):1}
#		inputResult = CountResult([counts], teleportED.inBlocks())
#		
#		for pauli in (Pauli.X, Pauli.Z):
#			result = teleportED.count(self.countingNoiseModels, pauli, inputResult)
#			expected = [{(2,): 1}, {(2,): 9}, {(2,): 30}]
#			print result.counts
#			assert expected == result.counts
#	
#	def _getComponent(self, kGood, code):
#		edif = EDInputFilter(TestTeleportED.TeleportED(kGood, code))
#		return edif
#	
if __name__ == "__main__":

	import util.cache
	util.cache.enableFetch(False)
	util.cache.enableMemo(False)

	complog = logging.getLogger('counting.component')
	complog.setLevel(logging.DEBUG)
	
	unittest.main()