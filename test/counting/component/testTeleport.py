'''
Created on May 3, 2011

@author: Adam
'''

from counting.component.adapter import InputAdapter
from counting.component.teleport import TeleportED, UncorrectedTeleport,\
	TeleportEDFilter
from qec.error import Pauli, PauliError, xType, zType
from unittest.case import SkipTest
import logging
import testComponent
import testBell
import unittest
from counting.key import MultiBlockSyndromeKeyGenerator
from counting.block import Block
from qec import ed422
from counting.result import CountResult
import counting
import qec


#class TestUncorrectedTeleport(testComponent.ComponentTestCase):
#	
#	def _getComponent(self, kGood, code):
#		bellPair = testBell.TestBellPair.BellPair(kGood, code)
#		bellMeas = testBell.TestBellMeas.BellMeas(kGood, code)
#		return UncorrectedTeleport(kGood, bellPair, bellMeas)
#	
#	def testCount(self):
#		kGood = {Pauli.X: 1, Pauli.Z: 1}
#		teleport = self._getComponent(kGood, self.trivialCode)
#		
#		expected = {Pauli.X: [{(0, 0, 0): 1}, {(0, 1, 1): 1, (0, 1, 0): 4, (0, 0, 0): 1, (0, 0, 1): 3}],
#				    Pauli.Z: [{(0, 0, 0): 1}, {(2, 0, 2): 1, (2, 0, 0): 5, (0, 0, 0): 1, (0, 0, 2): 2}]}
#		
#		for pauli in (Pauli.X, Pauli.Z):
#			result = teleport.count(self.countingNoiseModels, pauli)
#			print result.counts
#			assert expected[pauli] == result.counts
#		
	
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
			
	def testPostselect(self):
		'''
		Check that error-detection and Pauli frame updates work correctly.
		'''
		code = ed422.ED412Code()
		teleportED = self.TeleportED({}, code)
		
		blocks = [Block(str(i), code) for i in range(3)]
		keyInGen = MultiBlockSyndromeKeyGenerator(blocks)
		keyOutGen = MultiBlockSyndromeKeyGenerator([blocks[2]])
		logical = code.logicalOperators()[0]
		
		# Check X and Z errors separately.
		# Each possible error on the relevant measurement block is tested
		# individually.
		for pauliType in (qec.error.xType, qec.error.zType):
			dualType = qec.error.dualType(pauliType)
			measBlockNum = 0
			if xType == pauliType: measBlockNum = 1
			
			for e in range(1<<4):
				# Construct the Pauli error
				ebits = {qec.error.xType: 0, qec.error.zType: 0}
				ebits[pauliType] = e
				pauliError = PauliError(xbits=ebits[xType], zbits=ebits[zType])
				
				syndrome = code.getSyndrome(pauliError)
				errors = {str(measBlockNum): pauliError, '2': logical[dualType]}
				keyIn = keyInGen.getKey(errors)
				
				# Construct the expected output key (error).  It should be the
				# same as block 2 of the input, unless a logical correction (Pauli
				# frame update) is required.
				if pauliError.commutesWith(logical[dualType]):
					keyOut = keyOutGen.getKey({'2': logical[dualType]})
				else:
					keyOut = keyOutGen.getKey({'2': logical[pauliType] * logical[dualType]})
				
				result = CountResult([{keyIn: 42}], keyInGen.keyMeta(), blocks)
				postResult = teleportED._postCount(result, {}, Pauli.X)
				
#				print 'keyIn=', keyIn
#				print 'keyOut=', keyOut
#				print 'accepted=',postResult.counts
#				print 'rejected=',postResult.rejected
				if syndrome:
					assert [{}] == postResult.counts
					assert [{counting.key.rejectKey: 42}] == postResult.rejected
				else:
					assert [{counting.key.rejectKey: 0}] == postResult.rejected
					assert [{keyOut: 42}] == postResult.counts

		
		
	
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
	
class TestTeleportEDFilter(testComponent.ComponentTestCase):
	
	def testCount(self):
		kGood = {Pauli.X: 2, Pauli.Z: 2}
			
		for inSyndrome in (0, 1, 2):
			teleportED = self._getComponent(kGood, self.trivialCode, inputSyndrome=inSyndrome)
			
			for pauli in (Pauli.X, Pauli.Z):
				result = teleportED.count(self.countingNoiseModels, pauli)
				expected = [{(inSyndrome,): 1}, {(inSyndrome,): 9}, {(inSyndrome,): 30}]
				print result.counts
				assert expected == result.counts
	
	def _getComponent(self, kGood, code, inputSyndrome=0):
		bellPair = testBell.TestBellPair.BellPair(kGood, code)
		bellMeas = testBell.TestBellMeas.BellMeas(kGood, code)
		teleport = TeleportEDFilter(kGood, bellPair, bellMeas)
		return InputAdapter(teleport, (inputSyndrome,))
	
if __name__ == "__main__":

	import util.cache
	util.cache.enableFetch(False)
	util.cache.enableMemo(False)

	complog = logging.getLogger('counting.component')
	complog.setLevel(logging.DEBUG)
	
	unittest.main()