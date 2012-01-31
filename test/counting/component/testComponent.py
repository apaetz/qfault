'''
Created on May 3, 2011

@author: Adam
'''

from counting import probability
from counting.result import TrivialResult, CountResult
from qec.error import Pauli
from qec.qecc import TrivialStablizerCode
from settings.noise import NoiseModelXSympy, NoiseModelZSympy, NoiseModelXZSympy, \
	CountingNoiseModelX, CountingNoiseModelZ
import logging
import unittest
from counting.block import Block



class ComponentTestCase(unittest.TestCase):
	
	depolarizingNoiseModels = {Pauli.X: NoiseModelXSympy(),
			 			  Pauli.Z: NoiseModelZSympy(),
			 			  Pauli.Y: NoiseModelXZSympy()}
	
	countingNoiseModels = {Pauli.X: CountingNoiseModelX(),
						   Pauli.Z: CountingNoiseModelZ()}

	trivialCode = TrivialStablizerCode()
	
	def setUp(self):
		import util.cache
		util.cache.enableFetch(False)
		#util.cache.enableMemo(False)
	
		complog = logging.getLogger('counting.component')
		#complog.setLevel(logging.DEBUG)
		#complog.setLevel(logging.INFO)
	
	def testProbabilities(self):
		for k in range(3):
			for pauli in (Pauli.X, Pauli.Z, Pauli.Y):
				kGood = {pauli: k}
				component = self._getComponent(kGood, self.trivialCode)
				
				inputResult = TrivialResult(component.inBlocks())
				result = component.count(self.depolarizingNoiseModels, pauli, inputResult=inputResult)
				prBad = component.prBad(self.depolarizingNoiseModels[pauli], pauli)
				#prAccept = component.prAccept(self.depolarizingNoiseModels)
				pr = probability.countsToPoly(result.counts, component.locations().getTotals(), self.depolarizingNoiseModels[pauli])
				#prTot = pr/prAccept + prBad
				
				for p in (0, 0.0001, 0.001, 0.01, 0.1, 1):
					gamma = p/15
					try:
						#assert .9999999 <= prTot(gamma)
						#assert 1 >= prAccept(gamma)
						assert 0 <= prBad(gamma)
					except AssertionError as e:
						print 'k=', k
						print 'counts=', result.counts
						print 'pr=', pr
						#print 'prAccept=', prAccept
						print 'prBad=', prBad
						print 'pr({0})={1}'.format(gamma, pr(gamma))
						#print 'prAccept({0})={1}'.format(gamma, prAccept(gamma))
						print 'prBad({0})={1}'.format(gamma, prBad(gamma))
						#print 'prTot({0})={1}'.format(gamma, prTot(gamma))
						raise e
	
	def testInBlocks(self):
		blocks = self._getComponent(0, self.trivialCode).inBlocks()
		assert type(blocks) == tuple
		assert all(type(block) == Block for block in blocks)					

	def testOutBlocks(self):
		blocks = self._getComponent(0, self.trivialCode).outBlocks()
		assert type(blocks) == tuple
		assert all(type(block) == Block for block in blocks)				

	def _getComponent(self, kGood, code):
		raise NotImplementedError


def trivialInput(component):
	inputs = tuple([0]*len(component.inBlocks()))
	inputCounts = [{inputs: 1}]
	inputResult = CountResult(inputCounts, component.inBlocks())
	return inputResult