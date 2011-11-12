'''
Created on May 3, 2011

@author: Adam
'''

from counting import probability
from qec.error import Pauli
from qec.qecc import TrivialStablizerCode
from settings.noise import NoiseModelXSympy, NoiseModelZSympy, \
	NoiseModelXZLowerSympy, NoiseModelXZUpperSympy, CountingNoiseModelX, \
	CountingNoiseModelZ
import logging
import unittest



class ComponentTestCase(unittest.TestCase):
	
	depolarizingNoiseModels = {Pauli.X: NoiseModelXSympy(),
			 			  Pauli.Z: NoiseModelZSympy(),
			 			  Pauli.Y: NoiseModelXZUpperSympy()}
	
	countingNoiseModels = {Pauli.X: CountingNoiseModelX(),
						   Pauli.Z: CountingNoiseModelZ()}

	trivialCode = TrivialStablizerCode()
	
	def setUp(self):
		import util.cache
		util.cache.enableFetch(False)
		#util.cache.enableMemo(False)
	
		complog = logging.getLogger('counting.component')
		#complog.setLevel(logging.DEBUG)
		complog.setLevel(logging.INFO)
	
	def testProbabilities(self):
		for k in range(3):
			for pauli in (Pauli.X, Pauli.Z, Pauli.Y):
				kGood = {pauli: k}
				component = self._getComponent(kGood, self.trivialCode)
				
				result = component.count(self.depolarizingNoiseModels, pauli)
				prBad = component.prBad(self.depolarizingNoiseModels[pauli], pauli)
				prAccept = component.prAccept(self.depolarizingNoiseModels)
				pr = probability.countsToPoly(result.counts, component.locations().getTotals(), self.depolarizingNoiseModels[pauli])
				prTot = pr/prAccept + prBad
				
				for p in (0, 0.0001, 0.001, 0.01, 0.1, 1):
					gamma = p/15
					try:
						assert .9999999 <= prTot(gamma)
						assert 1 >= prAccept(gamma)
						assert 0 <= prBad(gamma)
					except AssertionError as e:
						print 'k=', k
						print 'counts=', result.counts
						print 'pr=', pr
						print 'prAccept=', prAccept
						print 'prBad=', prBad
						print 'pr({0})={1}'.format(gamma, pr(gamma))
						print 'prAccept({0})={1}'.format(gamma, prAccept(gamma))
						print 'prBad({0})={1}'.format(gamma, prBad(gamma))
						raise e
						
				

	def _getComponent(self, kGood, code):
		raise NotImplementedError
