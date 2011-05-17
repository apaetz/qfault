'''
Created on May 1, 2011

@author: Adam
'''
from abc import abstractmethod, ABCMeta
from block import CountedBlock
from cache import fetchable
from countErrors import propagateAndReduceZero, countErrors1Block, convolveABB
from error import Pauli

class Component(object):
	'''
	classdocs
	'''

	def __init__(self, kGood, nickname=None, subcomponents=[]):
		'''
		Constructor
		'''
		self.kGood = kGood
		self._nickname = nickname
		self._subs = {}
		
	def __set__(self, name, component):
		self._subs[name] = component
	
	def __get__(self, name):
		return self._subs[name]
		
	@fetchable
	def count(self, noise):
		counts = self._count(noise)
		block = self._convolve(counts)
		return self._postCount(block)
	
	def _count(self, noise):
		'''
		Subclass hook.
		'''
		subcounts = [(name, sub.count(noise)) for name, sub in self._subs.iteritems()]
		return dict(subcounts)
	
	def _convolve(self, counts):
		'''
		Subclass hook
		'''
		return counts
	
	def _postCount(self, block):
		return block

	def nickname(self):
		return self._nickname
	
	def fullname(self):
		full = str(self.__class__.name)
		if None != self._nickname:
			full += '.' + self._nickname
		
		return full
	
	def __repr__(self):
		rep = str(self.__class__.__name__) + '('
		if None != self._nickname:
			rep += self._name
		rep += ')'
		return rep
	
class BestComponent(Component):
	
	def __init__(self, kGood, kBest, **kwargs):
		super(BestComponent, self).__init__(kGood, **kwargs)
		self.kBest = kBest
	
class PrepZero(BestComponent):
	
	def __init__(self, kGood, kBest, locations, code):
		super(PrepZero, self).__init__(kGood, kBest, nickname=str(locations))
		self.locations = locations
		self.code = code
		
	def _count(self, noise):				
		counts = {}
		for pauli, kGood in [(Pauli.X, self.kGood), (Pauli.Z, self.kGood), (Pauli.X+Pauli.Z, self.kBest)]:
			reduced = propagateAndReduceZero(self.locations, self.code, pauli)
			counts[pauli] = [countErrors1Block(reduced, pauli, noise, k) for k in range(kGood+1)]
		return CountedBlock(str(self.locations), self.code, counts)
		
class VerifyX(BestComponent):
	
	cmName = 'cnotMeas'
	
	def __init__(self, kGood, kBest, zeroPrepA, zeroPrepB):
		cnotMeas = 0 # TODO
		subcomponents = [zeroPrepA, zeroPrepB, cnotMeas]
		super(VerifyX, self).__init__(kGood, kBest, subcomponents=subcomponents)
		self.prepAName = zeroPrepA.nickname()
		self.prepBName = zeroPrepB.nickname()
		
	def _convolve(self, counts):
		results = {}
		
		countsA = counts[self.prepAName]
		countsB = counts[self.prepBName]
		countsCM = counts[self.cmName]
		
		# First, X-error counts.
		countsBCM = convolve(countsCM[Pauli.X], 
						     countsB[Pauli.X], 
						     convolveFcn=convolveABB, 
						     kMax=self.kGood)
		
		countsVerified = convolve(countsA[Pauli.X], 
								  countsBCM, 
								  convolveFcn=convolveCountsPostselectX, 
								  kMax=self.kGood)
		
		results[Pauli.X] = countsVerified
		

		# Now Z errors.
		
		# TODO.  Need to reduce CM counts over two blocks to
		# marginal counts over just block A.  This will require
		# constructing the 2-block counts a bit differently.
		# i.e., 2-block counts are indexed by [s1][s2] rather
		# than a single index for the entire syndrome.
		countsVerify = convolve(countsPrepB, countsC, kMax=kGood)
		countsVerified = convolve(countsPrepA, countsVerify, kMax=kGood)
	
		
		
		return CountedBlock(str(countsA), countsA.getCode, results)
		
		
	