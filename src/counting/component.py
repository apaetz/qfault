'''
Created on May 1, 2011

@author: Adam
'''
from abc import abstractmethod, ABCMeta
from block import CountedBlock
from util.cache import fetchable
from countErrors import propagateAndReduceZero, convolveABB
from qec.error import Pauli
from counting.countErrors import countErrors
from util import counterUtils
import operator

class Component(object):
	'''
	Abstract class for circuit "components" of the extended rectangle (exRec).
	To count errors, the exRec is divided up into a hierarchy of components.
	
	Each component contains a number of sub-components.  Errors are counted
	recursively by first counting small numbers of errors in each sub-component
	and then combining the error information.
	
	A component can be treated as a black box which takes as input a maximum
	number of faults k, and outputs a set of counts, k counts for each possible
	error that can occur inside of the component.
	
	The primary public method of interest is count().  This method counts errors
	in the component in three steps:
	  1. Errors in the sub-components are counted.
	  2. Errors from each sub-component are combined.
	  3. Optional post-processing is performed.
	  
	The count() method is a template.  So rather than overriding count(),
	subclasses should instead implement hook methods _count(), _convolve(),
	and _postCount().
	'''

	def __init__(self, kGood, kBest=0, nickname=None, subcomponents=[]):
		'''
		Constructor
		'''
		self.kGood = kGood
		self.kBest = kBest
		self._nickname = nickname
		self._subs = {}
		
	def __set__(self, name, component):
		self._subs[name] = component
	
	def __get__(self, name):
		return self._subs[name]
		
	@fetchable
	def count(self, noise):
		'''
		Counts errors in the component.
		Returns a CountedBlock. 
		'''
		counts = self._count(noise)
		block = self._convolve(counts)
		return self._postCount(block)
	
	def _count(self, noise):
		'''
		Subclass hook.
		Counts the errors in each sub-component and returns the counts
		as a dictionary indexed by sub-component name.  
		It is expected that most concrete components will not need to 
		implement this method. 
		'''
		subcounts = [(name, sub.count(noise)) for name, sub in self._subs.iteritems()]
		return dict(subcounts)
	
	def _convolve(self, counts):
		'''
		Subclass hook.
		Combines errors from each of the sub-components.
		The default implementation simply returns 'counts'.
		Any non-trivial component will need to implement this method.
		
		Returns a CountedBlock object.
		'''
		return counts
	
	def _postCount(self, block):
		'''
		Subclass hook.
		Performs (optional) error count post-processing.
		
		Returns a CountedBlock object.
		'''
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
			rep += self._nickname
		rep += ')'
		return rep
	
class SyndromeKeyGenerator(object):
	
	def __init__(self, code, paulis):
		self._code = code
		self._paulis = paulis
		
	def getKey(self,e):
		return self._code.getSyndrome(e, self._paulis)
	
	def __repr__(self):
		return str(self._code) + ''.join(str(p) for p in self._paulis)
		
class PrepZero(Component):
	
	def __init__(self, kGood, kBest, locations, code):
		super(PrepZero, self).__init__(kGood, kBest, nickname=str(locations))
		self.locations = locations
		self.code = code
		
	def _count(self, noise):				
		counts = {}
		for pauli, kGood in [(Pauli.X, self.kGood), (Pauli.Z, self.kGood), (Pauli.X*Pauli.Z, self.kBest)]:
			reduced = propagateAndReduceZero(self.locations, self.code, pauli)
			blocknames = reduced.blocknames()
			if 1 != len(blocknames):
				raise Exception('Logical |0> should be only a single block.')
			
			keyGenerator = SyndromeKeyGenerator(self.code, pauli.types())
			counts[pauli] = [countErrors(k, reduced, blocknames, 0, 0, noise[pauli], keyGenerator) 
							 for k in range(kGood+1)]
		return CountedBlock(str(self.locations), self.code, counts)
		
class VerifyX(Component):
	
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
		
		
	