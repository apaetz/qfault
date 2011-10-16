'''
Created on May 1, 2011

@author: Adam
'''
from abc import abstractmethod, ABCMeta
from block import CountedBlock
from util.cache import fetchable
from countErrors import filterAndPropagate, convolveABB
from qec.error import Pauli
from counting.countErrors import countErrors, countBlocksBySyndrome
from util import counterUtils
import operator
from counting.location import Locations
from counting.block import Block

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

	def __init__(self, kGood, kBest=0, nickname=None, subcomponents={}):
		'''
		Constructor
		'''
		self.kGood = kGood
		self.kBest = kBest
		self._nickname = nickname
		self._subs = subcomponents
		
	def __set__(self, name, component):
		self._subs[name] = component
	
	def __get__(self, name):
		return self._subs[name]
	
	def locations(self):
		'''
		Returns the set of locations contained in the component (and all sub-components).
		:rtype: :class:`Locations`
		'''
		locs = Locations()
		for sub in self._subs:
			locs += sub.locations()
		return locs + self.internalLocations()
	
	def internalLocations(self):
		'''
		Returns the set of locations contained in the component *excluding* any sub-components.
		:rtype: :class:`Locations`
		'''
		return Locations()
		
	#@fetchable
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
		rep = str(self.__class__.__name__) + '-'
		if None != self._nickname:
			rep += self._nickname
		rep += '-'
		rep += str(self.kGood)
		if 0 != self.kBest:
			rep += str(self.kBest)
		return rep
	
class CountableComponent(Component):
	'''
	This is a component which, in addition to (or instead of) having sub-components,
	also has its own physical locations that must be counted.
	'''
	
	def __init__(self, locations, codes, kGood, kBest=0, nickname=None, subcomponents={}):
		# The number of faulty locations cannot exceed the total
		# number of locations.
		n = len(locations)
		kBest = min(kBest, n)
		kGood = min(kGood, n)
		if None == nickname:
			nickname=str(locations)
		super(CountableComponent, self).__init__(kGood, kBest, nickname=nickname)
		
		self.locations = locations
		blocknames = list(locations.blocknames())
		self.blocks = tuple([Block(name, codes[name]) for name in blocknames])
		
	def internalLocations(self):
		return self.locations
		
	def _count(self, noise):
		# First, count the sub-components.
		subcounts = super(CountableComponent, self)._count(noise)
			
		# Now count the internal locations.
		counts = {}
		for pauli, kGood in [(Pauli.X, self.kGood), (Pauli.Z, self.kGood), (Pauli.X*Pauli.Z, self.kBest)]:
			counts[pauli] = countBlocksBySyndrome(self.locations, self.blocks, pauli, noise[pauli], kGood)
				
		cb = CountedBlock(counts, subblocks=self.blocks, name=self.nickname())
		if 0 == len(subcounts):
			# There are no sub-component counts. So we can just return the
			# internal counts.
			return cb
		
		subcounts[self.nickname()] = cb 
		return subcounts

class Prep(CountableComponent):
	
	def __init__(self, kGood, kBest, locations, code):
		blocknames = list(locations.blocknames())
		super(Prep, self).__init__(locations, {name: code for name in blocknames}, kGood, kBest)
		
class TransCNOT(CountableComponent):
	'''
	Transversal Controlled-NOT.
	'''
	
	def __init__(self, kGood, kBest, ctrlCode, targCode):
		n = ctrlCode.blockLength()
		if n != targCode.blockLength():
			raise Exception('Control ({0}) and target ({1}) blocklengths do not match.'.format(n, targCode.blockLength()))
		
		nickname='transCNOT.'+str(n)
		print 'nickname=', nickname
		locs = Locations([counterUtils.loccnot('ctrl', i, 'targ', i) for i in range(n)], nickname)
		codes = {'ctrl': ctrlCode, 'targ': targCode}
		super(TransCNOT, self).__init__(locs, codes, kGood, kBest, nickname)
		
class Bell(Component):
	
	def __init__(self, kGood, kBest, plus, zero, code):
		super(Bell, self).__init__(kGood, kBest, subcomponents={'|+>': plus, '|0>': zero})
		self._code = code
		
	def _count(self, noise):
		prepCounts = super(Bell, self)._count(noise)
		
		# Sanity check.  Both of the input blocks must be states of the
		# 
		prepCodes = [block.getCode() for block in prepCounts]
		if not all(issubclass(code, self._code) for code in prepCodes):
			raise Exception('One of {0} is not a state of {1}'.format(prepCodes, self._code))
	
		
		
		
		
class VerifyX(Component):
	
	cmName = 'cnotMeas'
	
	def __init__(self, kGood, kBest, zeroPrepA, zeroPrepB):
		cnotMeas = 0 # TODO
		subcomponents = [zeroPrepA, zeroPrepB, cnotMeas]
		super(VerifyX, self).__init__(kGood, kBest, subcomponents=subcomponents)
		self.prepAName = zeroPrepA.nickname()
		self.prepBName = zeroPrepB.nickname()
		
#	def _convolve(self, counts):
#		results = {}
#		
#		countsA = counts[self.prepAName]
#		countsB = counts[self.prepBName]
#		countsCM = counts[self.cmName]
#		
#		# First, X-error counts.
#		countsBCM = convolve(countsCM[Pauli.X], 
#						     countsB[Pauli.X], 
#						     convolveFcn=convolveABB, 
#						     kMax=self.kGood)
#		
#		countsVerified = convolve(countsA[Pauli.X], 
#								  countsBCM, 
#								  convolveFcn=convolveCountsPostselectX, 
#								  kMax=self.kGood)
#		
#		results[Pauli.X] = countsVerified
#		
#
#		# Now Z errors.
#		
#		# TODO.  Need to reduce CM counts over two blocks to
#		# marginal counts over just block A.  This will require
#		# constructing the 2-block counts a bit differently.
#		# i.e., 2-block counts are indexed by [s1][s2] rather
#		# than a single index for the entire syndrome.
#		countsVerify = convolve(countsPrepB, countsC, kMax=kGood)
#		countsVerified = convolve(countsPrepA, countsVerify, kMax=kGood)
#	
#		
#		
#		return CountedBlock(str(countsA), countsA.getCode, results)
		
		
	