'''
Created on May 1, 2011

@author: Adam
'''
from abc import abstractmethod, ABCMeta
from block import CountedBlock
from countErrors import filterAndPropagate, convolveABB
from counting.block import Block
from counting.countErrors import countErrors, countBlocksBySyndrome, mapKeys
from counting.countParallel import convolve
from counting.location import Locations
from qec.error import Pauli, PauliError, xType, zType
from qec.qecc import Codeword
from util import counterUtils, bits
from util.cache import fetchable
import operator
from copy import copy

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

	def __init__(self, kMax, nickname=None, subcomponents={}):
		'''
		Constructor
		'''
		self.kMax = {pauli: kMax.get(pauli, 0) for pauli in (Pauli.X, Pauli.Z, Pauli.Y)}
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
		print 'Counting', str(self)
		countedBlocks = self._count(noise)
		print 'Convolving', str(self)
		block = self._convolve(countedBlocks)
		print 'Post processing', str(self)
		return self._postCount(block)
	
	def _count(self, noise):
		'''
		Subclass hook.
		Counts the errors in each sub-component and returns the counts
		as a dictionary indexed by sub-component name.  
		It is expected that most concrete components will not need to 
		implement this method. 
		'''
		return {name: sub.count(noise) for name, sub in self._subs.iteritems()} 
	
	def _convolve(self, blocks):
		'''
		Subclass hook.
		Combines errors from each of the sub-components.
		The default implementation simply returns 'counts'.
		Any non-trivial component will need to implement this method.
		
		Returns a CountedBlock object.
		'''
		
		# The block names won't be used
		blocks = blocks.values()
		
		if 1 == len(blocks):
			return blocks[0]
		
		convolved = {}			
		counts = [block.counts() for block in blocks]
		for pauli, k in self.kMax.iteritems():
			convolved[pauli] = counts[0][pauli]
			for count in counts[1:]:
				convolved[pauli] = convolve(convolved[pauli], count[pauli], kMax=k)
				
		# Assume that the properties of all of the blocks are the same.
		propBlock = blocks[0]
		return CountedBlock(convolved, propBlock.keyGenerators(), code=propBlock.getCode(), name=str(self))
	
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
		rep = str(self.__class__.__name__)
		if None != self._nickname:
			rep += '-' + self._nickname + '-'
		rep += str(self.kMax)
		return rep
	
class CountableComponent(Component):
	'''
	This is a component which, in addition to (or instead of) having sub-components,
	also has its own physical locations that must be counted.
	'''
	
	def __init__(self, locations, codes, kMax, nickname=None, subcomponents={}):
		# The number of faulty locations cannot exceed the total
		# number of locations.
		if None == nickname:
			nickname=str(locations)
		super(CountableComponent, self).__init__(kMax, nickname=nickname)
		
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
		keyGens = {}
		for pauli,k in self.kMax.iteritems():
			counts[pauli], keyGens[pauli] = countBlocksBySyndrome(self.locations, 
																  self.blocks, 
																  pauli, 
																  noise[pauli], 
																  k)
				
		cb = CountedBlock(counts, keyGens, subblocks=self.blocks, name=self.nickname())
		
		subcounts[self.nickname()] = cb 
		return subcounts

class Prep(CountableComponent):
	
	def __init__(self, kMax, locations, code):
		blocknames = list(locations.blocknames())
		super(Prep, self).__init__(locations, {name: code for name in blocknames}, kMax)
		
class TransCNOT(CountableComponent):
	'''
	Transversal Controlled-NOT.
	'''
	
	def __init__(self, kMax, ctrlCode, targCode):
		n = ctrlCode.blockLength()
		if n != targCode.blockLength():
			raise Exception('Control ({0}) and target ({1}) blocklengths do not match.'.format(n, targCode.blockLength()))
		
		nickname='transCNOT.'+str(n)
		print 'nickname=', nickname
		locs = Locations([counterUtils.loccnot('ctrl', i, 'targ', i) for i in range(n)], nickname)
		codes = {'ctrl': ctrlCode, 'targ': targCode}
		super(TransCNOT, self).__init__(locs, codes, kMax, nickname)
		
	#TODO: the resulting QECC is still not correct.  Need to come up with a multi-block code representation.
	
class TransMeas(CountableComponent):
	
	def __init__(self, kMax, code, basis):
		n = code.blockLength()
		nickname = 'transMeas' + str(basis) + str(n)
		if Pauli.X == basis:
			loc = counterUtils.locXmeas
		elif Pauli.Z == basis:
			loc = counterUtils.locZmeas
		else:
			raise Exception('{0}-basis measurement is not supported'.format(basis))
		
		locs = Locations([loc('0', i) for i in range(n)], nickname)
		super(TransMeas, self).__init__(locs, {'0': code}, kMax, nickname)
			
		
class BellPair(Component):
	
	plusName = '|+>'
	zeroName = '|0>'
	cnotName = 'cnot'
	
	def __init__(self, kMax, plus, zero, kMaxCnot):
		super(BellPair, self).__init__(kMax, subcomponents={self.plusName: plus, self.zeroName: zero})
		self.kMaxCnot = kMaxCnot							
		
	def _count(self, noise):
		prepBlocks = super(BellPair, self)._count(noise)
		
		# Construct a transversal CNOT component from the two input codes.
		cnot = TransCNOT(self.kMaxCnot, prepBlocks[self.plusName].getCode(), prepBlocks[self.zeroName].getCode())
		prepBlocks[self.cnotName] = cnot.count(noise)
		
		return prepBlocks
	
	def _convolve(self, blocks):
		parityChecks = blocks[self.plusName].keyGenerators()
		
		# Masks that identify X stabilizers and operators, and Z stabilizers
		# and operators.
		xmask = bits.listToBits((0 == check[zType]) for check in parityChecks)
		zmask = bits.listToBits((0 == check[xType]) for check in parityChecks)
		
		blockShift = len(parityChecks)

		for pauli in self.kMax.keys(): 	
			# Propagate the preparation errors through the CNOT.
			blocks[self.plusName].counts()[pauli] = mapKeys(blocks[self.plusName].counts()[pauli], 
														    lambda e: (e << blockShift) + (e & xmask))
			blocks[self.zeroName].counts()[pauli] = mapKeys(blocks[self.zeroName].counts()[pauli], 
														    lambda e: ((e & zmask) << blockShift) + e)
			
		# Now convolve.
		convolved = super(BellPair, self)._convolve(blocks)
		
		# The default _convolve() may not assign the correct block properties (i.e., key generators, code).
		cnotBlock = blocks[self.cnotName]
		return CountedBlock(convolved.counts(), cnotBlock.keyGenerators(), code=cnotBlock.getCode(), name=str(convolved))
	
	@staticmethod
	def ConstructBellCode(plusState, zeroState):
		pass
	
class BellMeas(CountableComponent):

	cnotName = 'cnot'
	measXName = 'measX'
	measZName = 'measZ'
	
	def __init__(self, kMax, code, kMaxMeasX=None, kMaxMeasZ=None, kMaxCnot=None):
		if None == kMaxMeasX: kMaxMeasX = kMax
		if None == kMaxMeasZ: kMaxMeasZ = kMax
		if None == kMaxCnot: kMaxCnot = kMax
		
		subs = {self.cnotName: TransCNOT(kMaxCnot, code, code),
			    self.measXName: TransMeas(kMaxMeasX, code, Pauli.X),
			    self.measZName: TransMeas(kMaxMeasZ, code, Pauli.Z)}
		super(BellMeas, self).__init__(kMax, subcomponents=subs)
		
	def _convolve(self, blocks):
		parityChecks = blocks[self.measXName].keyGenerators()
		blockShift = len(parityChecks)
		
		for pauli in self.kMax.keys():
			# Extend single block X-basis measurement keys to both blocks.
			# Z-basis measurements do not need to be modified because they are represented by the LSBs.
			blocks[self.measXName].counts()[pauli] = mapKeys(blocks[self.measXName].counts()[pauli], 
															 lambda e: e << blockShift)
			
		# Now convolve.
		convolved = super(BellPair, self)._convolve(blocks)
		
		# The default _convolve() may not assign the correct block properties (i.e., key generators, code).
		cnotBlock = blocks[self.cnotName]
		return CountedBlock(convolved.counts(), cnotBlock.keyGenerators(), code=cnotBlock.getCode(), name=str(convolved))
			
		
	
class Teleport(Component):
	
	bpName = 'BP'
	bmName = 'BM'
	
	def __init__(self, kGood, kBest, data, bellPair, bellMeas):
		subs = {self.bpName: bellPair,
				self.bmName: bellMeas}
		
		super(Teleport, self).__init__(kGood, kBest, subcomponents=subs)
		self._data = data
		
	def _convolve(self, blocks):
		#TODO
		pass
		
	
		
		
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
#						     kMax=self.kMax)
#		
#		countsVerified = convolve(countsA[Pauli.X], 
#								  countsBCM, 
#								  convolveFcn=convolveCountsPostselectX, 
#								  kMax=self.kMax)
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
#		countsVerify = convolve(countsPrepB, countsC, kMax=kMax)
#		countsVerified = convolve(countsPrepA, countsVerify, kMax=kMax)
#	
#		
#		
#		return CountedBlock(str(countsA), countsA.getCode, results)
		
		
	