'''
Created on May 1, 2011

@author: Adam
'''
from abc import abstractmethod, ABCMeta
from block import CountedBlock
from counting.block import Block
from counting.countErrors import countBlocksBySyndrome, mapKeys
from counting.countParallel import convolve
from counting.location import Locations
from qec.error import Pauli, xType, zType
from util import counterUtils, bits
from util.cache import fetchable
from counting import probability
from qec.qecc import StabilizerCode


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

	def __init__(self, kGood, nickname=None, subcomponents={}):
		'''
		Constructor
		'''
		self.kGood = {pauli: kGood.get(pauli, 0) for pauli in (Pauli.X, Pauli.Z, Pauli.Y)}
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
	
	def prBad(self, noise, kMax={pauli: None for pauli in [Pauli.X, Pauli.Z, Pauli.Y]}):
		prSelf = {pauli: probability.prBadPoly(self.kGood[pauli], self.internalLocations(), noiseModel, kMax[pauli])
				  for pauli, noiseModel in noise.iteritems()}
		
		prByPauli = [sub.prBad(noise, kMax=self.kGood) for sub in self._subs] + [prSelf]
		
		return {pauli: sum(pr[pauli] for pr in prByPauli) for pauli in noise}
	
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
		for pauli, k in self.kGood.iteritems():
			convolved[pauli] = counts[0][pauli]
			for count in counts[1:]:
				convolved[pauli] = convolve(convolved[pauli], count[pauli], kMax=k)
				
		# Assume that the properties of all of the blocks are the same.
		propBlock = blocks[0]
		return CountedBlock(convolved, propBlock.keyMeta(), code=propBlock.getCode(), name=str(self))
	
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
		rep += str(self.kGood)
		return rep
	
class CountableComponent(Component):
	'''
	This is a component which, in addition to (or instead of) having sub-components,
	also has its own physical locations that must be counted.
	'''
	
	def __init__(self, locations, blocknames, codes, kGood, nickname=None, subcomponents={}):
		# The number of faulty locations cannot exceed the total
		# number of locations.
		if None == nickname:
			nickname=str(locations)
		super(CountableComponent, self).__init__(kGood, nickname=nickname)
		
		self.locations = locations
		self.blocks = tuple([Block(name, codes[name]) for name in blocknames])
		
	def internalLocations(self):
		return self.locations
		
	def _count(self, noise):
		# First, count the sub-components.
		subcounts = super(CountableComponent, self)._count(noise)
			
		# Now count the internal locations.
		counts = {}
		keyMetas = set()
		for pauli,k in self.kGood.iteritems():
			counts[pauli], meta = countBlocksBySyndrome(self.locations, 
														self.blocks, 
														pauli, 
														noise[pauli], 
														k)
			keyMetas.add(meta)
				
		if len(keyMetas) > 1:
			raise Exception('Key metadata mismatch. {0}'.format(keyMetas))
		
		newCode = self._propagateCodes({block.name: block.getCode() for block in self.blocks})
		cb = CountedBlock(counts, keyMetas.pop(), code=newCode, name=self.nickname())
		
		subcounts[self.nickname()] = cb 
		return subcounts
	
	def _propagateCodes(self, codes):
		
		# In the case that there is only a single block, then a sensible
		# default is to simply return the code for that block.  It is not
		# clear how to combine codes from multiple blocks, however.
		if 1 == len(codes):
			return [c for c in codes.values()][0]
		
		raise NotImplementedError
	
class Empty(CountableComponent):
	
	def __init__(self, code, blockname='empty'):
		super(Empty, self).__init__(Locations([], blockname), [blockname], {blockname: code}, {})

class Prep(CountableComponent):
	
	def __init__(self, kGood, locations, code):
		blocknames = list(locations.blocknames())
		super(Prep, self).__init__(locations, blocknames, {name: code for name in blocknames}, kGood)
		
class TransCnot(CountableComponent):
	'''
	Transversal Controlled-NOT.
	'''
	
	ctrlName = 'ctrl'
	targName = 'targ'
	
	def __init__(self, kGood, ctrlCode, targCode, blockorder=[ctrlName, targName]):
		n = ctrlCode.blockLength()
		if n != targCode.blockLength():
			raise Exception('Control ({0}) and target ({1}) blocklengths do not match.'.format(n, targCode.blockLength()))
		
		nickname='transCNOT.'+str(n)
		print 'nickname=', nickname
		locs = Locations([counterUtils.loccnot(self.ctrlName, i, self.targName, i) for i in range(n)], nickname)
		codes = {self.ctrlName: ctrlCode, self.targName: targCode}
		super(TransCnot, self).__init__(locs, blockorder, codes, kGood, nickname)
		
	def _propagateCodes(self, codes):
		code = TransCnotCode(codes)
		return code
		
class TransCnotCode(StabilizerCode):
	
	def __init__(self, codes):
		self._codes = codes
		
		lengths = [code.n for code in codes.values()]
		if not all(lengths[0] == n for n in lengths):
			raise Exception('Codes are not all of the same length.')
		
		subblockLength = lengths[0]
		
		n = sum(code.n for code in codes.values())
		k = sum(code.k for code in codes.values())
		d = min(code.d for code in codes.values())
		name = ''.join(str(code) for code in codes.values())
		super(TransCnotCode, self).__init__(name, n, k, d)
		
		self.subblockLength = subblockLength

	def stabilizerGenerators(self):
		stabs  = {name: self._codes[name].stabilizerGenerators() 
				  for name in [TransCnot.ctrlName, TransCnot.targName]}
		
		# We assume that the transversal CNOT operation is a valid operation
		# in the code, and therefore the stabilizer generators do not change.
		# Additionally, we assume that the result of the CNOT leaves the
		# two blocks unentangled (as in, e.g., error correction).
		# Thus, we need only extend the operators into the larger space.
		
		I = Pauli.I ** self.subblockLength
		
		block0Stabs = [I + stab for stab in stabs[self._codes.keys()[0]]]
		block1Stabs = [stab + I for stab in stabs[self._codes.keys()[1]]]
		
		return tuple(block0Stabs + block1Stabs)
	
	def normalizerGenerators(self):
		ctrlNorms = self._codes[TransCnot.ctrlName].normalizerGenerators()
		targNorms = self._codes[TransCnot.targName].normalizerGenerators()
		
		return self.PropagateOperators(ctrlNorms, targNorms, self.subblockLength, self._codes.keys())
	
	@staticmethod
	def PropagateOperators(ctrlOps, targOps, blocklength, blockorder):
		newOps = []
		
		I = Pauli.I ** blocklength
		
		# For now, assume that the block ordering is [ctrl, targ].
		for stab in ctrlOps:
			if 0 == stab[zType]:
				# This is an X stabilizer. (X -> XX)
				newOps.append(stab+stab)
			else:
				# This is a Z stabilizer. (Z -> ZI)
				newOps.append(stab+I)
				
		for stab in targOps:
			if 0 == stab[zType]:
				# This is an X stabilizer. (X -> IX)
				newOps.append(I+stab)
			else:
				# This is a Z stabilizer. (Z -> ZZ)
				newOps.append(stab+stab)
				
		# Check our original assumption.  Block ordering is
		# big endian, so blockorder[0] is the MSBs.
		if TransCnot.ctrlName != blockorder[0]:
			# Swap the blocks around.
			shift = blocklength
			mask = (1 << shift) - 1
			newOps = [(stab >> shift) ^ ((stab & mask) << shift) for stab in newOps]
				
		return newOps
	
class TransMeas(CountableComponent):
	
	def __init__(self, kGood, code, basis, blockname='0'):
		n = code.blockLength()
		nickname = 'transMeas' + str(basis) + str(n)
		if Pauli.X == basis:
			loc = counterUtils.locXmeas
		elif Pauli.Z == basis:
			loc = counterUtils.locZmeas
		else:
			raise Exception('{0}-basis measurement is not supported'.format(basis))
		
		locs = Locations([loc('0', i) for i in range(n)], nickname)
		super(TransMeas, self).__init__(locs, [blockname], {blockname: code}, kGood, nickname)
		
	
class BellPairCnot(TransCnot):

	def _propagateCodes(self, codes):
		return BellPairCode(codes)
	
class BellPairCode(TransCnotCode):
	
	def stabilizerGenerators(self):
		raise NotImplementedError
		return TransCnotCode.stabilizerGenerators(self)
	


			
		
class BellPair(Component):
	
	plusName = '|+>'
	zeroName = '|0>'
	cnotName = 'cnot'
	
	def __init__(self, kGood, plus, zero, kGoodCnot):
		super(BellPair, self).__init__(kGood, subcomponents={self.plusName: plus, self.zeroName: zero})
		self.kGoodCnot = kGoodCnot							
		
	def _count(self, noise):
		prepBlocks = super(BellPair, self)._count(noise)
		
		# Construct a transversal CNOT component from the two input codes.
		cnot = BellPairCnot(self.kGoodCnot, prepBlocks[self.plusName].getCode(), prepBlocks[self.zeroName].getCode())
		prepBlocks[self.cnotName] = cnot.count(noise)
		
		return prepBlocks
	
	def _convolve(self, blocks):
		# Assume that key metadata is identical for all blocks.
		# This is a safe assumption since each block should be encoded in the
		# same code, and the parity checks for a given code are fixed.
		parityChecks = blocks[self.plusName].keyMeta().parityChecks()
		
		# Masks that identify X stabilizers and operators, and Z stabilizers
		# and operators.
		xmask = bits.listToBits((0 == check[zType]) for check in parityChecks)
		zmask = bits.listToBits((0 == check[xType]) for check in parityChecks)
		
		blockShift = len(parityChecks)

		for pauli in self.kGood.keys(): 	
			# Propagate the preparation errors through the CNOT.
			# TODO: this assumes a particular ordering of the blocks in the CNOT error keys.
			# Should use CNOT key metadata to make sure the bit shifting is correct.
			blocks[self.plusName].counts()[pauli] = mapKeys(blocks[self.plusName].counts()[pauli], 
														    lambda e: (e << blockShift) + (e & xmask))
			blocks[self.zeroName].counts()[pauli] = mapKeys(blocks[self.zeroName].counts()[pauli], 
														    lambda e: ((e & zmask) << blockShift) + e)
			
		# Now convolve.
		convolved = super(BellPair, self)._convolve(blocks)
		
		# The default _convolve() may not assign the correct block properties (i.e., key generators, code).
		cnotBlock = blocks[self.cnotName]
		return CountedBlock(convolved.counts(), cnotBlock.keyMeta(), code=cnotBlock.getCode(), name=str(convolved))
	
	
class BellMeas(Component):

	cnotName = 'cnot'
	measXName = 'measX'
	measZName = 'measZ'
	
	def __init__(self, kGood, code, kGoodMeasX=None, kGoodMeasZ=None, kGoodCnot=None):
		if None == kGoodMeasX: kGoodMeasX = kGood
		if None == kGoodMeasZ: kGoodMeasZ = kGood
		if None == kGoodCnot: kGoodCnot = kGood
		
		subs = {self.cnotName: TransCnot(kGoodCnot, code, code),
			    self.measXName: TransMeas(kGoodMeasX, code, Pauli.X),
			    self.measZName: TransMeas(kGoodMeasZ, code, Pauli.Z)}
		
		super(BellMeas, self).__init__(kGood, subcomponents=subs)
		
	def _convolve(self, blocks):
		blockShift = blocks[self.measXName].keyMeta().length()
		
		for pauli in self.kGood.keys():
			# Extend single block X-basis measurement keys to both blocks.
			# Z-basis measurements do not need to be modified because they are represented by the LSBs.
			# TODO: this assumes a particular ordering of the blocks in the CNOT error keys.
			# Should use CNOT key metadata to make sure the bit shifting is correct.
			blocks[self.measXName].counts()[pauli] = mapKeys(blocks[self.measXName].counts()[pauli], 
															 lambda e: e << blockShift)
			
		# Now convolve.
		convolved = super(BellMeas, self)._convolve(blocks)
		
		# The default _convolve() may not assign the correct code.
		cnotBlock = blocks[self.cnotName]
		return CountedBlock(convolved.counts(), cnotBlock.keyMeta(), code=cnotBlock.getCode(), name=str(convolved))
			
		
	
class Teleport(Component):
	
	bpName = 'BP'
	bmName = 'BM'
	
	def __init__(self, kGood, data, bellPair, bellMeas):
		subs = {self.bpName: bellPair,
				self.bmName: bellMeas}
		
		super(Teleport, self).__init__(kGood, subcomponents=subs)
		self._data = data
		
	def _convolve(self, blocks):
		parityChecks = self._data.keyMeta().parityChecks()
		
		# Masks that identify X stabilizers and operators, and Z stabilizers
		# and operators.
		xmask = bits.listToBits((0 == check[zType]) for check in parityChecks)
		
		blockShift = len(parityChecks)

		for pauli in self.kGood.keys(): 	
			# Propagate the input errors through the CNOT.
			# TODO: this assumes a particular ordering of the blocks in the bell-pair
			# and bell-measurement error keys.
			# Should use CNOT key metadata to make sure the bit shifting is correct.
			self._data.counts()[pauli] = mapKeys(self._data.counts()[pauli], 
												 lambda e: (e << blockShift) + (e & xmask))
			
		blocks['input'] = self._data
		convolved = super(Teleport, self)._convolve(blocks)
		
		return CountedBlock(convolved.counts(), convolved.keyMeta(), code=self._data.getCode(), name=str(convolved))
	
class TeleportED(Teleport):
	
	def _postCount(self, block):
		'''
		Post-select for the trivial syndrome on the two measured blocks.
		'''
		counts = block.counts()
		keyMeta = block.keyMeta()
		stabilizers = set(block.getCode().stabilizerGenerators())
		blockChecks = keyMeta.parityChecks()
		parityChecks = []
		shifts = [keyMeta.blockRange(name)[0] for name in keyMeta.blocknames()]
		for shift in shifts:
			checks = [(check << shift) for check in blockChecks]
			parityChecks = checks + parityChecks
		
		syndromeBits = [i for i,check in enumerate(parityChecks) if check in stabilizers]
		rejectMask = bits.listToBits(syndromeBits)
		
		for pauli in counts:
			for k,count in enumerate(counts[pauli]):
				accepted = {key: c for key,c in count.iteritems() if 0 == (c & rejectMask)}
				counts[pauli][k] = accepted
			
		# TODO: rejected counts, prAccept, code, etc.
		return block
		
		
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
#						     kGood=self.kGood)
#		
#		countsVerified = convolve(countsA[Pauli.X], 
#								  countsBCM, 
#								  convolveFcn=convolveCountsPostselectX, 
#								  kGood=self.kGood)
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
#		countsVerify = convolve(countsPrepB, countsC, kGood=kGood)
#		countsVerified = convolve(countsPrepA, countsVerify, kGood=kGood)
#	
#		
#		
#		return CountedBlock(str(countsA), countsA.getCode, results)
		
		
	