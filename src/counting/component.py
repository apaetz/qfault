'''
Created on May 1, 2011

@author: Adam
'''
from abc import abstractmethod, ABCMeta
from copy import copy
from counting import probability, block, key, countErrors
from counting.block import Block
from counting.countErrors import countBlocksBySyndrome, extendCounts, mapCounts
from counting.countParallel import convolve
from counting.key import copyKeys, extendKeys, keyForBlock, keyCopier, \
	keySplitter, keyConcatenator, keyExtender
from counting.location import Locations
from qec.error import Pauli, xType, zType
from qec.qecc import StabilizerCode, StabilizerState
from result import CountResult
from util import counterUtils, bits, listutils
from util.cache import fetchable


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
	
	def outBlocks(self):
		raise NotImplementedError
	
	def subcomponents(self):
		return self._subs
		
	#@fetchable
	def count(self, noise, pauli):
		'''
		Counts errors in the component.
		Returns a CountResult. 
		'''
		print 'Counting', str(self)
		countedBlocks = self._count(noise, pauli)
		print 'Convolving', str(self)
		block = self._convolve(countedBlocks, pauli)
		print 'Post processing', str(self)
		return self._postCount(block)
	
	def prBad(self, noise, pauli, kMax={pauli: None for pauli in [Pauli.X, Pauli.Z, Pauli.Y]}):
		prSelf = probability.prBadPoly(self.kGood[pauli], self.internalLocations(pauli), noise, kMax[pauli])
		prSubs = [sub.prBad(noise, pauli, kMax=self.kGood) for sub in self._subs] + [prSelf]
		return sum(prSubs) + prSelf
	
	def propagateCounts(self, counts, keyMeta, blockname):
		propagator, keyMeta = self.keyPropagator(keyMeta, blockname)
		propagated = mapCounts(counts, propagator)		
		return propagated, keyMeta
	
	def keyPropagator(self, keyMeta, blockname):
		return (lambda key: key, keyMeta)
	
	def _count(self, noise, pauli):
		'''
		Subclass hook.
		Counts the errors in each sub-component and returns the counts
		as a dictionary indexed by sub-component name.  
		It is expected that most concrete components will not need to 
		implement this method. 
		'''
		return {name: sub.count(noise, pauli) for name, sub in self._subs.iteritems()} 
	
	def _convolve(self, results, pauli):
		'''
		Subclass hook.
		Combines errors from each of the sub-components.
		The default implementation simply returns 'counts'.
		Any non-trivial component will need to implement this method.
		
		Returns a CountResult object.
		'''
		
		# The sub-component names won't be used
		results = results.values()
		
		if 1 == len(results):
			return results[0]
		
		keyMeta = results[0].keyMeta
		if not all(r.keyMeta == keyMeta for r in results):
			raise Exception('Key metadatas are not all identical. {0}'.format([r.keyMeta for r in results]))
		
		blocks = results[0].blocks
		if not all(r.blocks == blocks for r in results):
			raise Exception('Result blocks are not all identical. {0}'.format([r.blocks for r in results]))
		
		counts = [result.counts for result in results]
		convolved = counts[0]
		k = self.kGood[pauli]
		for count in counts[1:]:
			convolved = convolve(convolved, count, kMax=k, convolveFcn=key.convolveKeyCounts, extraArgs=[keyMeta])
			
		return CountResult(convolved, keyMeta, blocks)
	
	def _postCount(self, result):
		'''
		Subclass hook.
		Performs (optional) error count post-processing.
		
		Returns a CountResult object.
		'''
		return result

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
		
	def internalLocations(self, pauli):
		return countErrors.pauliFilter(copy(self.locations), pauli)
	
	def outBlocks(self):
		return self.blocks
		
	def _count(self, noise, pauli):
		# First, count the sub-components.
		subcounts = super(CountableComponent, self)._count(noise, pauli)
			
		# Now count the internal locations.
		locations = self.internalLocations(pauli)
		counts, meta = countBlocksBySyndrome(locations, self.blocks, noise, self.kGood[pauli])
		
		checks = self._logicalChecks({block.name: block.getCode() for block in self.blocks})
		cb = CountResult(counts, meta, self.blocks, logicalChecks=checks, name=self.nickname())
		
		subcounts[self.nickname()] = cb
		return subcounts
	
	def _logicalChecks(self, codes):
		return []
	
class Empty(CountableComponent):
	
	def __init__(self, code, blockname='empty'):
		# We need at least one location to make the counting functions work properly.
		locs = Locations([counterUtils.locrest(blockname, 0)], blockname)
		super(Empty, self).__init__(locs, [blockname], {blockname: code}, {})

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
		
		if isinstance(ctrlCode, StabilizerState): ctrlCode = ctrlCode.getCode()
		if isinstance(targCode, StabilizerState): targCode = targCode.getCode()
		
		codes = {self.ctrlName: ctrlCode, self.targName: targCode}
		
		super(TransCnot, self).__init__(locs, blockorder, codes, kGood, nickname)
		self.blockorder = blockorder
		
	def _logicalChecks(self, codes):
		code = TransCnotCode(codes)
		return code
	
	def propagateOperator(self, op):
		ctrl = self.ctrlName
		targ = self.targName
		newOp = {ctrl: op[ctrl] ^ op[targ][zType],
				 targ: op[targ] ^ op[ctrl][xType]}
		
		return newOp
	
	def keyPropagator(self, keyMeta, blockname):
		parityChecks = keyMeta.parityChecks()
		
		blocknum = self.blockorder.index(blockname)
		
		# Extend the single input block to both blocks
		before = 1 * (1 == blocknum)
		after = 1 * (0 == blocknum)
		extender, keyMeta = keyExtender(keyMeta, blocksBefore=before, blocksAfter=after)
		
		if self.ctrlName == blockname:
			# We're on the control input.  X errors propagate through to the target
			# block.
			mask = bits.listToBits((0 == check[xType]) for check in parityChecks)
		else:
			# We're on the target input.  Z errors propagate through to the control
			# block.
			mask = bits.listToBits((0 == check[zType]) for check in parityChecks)
		
		copier = keyCopier(keyMeta, blocknum, not blocknum, mask=mask)
	
		def propagator(key):
			return copier(extender(key))
		
		return propagator, keyMeta
		
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
		
		locs = Locations([loc(blockname, i) for i in range(n)], nickname)
		super(TransMeas, self).__init__(locs, [blockname], {blockname: code}, kGood, nickname)
		
	
class BellPairCnot(TransCnot):

	def _logicalChecks(self, codes):
		return BellPairCode(codes)
	
class BellPairCode(TransCnotCode):
	
	def stabilizerGenerators(self):
		raise NotImplementedError
		return TransCnotCode.stabilizerGenerators(self)
	

class CnotConvolver(Component):
	
	ctrlName = 'ctrl'
	targName = 'targ'
	cnotName = 'cnot'
	
	def __init__(self, kGood, kGoodCnot, ctrlInput, targInput, ctrlName=ctrlName, targName=targName):
		
		# Construct a transversal CNOT component from the two input codes.
		ctrlCode = ctrlInput.outBlocks()[0].getCode()
		targCode = targInput.outBlocks()[0].getCode()
		cnot = TransCnot(kGoodCnot, ctrlCode, targCode)
		
		super(CnotConvolver, self).__init__(kGood, subcomponents={ctrlName: ctrlInput, 
																  targName: targInput,
																  self.cnotName: cnot})
		self.ctrlName = ctrlName
		self.targName = targName
		
	def outBlocks(self):
		self.subcomponents()[self.cnotName].outBlocks()
	
	def _convolve(self, results, pauli):
		
		cnot = self.subcomponents()[self.cnotName]
		ctrlResult = results[self.ctrlName]
		targResult = results[self.targName]
		cnotResult = results[self.cnotName]
		
		# First, propagate the input results through the CNOT	
		ctrlResult.counts, ctrlResult.keyMeta = cnot.propagateCounts(ctrlResult.counts,
															 ctrlResult.keyMeta,
															 cnot.ctrlName)
		targResult.counts, targResult.keyMeta = cnot.propagateCounts(targResult.counts,
															 targResult.keyMeta,
															 cnot.targName)
		
		ctrlResult.blocks = cnotResult.blocks
		targResult.blocks = cnotResult.blocks
					
		# Now convolve.
		return super(CnotConvolver, self)._convolve(results, pauli)			
		
class BellPair(CnotConvolver):
	
	ctrlName = '|+>'
	targName = '|0>'

	def __init__(self, kGood, plus, zero, kGoodCnot):
		super(BellPair, self).__init__(kGood, kGoodCnot, plus, zero, self.ctrlName, self.targName)
	
	
class BellMeas(Component):

	cnotName = 'cnot'
	measXName = 'measX'
	measZName = 'measZ'
	
	def __init__(self, kGood, code, kGoodMeasX=None, kGoodMeasZ=None, kGoodCnot=None):
		if None == kGoodMeasX: kGoodMeasX = kGood
		if None == kGoodMeasZ: kGoodMeasZ = kGood
		if None == kGoodCnot: kGoodCnot = kGood
		
		subs = {self.cnotName: TransCnot(kGoodCnot, code, code),
			    self.measXName: TransMeas(kGoodMeasX, code, Pauli.X, self.measXName),
			    self.measZName: TransMeas(kGoodMeasZ, code, Pauli.Z, self.measZName)}
		
		super(BellMeas, self).__init__(kGood, subcomponents=subs)
		self.code = code
		
	def outBlocks(self):
		measX = Block(self.measXName, self.code)
		measZ = Block(self.measZName, self.code)
		return (measX, measZ)

	def keyPropagator(self, keyMeta, blockname):
		cnot = self.subcomponents()[self.cnotName]
		namemap = {self.measXName: cnot.ctrlName, self.measZName: cnot.targName}
		return cnot.keyPropagator(keyMeta, namemap[blockname])
		
	def _convolve(self, results, pauli):
		
		cnot = results[self.cnotName]
		measX = results[self.measXName]
		measZ = results[self.measZName]
		
		measX.counts, measX.keyMeta = extendCounts(measX.counts, measX.keyMeta, blocksAfter=1)
		measZ.counts, measZ.keyMeta = extendCounts(measZ.counts, measZ.keyMeta, blocksBefore=1)
			
		measX.blocks = cnot.blocks
		measZ.blocks = cnot.blocks
			
		# Now convolve.
		return super(BellMeas, self)._convolve(results, pauli)
			
		
	
class Teleport(Component):
	
	bpName = 'BP'
	bmName = 'BM'
	
	def __init__(self, kGood, data, bellPair, bellMeas, bpMeasBlock=0):
		subs = {self.bpName: bellPair,
				self.bmName: bellMeas}
		
		super(Teleport, self).__init__(kGood, subcomponents=subs)
		self._data = data
		self._bpMeasBlock = bpMeasBlock
		
	def outBlocks(self):
		return self._data.blocks
		
	@staticmethod
	def _convolveBP(bellMeas, bellPairResult, bpMeasBlock):
		# TODO: this implementation is wrong.
		# The structure is correct, but need to compute which blocks
		# go where.
		
		bpMeta = bellPairResult.keyMeta
		splitter, meta0, meta1 = keySplitter(bpMeta, 1)
		bpSplitMeta = (meta0, meta1)
		bmPropagator, bmMeta = bellMeas.keyPropagator(bpSplitMeta[bpMeasBlock], bellMeas.measZName)
		concatenator, propMeta = keyConcatenator(bmMeta, bpSplitMeta[not bpMeasBlock])
		
		def propagator(key):
			bpKeys = splitter(key)
			key = concatenator(bmPropagator(bpKeys[bpMeasBlock]), bpKeys[not bpMeasBlock])
			return key

		bellPairResult.counts = mapCounts(bellPairResult.counts, propagator)
		bellPairResult.keyMeta = propMeta
		
		return bellPairResult

		
		
	def _convolve(self, results, pauli):
		results[self.bpName] = self._convolveBP(self.subcomponents()[self.bmName],
											    results[self.bpName], 
											    self._bpMeasBlock)
		
		# Propagate the input through the CNOT of the Bell measurement.
		bellMeas = self.subcomponents()[self.bmName]
		self._data.counts, self._data.keyMeta = bellMeas.propagateCounts(self._data.counts, 
																		 self._data.keyMeta,
																		 bellMeas.measXName)
		#self._data.blocks = results[self.bmName].blocks			
		results['input'] = self._data
		
		# Extend each of the results over all three blocks.
		bpRes = results[self.bpName]
		bpBlocksAfter = self._bpMeasBlock

		bmRes = results[self.bmName]
		bmRes.counts, bmRes.keyMeta = extendCounts(bmRes.counts, bmRes.keyMeta, not bpBlocksAfter, bpBlocksAfter)

		inRes = results['input']
		inRes.counts, inRes.keyMeta = extendCounts(inRes.counts, inRes.keyMeta, not bpBlocksAfter, bpBlocksAfter)
		
		if bpBlocksAfter:
			blocks = bpRes.blocks + inRes.blocks
		else:
			blocks = inRes.blocks + bpRes.blocks
			
		bpRes.blocks = bmRes.blocks = inRes.blocks = blocks

		return super(Teleport, self)._convolve(results, pauli)

	
class TeleportED(Teleport):
	
	def _postCount(self, result):
		'''
		Post-select for the trivial syndrome on the two measured blocks.
		'''
		counts = result.counts
		keyMeta = result.keyMeta
		stabilizers = set(result.blocks[0].getCode().stabilizerGenerators())
		blockChecks = keyMeta.parityChecks()
		parityChecks = blockChecks
		
		syndromeBits = [i for i,check in enumerate(parityChecks) if check in stabilizers]
		rejectMask = bits.listToBits(syndromeBits)
		
		for k,count in enumerate(counts):
			# TODO compute the actual block numbers to check.
			# TODO also check logical operators that are in the stabilizer.
			accepted = {keyForBlock(key, keyMeta, 2): c for key,c in count.iteritems() if 0 == (key[1] & rejectMask)}
			counts[k] = accepted
			
		# TODO: rejected counts, prAccept, code, etc.
		return result
		
		
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
#		return CountResult(str(countsA), countsA.getCode, results)
		
		
	