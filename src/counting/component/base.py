'''
Created on 2011-10-25

This file contains base classes for fault-tolerant "components".  A component
is a part of a rectangle or extended rectangle (exRec) in which combinations of
errors can be counted.

This file contains both abstract base classes, such as Component, CountableComponent,
etc., and a few simple concrete classes such as Empty, and Prep.

@author: adam
'''
from copy import copy
from counting import probability, key, countErrors
from counting.block import Block
from counting.countErrors import mapCounts, countBlocksBySyndrome
from counting.countParallel import convolve
from counting.key import KeyManipulator, SyndromeKeyGenerator, \
	IdentityManipulator, KeyMerger, KeyPermuter
from counting.location import Locations
from counting.result import CountResult
from qec.error import Pauli
from qec.qecc import ConcatenatedCode
from util import listutils
from util.cache import fetchable
import hashlib
import logging

logger = logging.getLogger('counting.component')


class Component(object):
	'''
	Abstract class for circuit "components" of the extended rectangle (exRec).
	
	To count errors, the exRec is divided up into a hierarchy of components.
	Each component contains a number of sub-components.  Errors are counted
	recursively by first counting small numbers of errors in each sub-component
	and then combining the error information.
	
	Most components can be treated as a black boxes which take as input a maximum
	number of faults k, and outputs a set of counts, k counts for each possible
	error that can occur inside of the component.  In general, the output of 
	a component may also depend on errors on its input blocks. Error counts
	are obtained by calling the count() method.
	 	
	Concrete subclasses should, at a minimum, define the input blocks of the component
	(the inBlocks() method) and probably an input key propagator (the keyPropagator())
	method.  The key propagator transforms errors on the input blocks into the corresponding
	errors on the output blocks.
	
	In many cases, the input error counts may contain more blocks than are specified by the
	component.  Consider, for example, a circuit with three components.  Component 1 outputs
	3 blocks, component 2 takes two input blocks, and component three takes 3 input blocks.
	In this case, the input to component 2 is three blocks, even though it only accepts two
	input blocks.  However, the third block must be maintained for input eventual input to
	component 3.  As a result, all components must be able to accept inputs with more blocks
	than specified by inBlocks().  Extra blocks will be preserved (as if the identity was applied
	on that subspace).
	'''
	
	def __init__(self, kGood, subcomponents=[]):
		'''
		:param dict kGood: Specifies the maximum number of faults to count, for each Pauli type.
		:param list subcomponents: (optional) A list of sub-components. 
		'''
		self.kGood = {pauli: kGood.get(pauli, 0) for pauli in (Pauli.X, Pauli.Z, Pauli.Y)}
		self._subs = subcomponents
		
		self._id = hashlib.md5(self._hashStr())
	
		
#################
# Public methods
#################

	#@fetchable
	def count(self, noiseModels, pauli, inputResult=None, kMax=None):
		'''
		Counts errors in the component.
		Returns a CountResult. 
		
		:param dict noiseModels: A dictionary, indexed by Pauli error, of noise models.
		:param pauli: The error type to count.  Use Pauli.Y to count X and Z errors together.
		'''
		
		raise NotImplementedError
		
		if None == inputResult:
			inputs = tuple([0]*len(self.inBlocks()))
			inputCounts = [{inputs: 1}]
			inputResult = CountResult(inputCounts, self.inBlocks())
			
		k = self.kGood[pauli]
		if None != kMax:
			k = min(k, kMax)
			
		try:
			self._log(logging.INFO, 'Counting: ' + str(pauli))
			
			result = inputResult
			for sub in self.subcomponents():
				result = sub.count(noiseModels, pauli, result, k)
			
			self._log(logging.DEBUG, 'counts=%s', result.counts)		
		except:
			self._log(logging.ERROR, 'Error while counting')
			raise
		
		return result
	
	def prBad(self, noise, pauli, kMax=None):
		'''
		Returns polynomial representing an upper bound on the probability that the component is
		bad.
		:param noise: The noise model.
		:param pauli: The error type
		:param kmax: The maximum number of faults to consider.
		'''
		prSelf = probability.prBadPoly(self.kGood[pauli], self.locations(pauli), noise, kMax)
		self._log(logging.DEBUG, 'Pr[bad] (self)=%s', prSelf)
		prSubs = [sub.prBad(noise, pauli, kMax=self.kGood[pauli]) for sub in self.subcomponents()]
		return sum(prSubs) + prSelf
	
#	def prAccept(self, noiseModels, kMax=None):
#		'''
#		Returns a polynomial representing a lower bound on the probability that the component
#		accepts.  (This is just 1 for components without postselection.)
#		
#		:param dict noiseModels: A dictionary of noise models, indexed by Pauli error type.
#		:param kMax: (optional) The maximum number of faults to consider. (i.e., Pr[accept, K <= kMax])
#		'''
#		prSubs = [sub.prAccept(noiseModels, kMax=self.kGood) for sub in self.subcomponents().values()]
#		return reduce(operator.mul, prSubs, SymPolyWrapper(sympoly1d([1])))
		
	def locations(self, pauli=Pauli.Y):
		'''
		Returns the set of locations contained in the component (and all sub-components).
		Note that the order of returned set of locations need not match the physical order
		of locations in the circuit.
		
		:rtype: :class:`Locations`
		'''
		locs = Locations()
		for sub in self.subcomponents():
			locs += sub.locations(pauli)
		return locs

	def inBlocks(self):
		'''
		Returns the list of input blocks.
		'''
		return tuple()
	
	def outBlocks(self):
		'''
		Returns the list of output blocks.
		'''
		return self.inBlocks()
	
#	def logicalStabilizers(self):
#		'''
#		TODO: necessary?
#		'''
#		return tuple([])
	
	def subcomponents(self):
		'''
		Returns the list of sub-components.
		'''
		return self._subs
	
	def propagateCounts(self, inputResult):
		'''
		Propagates the given counts through the component.
		
		:param dict counts: Counts indexed by [k][key][count]
		:param keyMeta: The key metadata for counts.
		:rtype tuple: The propagated counts and corresponding key metadata
		'''
		propagator = self.keyPropagator()
		# TODO: could/should verify that the inputResult blocks match the
		# expected input blocks.
		propagated = mapCounts(inputResult.counts, propagator)
		
		# Blocks in the input space of the component are converted to the output space.
		# Other input blocks are unchanged.
		blocks = self.outBlocks() + tuple(inputResult.blocks[len(self.inBlocks()):])
				
		return CountResult(propagated, blocks)
	
	def keyPropagator(self, subPropagator=IdentityManipulator()):
		'''
		Returns an object that may be used to propagate keys through the component.
		:param keyMeta: The key metadata (or another KeyPropagator)
		:rtype: :class:`KeyPropagator`
		'''
		return subPropagator
	
	def identifier(self):
		return self._id
	
###################################
# Private methods (subclass hooks)
###################################
	
#	def _count(self, *args):
#		'''
#		Subclass hook.
#		Counts the errors in each sub-component and returns the counts
#		as a dictionary indexed by sub-component name.  
#		It is expected that most concrete components will not need to 
#		implement this method. 
#		'''
#		return {name: sub.count(*args) for name, sub in self._subs.iteritems()}
#	
#	def _convolve(self, results, noiseModels, pauli):
#		'''
#		Subclass hook.
#		Combines errors from each of the sub-components.
#		The default implementation assumes that all subcomponents contain the same
#		number of blocks and that those blocks line up correctly.
#		
#		:rtype: :class:`CountResult`
#		'''
#		
#		# The sub-component names won't be used
#		results = results.values()
#		
#		if 1 == len(results):
#			return results[0]
#		
#		keyMeta = results[0].keyMeta
#		if not all((r.keyMeta == keyMeta) for r in results):
#			raise Exception('Key metadatas are not all identical. {0}'.format([r.keyMeta for r in results]))
#		
#		blocks = results[0].blocks
#		if not all(len(r.blocks) == len(blocks) for r in results):
#			raise Exception('Block count mismatch. {0}'.format([len(r.blocks) for r in results]))
#		
#		counts = [result.counts for result in results]
#		convolved = counts[0]
#		k = self.kGood[pauli]
#		for count in counts[1:]:
#			convolved = convolve(convolved, count, kMax=k, convolveFcn=key.convolveKeyCounts, extraArgs=[keyMeta])
#			
#		return CountResult(convolved, blocks)
#	
#	def _postCount(self, result, noiseModels, pauli):
#		'''
#		Subclass hook.
#		Performs (optional) error count post-processing.
#		
#		:rtype: :class:`CountResult`
#		'''
#		# TODO: this could be a good place to check for, and compute XZ corrections.
#		return result
	

#################################
	
	def __setitem__(self, name, component):
		self._subs[name] = component
	
	def __getitem__(self, name):
		return self._subs[name]
		
	def _log(self, level, msg, *args, **kwargs):
		classname = self.__class__.__name__
		logger.log(level, ''.join([classname, ': ', msg]), *args, **kwargs)
	
#	def nickname(self):
#		return self._nickname
#	
#	def fullname(self):
#		full = str(self.__class__.name)
#		if None != self._nickname:
#			full += '.' + self._nickname
#		
#		return full
	
	def descriptor(self):
		rep = str(self.__class__.__name__)
		rep += str(self.kGood)
		return rep
	
	def _hashStr(self):
		hashStr = self.descriptor() + ''.join([sub.identifier().hexdigest() for sub in self.subcomponents()])
		return hashStr
	
	def __str__(self):
		return self.__repr__()
	
	def __repr__(self):
		return ''.join([self.descriptor(), '.', self.identifier().hexdigest()])
	
class CountableComponent(Component):
	'''
	A component which, in instead of having sub-components,
	has its own physical locations that must be counted.
	'''
	
	def __init__(self, kGood, locations):
		# The number of faulty locations cannot exceed the total
		# number of locations.
				
		self._locations = locations
		
		super(CountableComponent, self).__init__(kGood)

	def count(self, noiseModels, pauli, inputResult=None, kMax=None):		
		result = self._count(noiseModels, pauli)
		
		if None == inputResult:
			return result
		
		# Avoid altering caller's data.
		inputResult = copy(inputResult)
		inputResult = self.propagateCounts(inputResult)
		# TODO: more robust way of getting key lengths?
		keyLengths = [len(SyndromeKeyGenerator(block.getCode(), None).parityChecks()) for block in result.blocks]
		inKeyLengths = [len(SyndromeKeyGenerator(block.getCode(), None).parityChecks()) for block in inputResult.blocks]
		result.counts = convolve(inputResult.counts, 
								result.counts, 
								kMax=kMax, 
								convolveFcn=key.convolveKeyCounts, 
								extraArgs=[inKeyLengths, keyLengths])
		
		result.blocks = inputResult.blocks
		
		return result
	
	@fetchable
	def _count(self, noiseModels, pauli):
		# Count the internal locations.
		locations = self.locations(pauli)
		blocks = self.outBlocks()
		counts = countBlocksBySyndrome(locations, blocks, noiseModels[pauli], self.kGood[pauli])
		
		cb = CountResult(counts, blocks, name=str(self))
		return cb
				
	def locations(self, pauli=Pauli.Y):
		return countErrors.pauliFilter(copy(self._locations), pauli)

	def _hashStr(self):
		return super(CountableComponent, self)._hashStr() + str(self._locations.list)
	
	
class CompositeComponent(Component):
	
	def count(self, noiseModels, pauli, inputResult=None, kMax=None):
		'''
		Counts errors in the component.
		Returns a CountResult. 
		
		:param dict noiseModels: A dictionary, indexed by Pauli error, of noise models.
		:param pauli: The error type to count.  Use Pauli.Y to count X and Z errors together.
		'''
		
		if None == inputResult:
			inputs = tuple([0]*len(self.inBlocks()))
			inputCounts = [{inputs: 1}]
			inputResult = CountResult(inputCounts, self.inBlocks())
			
		kIn = len(inputResult.counts) - 1
		
		kGood = self.kGood[pauli]
		if None == kMax:
			kMax = kGood + kIn
		else:
			kGood = min(kGood, kMax)
			
		try:
			self._log(logging.INFO, 'Counting: ' + str(pauli))
			
			# Convolve the input with sub-component counts by counting the
			# sub-components with a single order of the input counts at a time.
			# That is, take order-k input counts and treat them as order zero.
			# Then convolve with the sub-component counts.  Once convolved,
			# shift the result by k so that the orders of the result are
			# correct.  Counting in this way allows the sub-components to be
			# counted up to the correct fault order (no overcounting).
			
			# TODO: optimize
			results = []
			for k in range(min(kMax, kIn) + 1):
				counts = [inputResult.counts[k]]
				result = CountResult(counts, inputResult.blocks)
				
				result = self._countInputOrderZero(noiseModels, pauli, result, max(kMax-k, 0))
					
				result.counts = [{} for _ in range(k)] + result.counts + [{} for _ in range(kMax+1 - len(result.counts) - k)]
				results.append(result)
	
			resultCounts = [listutils.addDicts(*[r.counts[k] for r in results]) for k in range(kMax+1)]
			result = CountResult(resultCounts, results[0].blocks)
			
			self._log(logging.DEBUG, 'counts=%s', result.counts)		
		except:
			self._log(logging.ERROR, 'Error while counting')
			raise

		return result
	
	def _countInputOrderZero(self, noiseModels, pauli, inputResult, kMax):
		'''
		Count the sub-components with an input that has only order-zero counts.
		'''
		raise NotImplementedError
	
class SequentialComponent(CompositeComponent):
	'''
	A component for which sub-components are ordered sequentially in time.
	'''
	
	def __init__(self, kGood, subcomponents=[]):
		super(SequentialComponent, self).__init__(kGood, subcomponents=subcomponents)
		self._ValidateSubcomponents(subcomponents)
	
	@staticmethod
	def _ValidateSubcomponents(subcomponents):
		'''
		Validates the ordered list of subcomponents by checking the output blocks of
		component i against the input blocks of component (i+1), for all components in
		the list.
		'''
		if len(subcomponents) <= 1:
			return
		
		outSub = subcomponents[0]
		for inSub in subcomponents[1:]:
			outBlocks = outSub.outBlocks()
			inBlocks = inSub.inBlocks()
			if len(outBlocks) != len(inBlocks) or\
			   not all(outBlocks[i] == inBlocks[i] for i in range(len(outBlocks))):
				raise Exception('Subcomponent input/output mismatch. {0}: {1}, {2}: {3}'.format(outSub, outBlocks, inSub, inBlocks))
			
			outSub = inSub
			
	def inBlocks(self):
		return self[0].inBlocks()
	
	def outBlocks(self):
		return self[-1].outBlocks()
	
	def _countInputOrderZero(self, noiseModels, pauli, inputResult, kMax):
		kMaxSub = min(self.kGood[pauli], kMax)
		result = inputResult
		for sub in self.subcomponents():
			result = sub.count(noiseModels, pauli, result, kMaxSub)
		return result
			
#	def count(self, noiseModels, pauli, inputResult=None, kMax=None):
#		'''
#		Counts errors in the component.
#		Returns a CountResult. 
#		
#		:param dict noiseModels: A dictionary, indexed by Pauli error, of noise models.
#		:param pauli: The error type to count.  Use Pauli.Y to count X and Z errors together.
#		'''
#		
#		if None == inputResult:
#			inputs = tuple([0]*len(self.inBlocks()))
#			inputCounts = [{inputs: 1}]
#			inputResult = CountResult(inputCounts, self.inBlocks())
#			
#		kIn = len(inputResult.counts) - 1
#		
#		kGood = self.kGood[pauli]
#		if None == kMax:
#			kMax = kGood + kIn
#		else:
#			kGood = min(kGood, kMax)
#			
#		try:
#			self._log(logging.INFO, 'Counting: ' + str(pauli))
#			
#			# TODO: generalize/abstract/optimize this process.
#			results = []
#			for k in range(min(kMax, kIn) + 1):
#				counts = [inputResult.counts[k]]
#				result = CountResult(counts, inputResult.blocks)
#				
#				for sub in self.subcomponents():
#					kMaxSub = max(min(kGood, kMax-k), 0)
#					result = sub.count(noiseModels, pauli, result, kMaxSub)
#					
#				result.counts = [{} for _ in range(k)] + result.counts + [{} for _ in range(kMax+1 - len(result.counts) - k)]
#				results.append(result)
#	
#			try:
#				resultCounts = [listutils.addDicts(*[r.counts[k] for r in results]) for k in range(kMax+1)]
#				result = CountResult(resultCounts, results[0].blocks)
#			except Exception:
#				raise
#			
#			self._log(logging.DEBUG, 'counts=%s', result.counts)		
#		except:
#			self._log(logging.ERROR, 'Error while counting')
#			raise
#		
#		if 0 == len(result.counts):
#			print 'foo?'
#		return result
	
	def keyPropagator(self, subPropagator=IdentityManipulator()):
		propagator = subPropagator
		
		# Propagate sequentially through each sub-component.
		for sub in self:
			propagator = sub.keyPropagator(propagator)
			
		return propagator
					
	
class Empty(CountableComponent):
	'''
	A completely empty component.
	'''
	
	def __init__(self, code, blockname='0'):
		locs = Locations([])
		kGood = {}
		super(Empty, self).__init__(kGood, locs)
		self._block = Block(blockname, code)
		
	def inBlocks(self):
		return (self._block,)
	
class Prep(CountableComponent):
	'''
	Codeword preparation.
	'''
	
	def __init__(self, kGood, locations, code):
		blocknames = list(locations.blocknames())
		self._outBlocks = tuple(Block(name, code) for name in blocknames)
		super(Prep, self).__init__(kGood, locations)
		
	def inBlocks(self):
		# TODO: Strictly speaking, the input is a set of 1-qubit blocks (i.e. an unencoded block).
		# But returning the output blocks here allows count propagation to work correctly.
		return self.outBlocks()
		
	def outBlocks(self):
		return self._outBlocks
	
class Filter(Component):
	'''
	A component which acts on its input only.
	'''
	
	def __init__(self):
		super(Filter, self).__init__({})
		
	def count(self, noiseModels, pauli, inputResult=None, kMax=None):
		return self.propagateCounts(inputResult)

	def keyPropagator(self, subPropagator=IdentityManipulator()):
		'''
		Subclasses must implement a post count method.
		'''
		raise NotImplementedError	
	
class PostselectionFilter(Filter):
	'''
	Special case of a Filter component in which postselection is used.  That is, some inputs
	are 'accepted' and other inputs are 'rejected', i.e., removed.
	'''
	
	def __init__(self, pauliDependency=Pauli.Y):
		super(PostselectionFilter, self).__init__()
		self._pauliDependency = pauliDependency

	def count(self, noiseModels, pauli, inputResult=None, kMax=None):
		result = super(PostselectionFilter, self).count(noiseModels, pauli, inputResult, kMax)
		
		# Remove the rejected counts
		for count in result.counts:
			count.pop(None, None)
			
		return result

	def prAccept(self, noiseModels, inputResult, kMax):
		'''
		Computes a lower bound on the acceptance probability using upper bounds on the
		rejection probability.
		'''
		
		# First, compute the acceptance probabililty of the subcomponents.
		#prSubs = super(PostselectionFilter, self).prAccept(noiseModels, inputResult, kMax)
	
		# TODO: does it make sense to have a Pauli dependency?	
		pauli = self._pauliDependency
	
		propagated = self.propagateCounts(inputResult)
		
		rejected = [{None: count.get(None, 0)} for count in propagated.counts]
		
		prBad = self.prBad(noiseModels[pauli], pauli, kMax)
		locTotals = self.locations(pauli).getTotals()
		
		prAccept = 1 - probability.upperBoundPoly(rejected, prBad, locTotals, noiseModels[pauli])
		self._log(logging.DEBUG, 'Pr[accept] = %s', prAccept)
		
		return prAccept# * prSubs
	
	
class ConcatenationFilter(Filter):
	
	def __init__(self, topCode, bottomCode):
		super(ConcatenationFilter, self).__init__()
		self.top = topCode
		self.bottom = bottomCode
	
	def inBlocks(self):
		return tuple([Block('Subblock', self.bottom)]*self.top.n)
	
	def outBlocks(self):
		return (Block('Cat', self.top), )

#	def outBlocks(self):
#		'''
#		Output a single level-2 block rather than four level-1 blocks.
#		'''
#		subblocks = super(BellPairLevel2, self).outBlocks()
#		code = subblocks[0].getCode()
#		catCode = ConcatenatedCode(code, code)
#		return tuple([Block('2-Prep', catCode)])


	def keyPropagator(self, subPropagator=IdentityManipulator()):
		'''
		Convert the four level-1 blocks into a single level-2 block.
		'''
		subblockLength = len(SyndromeKeyGenerator(self.bottom, None).parityChecks())
		keyLengths = [subblockLength] * self.top.n

		return KeyMerger(subPropagator, keyLengths)
	
class ParallelComponent(CompositeComponent):
	'''
	A component for which sub-components are parallel in time (and sequential in space).
	
	Wraps multiple independent components so that they act as a single component.
	This is useful, for example, when a component with many output blocks
	is connected to several components with smaller numbers of inputs.
	'''
	
	def __init__(self, kGood, *components):
		subs = tuple(components)
		super(ParallelComponent, self).__init__(kGood, subcomponents=subs)
		
	def inBlocks(self):
		return sum((self[i].inBlocks() for i in range(len(self.subcomponents()))), tuple())
	
	def outBlocks(self):
		return sum((self[i].outBlocks() for i in range(len(self.subcomponents()))), tuple())

	def keyPropagator(self, subPropagator=IdentityManipulator()):
		propagator = subPropagator
		for sub in self:
			propagator = sub.keyPropagator(propagator)
			
			# Shift the input blocks for the next component into place.
			propagator = self.TupleRotator(len(sub.outBlocks()), propagator)
			
		# It is possible that the inputResult space is larger than the output space
		# of the parallel component.
		propagator = self.TupleRotator(-len(self.outBlocks()), propagator)

		return propagator
	
	def _countInputOrderZero(self, noiseModels, pauli, inputResult, kMax):
		kMaxSub = min(self.kGood[pauli], kMax)
		result = inputResult
		for sub in self:
			result = sub.count(noiseModels, pauli, result, kMaxSub)
			
			# Shift the input blocks for the next component into place.
			rotator = self.TupleRotator(len(sub.outBlocks()))
			result.counts = mapCounts(result.counts, rotator)
			result.blocks = rotator(result.blocks)
			
		# It is possible that the inputResult space is larger than the output space
		# of the parallel component.
		rotator = self.TupleRotator(-len(self.outBlocks()))
		result.counts = mapCounts(result.counts, rotator)
		result.blocks = rotator(result.blocks)
		
		return result
		
	def count(self, noiseModels, pauli, inputResult=None, kMax=None):

		# The idea here is to count each of the sub-components sequentially, but
		# permuting the input blocks at each step.
		
		
		k = self.kGood[pauli]
		if None != kMax:
			k = min(k, kMax)
		
		result = copy(inputResult)
		for sub in self:
			result = sub.count(noiseModels, pauli, result, k)
			if 0 == len(result.counts):
				print 'foo?'
			
			# Shift the input blocks for the next component into place.
			rotator = self.TupleRotator(len(sub.outBlocks()))
			result.counts = mapCounts(result.counts, rotator)
			result.blocks = rotator(result.blocks)
			
		# It is possible that the inputResult space is larger than the output space
		# of the parallel component.
		rotator = self.TupleRotator(-len(self.outBlocks()))
		result.counts = mapCounts(result.counts, rotator)
		result.blocks = rotator(result.blocks)
		
		if 0 == len(result.counts):
			print 'foo?'
		return result
	
#	def __repr__(self):
#		subStr = ''.join(sub.__class__.__name__ for sub in self)
#		return super(ParallelComponent, self).__repr__() + '.' + subStr
	
	class TupleRotator(KeyManipulator):
		'''
		Performs a left rotation of a tuple by a specified number of indices.
		'''
		
		def __init__(self, rotation, manipulator=IdentityManipulator()):
			super(ParallelComponent.TupleRotator, self).__init__(manipulator)
			self.rotation = rotation
			
		def _manipulate(self, tup):
			return tup[self.rotation:] + tup[:self.rotation]
		

	
	