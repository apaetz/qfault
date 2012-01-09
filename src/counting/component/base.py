'''
Created on 2011-10-25

@author: adam
'''
from copy import copy
from counting import probability, key, countErrors
from counting.block import Block
from counting.countErrors import mapCounts, countBlocksBySyndrome
from counting.countParallel import convolve
from counting.key import KeySplitter, KeyManipulator, KeyConcatenator, KeyExtender,\
	SyndromeKeyGenerator, MultiBlockSyndromeKeyGenerator, IdentityManipulator
from counting.location import Locations
from counting.result import CountResult
from qec.error import Pauli
from util.cache import fetchable, memoize
import logging
import operator
from util.polynomial import sympoly1d, SymPolyWrapper
import hashlib

logger = logging.getLogger('counting.component')


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
	
	Other useful public methods include prBad(), and prAccept(), which return
	polynomials representing the probability that the component is "bad" and
	that the component accepts, respectively.
	'''
	
	inputName = 'input'

	def __init__(self, kGood, subcomponents=[]):
		'''
		Constructor
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
		return locs + self.internalLocations(pauli)
	
	def internalLocations(self, pauli=Pauli.Y):
		'''
		Returns the set of locations contained in the component *excluding* any sub-components.
		:rtype: :class:`Locations`
		'''
		return Locations()
	
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
	
	def logicalStabilizers(self):
		'''
		TODO: necessary?
		'''
		return tuple([])
	
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
		
		if len(propagated[0].keys()[0]) != len(blocks):
			print self.outBlocks(), self.inBlocks()
			print inputResult.blocks
			print 'foo'
				
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
	A component which, in addition to (or instead of) having sub-components,
	also has its own physical locations that must be counted.
	'''
	
	locsName = 'locations'
	
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
		locations = self.internalLocations(pauli)
		blocks = self.outBlocks()
		counts = countBlocksBySyndrome(locations, blocks, noiseModels[pauli], self.kGood[pauli])
		
		cb = CountResult(counts, blocks, name=str(self))
		return cb
				
	def internalLocations(self, pauli=Pauli.Y):
		return countErrors.pauliFilter(copy(self._locations), pauli)

	def _hashStr(self):
		return super(CountableComponent, self)._hashStr() + str(self._locations.list)
	
	
class SequentialComponent(Component):
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
	
#class FixedOutput(Empty):
#	
#	def __init__(self, code, outputCounts, keyMeta):
#		self._outputCounts = outputCounts
#		self._keyMeta = keyMeta
#		super(FixedOutput, self).__init__(code)
#	
#	def count(self, noiseModels, pauli):
#		return CountResult(self._outputCounts, self.outBlocks())
#
#	def _hashStr(self):
#		return super(CountableComponent, self)._hashStr() + str([self._outputCounts, self._keyMeta])
	
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
		
#class InputDependentComponent(Component):
#	'''
#	Abstract class for components which require input.
#	
#	Many components can be counted by specifying only the number of faults. Input to those components
#	can be specified later on and convolved with the count of the original component.
#	The behavior of some components, however, depends on the input.  This class accounts for the input
#	by adding an additional input parameter to the necessary Component methods.  See also :class:`InputAdapter`.
#	'''
#
#	def count(self, noiseModels, pauli, inputResult):
#		# Counting and convolving may manipulate the input directly.
#		# Copy to avoid changing the user's input.
#		inputResult = copy(inputResult)
#		
#		logger.info('Counting : ' + str(pauli))
#		results = self._count(noiseModels, pauli)
#		
#		logger.debug('Convolving')
#		result = self._convolve(results, noiseModels, pauli, inputResult)
#		
#		logger.debug('Post processing')
#		return self._postCount(result, noiseModels, pauli)
#	
#	def prAccept(self, noiseModels, inputResult, kMax=None):
#		return super(InputDependentComponent, self).prAccept(noiseModels, kMax)
#		
#	@memoize
#	def _count(self, noiseModels, pauli):
#		return super(InputDependentComponent, self)._count(noiseModels, pauli)
#	
#	def _convolve(self, results, noiseModels, pauli, inputResult):
#		inputResult = self._prepareInput(inputResult)
#		if len(results):
#			inputResult.blocks = results.values()[0].blocks
#		results['input'] = inputResult
#		return super(InputDependentComponent, self)._convolve(results, noiseModels, pauli)
#
#	def _prepareInput(self, inputResult):
#		raise NotImplementedError
#	
#	def keyPropagator(self, keyMeta):
#		raise Exception('Cannot propagate keys through an input dependent component.')
	
#class EmptyIDC(InputDependentComponent):
#	'''
#	A completely empty component, that takes an input.
#	'''
#	
#	def __init__(self, code):
#		kGood = {}
#		self.blocks = (Block('0', code), )
#		super(EmptyIDC, self).__init__(kGood)
#		
#	def inBlocks(self):
#		return self.blocks
#	
#	def _prepareInput(self, inputResult):
#		return inputResult
	
class Filter(Component):
	
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
	Special case of an input-dependent component in which postselection is used.
	Postselecting components accept some errors, passing them to the output, and reject other errors.
	'''
	
	def __init__(self, pauliDependency=Pauli.Y):
		super(PostselectionFilter, self).__init__()
		self._pauliDependency = pauliDependency

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
		
		rejected = [count.get(None, 0) for count in propagated.counts]
		
		prBad = self.prBad(noiseModels[pauli], pauli, kMax)
		locTotals = self.locations(pauli).getTotals()
		
		prAccept = 1 - probability.upperBoundPoly(rejected, prBad, locTotals, noiseModels[pauli])
		self._log(logging.DEBUG, 'Pr[accept] = %s', prAccept)
		
		return prAccept# * prSubs
	
class ParallelComponent(Component):
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
		# Split the incoming keys according to each of the components
		nsubs = len(self.subcomponents())
		splits = [0] * (nsubs - 1)
		for i in range(nsubs-1):
			sub = self[i]
			nblocks = len(sub.inBlocks())
			splits[i] = splits[i-1] + nblocks
		
		splitter = KeySplitter(subPropagator, splits)
		# Then propagate the keys through the components
		propagator = self.Propagator(self.subcomponents(), splitter)
		
		# Then glue them back together.
		mapper = KeyConcatenator(propagator)
		return mapper
		
	def count(self, noiseModels, pauli, inputResult=None, kMax=None):

		# The idea here is to count each of the sub-components sequentially, but
		# permuting the input blocks at each step.
		
		
		k = self.kGood[pauli]
		if None != kMax:
			k = min(k, kMax)
		
		result = inputResult
		for sub in self:
			result = sub.count(noiseModels, pauli, result, k)
			
			# Shift the input blocks for the next component into place.
			rotator = self.TupleRotator(len(sub.outBlocks()))
			result.counts = mapCounts(result.counts, rotator)
			result.blocks = rotator(result.blocks)
			
		# It is possible that the inputResult space is larger than the output space
		# of the parallel component.
		rotator = self.TupleRotator(len(inputResult.blocks) - len(self.outBlocks()))
		result.counts = mapCounts(result.counts, rotator)
		result.blocks = rotator(result.blocks)
			
		if len(result.counts[0].keys()[0]) != len(result.blocks):
			print self.outBlocks(), self.inBlocks()
			print result.blocks
			print 'foo'
		
		return result
	
#	def __repr__(self):
#		subStr = ''.join(sub.__class__.__name__ for sub in self)
#		return super(ParallelComponent, self).__repr__() + '.' + subStr
	
	class TupleRotator(object):
		'''
		Performs a left rotation of a tuple by a specified number of indices.
		'''
		
		def __init__(self, rotation):
			self.rotation = rotation
			
		def __call__(self, tup):
			return tup[self.rotation:] + tup[:self.rotation]
	