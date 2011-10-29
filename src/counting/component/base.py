'''
Created on 2011-10-25

@author: adam
'''
from copy import copy
from counting import probability, key, countErrors
from counting.block import Block
from counting.countErrors import mapCounts, countBlocksBySyndrome
from counting.countParallel import convolve
from counting.key import SyndromeKeyGenerator, SyndromeKeyMeta,\
    KeySplitter, KeyManipulator, KeyConcatenator, KeyExtender
from counting.location import Locations
from counting.result import CountResult
from qec.error import Pauli
from util import counterUtils
from util.cache import fetchable, memoize
import logging
import operator

logger = logging.getLogger('component')


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
        
    def __setitem__(self, name, component):
        self._subs[name] = component
    
    def __getitem__(self, name):
        return self._subs[name]
    
    def locations(self, pauli=Pauli.Y):
        '''
        Returns the set of locations contained in the component (and all sub-components).
        :rtype: :class:`Locations`
        '''
        locs = Locations()
        for sub in self._subs.values():
            locs += sub.locations(pauli)
        return locs + self.internalLocations(pauli)
    
    def internalLocations(self, pauli=Pauli.Y):
        '''
        Returns the set of locations contained in the component *excluding* any sub-components.
        :rtype: :class:`Locations`
        '''
        return Locations()
    
    def inBlocks(self):
        return tuple()
    
    def outBlocks(self):
        raise NotImplementedError
    
    def logicalStabilizers(self):
        return tuple([])
    
    def subcomponents(self):
        return self._subs
        
    @fetchable
    def count(self, noiseModels, pauli):
        '''
        Counts errors in the component.
        Returns a CountResult. 
        '''
        logger.info('Counting ' + str(self) + ': ' + str(pauli))
        results = self._count(noiseModels, pauli)
        
        logger.debug('Convolving ' + str(self))
        result = self._convolve(results, noiseModels, pauli)
        
        logger.debug('Post processing ' + str(self))
        return self._postCount(result, noiseModels, pauli)
    
    def prBad(self, noise, pauli, kMax=None):
        prSelf = probability.prBadPoly(self.kGood[pauli], self.internalLocations(pauli), noise, kMax)
        prSubs = [sub.prBad(noise, pauli, kMax=self.kGood[pauli]) for sub in self.subcomponents().values()] + [prSelf]
        return sum(prSubs) + prSelf
    
    def prAccept(self, noiseModels, kMax=None):
        prSubs = [sub.prAccept(noiseModels, kMax) for sub in self.subcomponents().values()]
        return reduce(operator.mul, prSubs, 1)
    
    def propagateCounts(self, counts, keyMeta):
        propagator = self.keyPropagator(keyMeta)
        propagated = mapCounts(counts, propagator)        
        return propagated, propagator.meta()
    
    def keyPropagator(self, keyMeta):
        return keyMeta
    
    def _count(self, *args):
        '''
        Subclass hook.
        Counts the errors in each sub-component and returns the counts
        as a dictionary indexed by sub-component name.  
        It is expected that most concrete components will not need to 
        implement this method. 
        '''
        return {name: sub.count(*args) for name, sub in self._subs.iteritems()} 
    
    def _convolve(self, results, noiseModels, pauli):
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
        if not all(len(r.blocks) == len(blocks) for r in results):
            raise Exception('Block count mismatch. {0}'.format([len(r.blocks) for r in results]))
        
        counts = [result.counts for result in results]
        convolved = counts[0]
        k = self.kGood[pauli]
        for count in counts[1:]:
            convolved = convolve(convolved, count, kMax=k, convolveFcn=key.convolveKeyCounts, extraArgs=[keyMeta])
            
        return CountResult(convolved, keyMeta, blocks)
    
    def _postCount(self, result, noiseModels, pauli):
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
        
        self._locations = locations
        self.blocks = tuple([Block(name, codes[name]) for name in blocknames])
        
    def internalLocations(self, pauli=Pauli.Y):
        return countErrors.pauliFilter(copy(self._locations), pauli)
    
    def inBlocks(self):
        return self.blocks
    
    def outBlocks(self):
        return self.inBlocks()
        
    def _count(self, noiseModels, pauli):
        # First, count the sub-components.
        subcounts = super(CountableComponent, self)._count(noiseModels, pauli)
            
        # Now count the internal locations.
        locations = self.internalLocations(pauli)
        counts, meta = countBlocksBySyndrome(locations, self.blocks, noiseModels[pauli], self.kGood[pauli])
        
        cb = CountResult(counts, meta, self.blocks, name=self.nickname())
        
        subcounts[self.nickname()] = cb
        return subcounts
    
class Empty(CountableComponent):
    
    def __init__(self, code, nblocks=1):
        locs = Locations([])
        blocknames = [str(i) for i in range(nblocks)]
        codes = {name: code for name in blocknames}
        super(Empty, self).__init__(locs, blocknames, codes, {})

class Prep(CountableComponent):
    
    def __init__(self, kGood, locations, code):
        blocknames = list(locations.blocknames())
        super(Prep, self).__init__(locations, blocknames, {name: code for name in blocknames}, kGood)
        
class InputDependentComponent(Component):

    def count(self, noiseModels, pauli, inputResult):
        # Counting and convolving may manipulate the input directly.
        # Copy to avoid changing the user's input.
        inputResult = copy(inputResult)
        
        logger.info('Counting ' + str(self) + ': ' + str(pauli))
        results = self._count(noiseModels, pauli)
        
        logger.debug('Convolving ' + str(self))
        result = self._convolve(results, noiseModels, pauli, inputResult)
        
        logger.debug('Post processing ' + str(self))
        return self._postCount(result, noiseModels, pauli)
        
    @memoize
    def _count(self, noiseModels, pauli):
        return super(InputDependentComponent, self)._count(noiseModels, pauli)
    
    def _convolve(self, results, noiseModels, pauli, inputResult):
        inputResult.counts, inputResult.keyMeta = self._propagateInput(inputResult)
        inputResult.blocks = results.values()[0].blocks
        results['input'] = inputResult
        return super(InputDependentComponent, self)._convolve(results, noiseModels, pauli)

    def _propagateInput(self, inputResult):
        raise NotImplementedError
    
class InputAdapter(Component):
    
    def __init__(self, component, inputKey):
        blocks = component.inBlocks()
        code = blocks[0].getCode()
        parityChecks = SyndromeKeyGenerator.ParityChecks(code)
        keyMeta = SyndromeKeyMeta(parityChecks, len(blocks))
        inputResult = CountResult([{inputKey: 1}], keyMeta, blocks)
        
        super(InputAdapter, self).__init__({})
        
        self._component = component
        self._inResult = inputResult
        
    def count(self, noiseModels, pauli):
        return self._component.count(noiseModels, pauli, self._inResult)
        
class ConcatenatedComponent(Component):
    
    def __init__(self, kGood, *components):
        subs = {i: comp for i,comp in enumerate(components)}
        super(ConcatenatedComponent, self).__init__(kGood, subcomponents=subs)
        
    def inBlocks(self):
        return sum((self[i].inBlocks() for i in range(len(self.subcomponents()))), tuple())
    
    def outBlocks(self):
        return sum((self[i].outBlocks() for i in range(len(self.subcomponents()))), tuple())
    
        
    def propagateCounts(self, counts, keyMeta):
        # Split the incoming keys according to each of the components
        nsubs = len(self.subcomponents())
        splits = [0] * (nsubs - 1)
        for i in range(nsubs-1):
            sub = self[i]
            nblocks = len(sub.inBlocks())
            splits[i] = splits[i-1] + nblocks
        
        splitter = KeySplitter(keyMeta, splits)
        # Then propagate the keys through the components
        propagator = self.Propagator(self.subcomponents(), splitter)
        
        # Then glue them back together.
        mapper = KeyConcatenator(propagator)
        
        counts = mapCounts(counts, mapper)
        return counts, mapper.meta()
        
    def _convolve(self, results, noiseModels, pauli):
        
        # Extend each of the component results so that they each contain the same
        # number of output blocks.  Then they can be convolved with the default
        # algorithm.
        blocks = sum([results[i].blocks for i in range(len(results))], tuple())
        n = len(blocks)
        blocknum = 0
        for i in range(len(results)):
            result = results[i]
            before = blocknum
            after = n - len(result.blocks) - before
            blocknum += len(result.blocks)
            
            extender = KeyExtender(result.keyMeta, blocksBefore=before, blocksAfter=after)
            result.counts = mapCounts(result.counts, extender)
            result.keyMeta = extender.meta()
            result.blocks = blocks
            
        return super(ConcatenatedComponent, self)._convolve(results, noiseModels, pauli)
    
    class Propagator(KeyManipulator):
        
        def __init__(self, subcomponents, splitter):
            super(ConcatenatedComponent.Propagator, self).__init__(splitter)
            
            self._metas = splitter.meta()
            self._propagators = [subcomponents[i].keyPropagator(self._metas[i]) 
                                 for i in range(len(subcomponents))]
            
        def meta(self):
            return self._metas
        
        def _manipulate(self, keys):
            key = tuple(self._propagators[i](key) for i,key in enumerate(keys))
            return key
    