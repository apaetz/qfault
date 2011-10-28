'''
Created on 2011-10-25

@author: adam
'''
from copy import copy
from counting import probability
from counting.component.base import Component, InputDependentComponent
from counting.countErrors import mapCounts, extendCounts
from counting.key import keyForBlock, KeyExtender, KeySplitter, KeyManipulator,\
    KeyConcatenator
from qec.error import Pauli
from util import bits
from util.cache import fetchable
import logging
from counting.block import Block

logger = logging.getLogger('component')

class UncorrectedTeleport(Component):
    
    inName = 'input'
    bpName = 'BP'
    bmName = 'BM'
    measXName = 'measX'
    measZName = 'measZ'
    outName = 'output'
    
    def __init__(self, kGood, bellPair, bellMeas):
        subs = {self.bpName: bellPair,
                self.bmName: bellMeas}
        
        super(UncorrectedTeleport, self).__init__(kGood, subcomponents=subs)
        #self._outblockOrder = outblockOrder
        
    def inBlocks(self):
        return (self[self.bmName].inBlocks()[0], )
    
    def outBlocks(self):
        #TODO
        raise NotImplementedError
    
    def keyPropagator(self, keyMeta):
        # Propagate the input through the CNOT of the Bell measurement and
        # then extend to the output block.
        bmExtender = KeyExtender(keyMeta, blocksAfter=1)
        bellMeas = self[self.bmName]
        bmPropagator = bellMeas.keyPropagator(bmExtender)
        outExtender = KeyExtender(bmPropagator, blocksAfter=1)
        
        return outExtender
        
#    @staticmethod
#    @fetchable
#    def _convolveBP(bellMeas, bellPairResult, bpMeasBlock):
#        bpMeta = bellPairResult.keyMeta
#        splitter, meta0, meta1 = keySplitter(bpMeta, 1)
#        bpSplitMeta = (meta0, meta1)
#        bmPropagator, bmMeta = bellMeas.keyPropagator(bpSplitMeta[bpMeasBlock], bellMeas.measZName)
#        concatenator, propMeta = keyConcatenator(bmMeta, bpSplitMeta[not bpMeasBlock])
#        
#        def propagator(key):
#            bpKeys = splitter(key)
#            key = concatenator(bmPropagator(bpKeys[bpMeasBlock]), bpKeys[not bpMeasBlock])
#            return key
#
#        bellPairResult.counts = mapCounts(bellPairResult.counts, propagator)
#        bellPairResult.keyMeta = propMeta
#        
#        return bellPairResult
        
    def _convolve(self, results, noiseModels, pauli):
#        results[self.bpName] = self._convolveBP(self.subcomponents()[self.bmName],
#                                                results[self.bpName], 
#                                                self._bpMeasBlock)

        bm = self[self.bmName]
        bpRes = results[self.bpName]
        bmRes = results[self.bmName]
        
        # Propagate the first half of the bell pair through the bell measurement.
        bpMeta = bpRes.keyMeta
        extender = KeyExtender(bpMeta, blocksBefore=1)
        splitter = KeySplitter(extender, [2])
        bmPropagator = self.BMPropagator(splitter, bm)
        concatenator = KeyConcatenator(bmPropagator)
        
#        def propagator(key):
#            bpKey0, bpKey1 = splitter(key)
#            key = concatenator(bmPropagator(extender(bpKey0)), bpKey1)
#            return key

        bpRes.counts = mapCounts(bpRes.counts, concatenator)
        bpRes.keyMeta = concatenator.meta()
        
        # Extend bell measurement results over all three blocks.
        extender = KeyExtender(bmRes.keyMeta, blocksAfter=1)
        bmRes.counts = mapCounts(bmRes.counts, extender)
        bmRes.keyMeta = extender.meta()
                    
        bpRes.blocks = bmRes.blocks = bmRes.blocks + tuple([bpRes.blocks[1]])

        return super(UncorrectedTeleport, self)._convolve(results, noiseModels, pauli)

#    def _postCount(self, result):
#        # TODO: permute the output blocks to match outblockOrder
#        raise NotImplementedError

    class BMPropagator(KeyManipulator):
        
        def __init__(self, splitter, bm):
            super(UncorrectedTeleport.BMPropagator, self).__init__(splitter)
            bmMeta, self._outMeta = splitter.meta()
            self._bmPropagator = bm.keyPropagator(bmMeta)
            
        def meta(self):
            return (self._bmPropagator.meta(), self._outMeta)
        
        def _manipulate(self, key):
            keyBM, keyOut = key
            return (self._bmPropagator(keyBM), keyOut)
    
class TeleportED(InputDependentComponent):
    
    ucTeleportName = 'uct'
    inName = 'input'
    
    def __init__(self, kGood, bellPair, bellMeas):
        ucTeleport = UncorrectedTeleport(kGood, bellPair, bellMeas)
        super(TeleportED, self).__init__(kGood, subcomponents={self.ucTeleportName: ucTeleport})
        self._outBlock = bellPair.outBlocks()[1]
    
    def prAccept(self, noiseModels, kMax=None):
        # TODO: Currently only using the full XZ counts.  This
        # is still valid since Pr[bad] is also calculated with XZ here.
        # However, we'll want to use the X-only and Z-only counts,
        # as well.
        rejected = self.count(noiseModels[Pauli.Y], Pauli.Y).rejected
        prBad = self.prBad(noiseModels[Pauli.Y], Pauli.Y, kMax)
        locTotals = self.locations(Pauli.Y).getTotals()
        
        prSelf = 1 - probability.upperBoundPoly(rejected, prBad, locTotals, noiseModels[Pauli.Y])
        
        return prSelf * super(TeleportED, self).prAccept(noiseModels, kMax)
    
    def inBlocks(self):
        return self[self.ucTeleportName].inBlocks()
    
    def outBlocks(self):
        return (self._outBlock,)
    
    def _propagateInput(self, inputResult):
        ucTeleport = self[self.ucTeleportName]
        counts, keyMeta = ucTeleport.propagateCounts(inputResult.counts, inputResult.keyMeta)
                
        return counts, keyMeta
    
    def _postCount(self, result, noiseModels, pauli):
        '''
        Post-select for the trivial syndrome on the two measured blocks.
        '''
        result.counts, result.rejected = self._postselect(result)
        result.keyMeta = KeySplitter(result.keyMeta, [2]).meta()[1]
        result.blocks = (result.blocks[2],)
        return result
        
    def _postselect(self, result):
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
        
        # The UncorrectedTeleport component blocks are setup as follows:
        # block 0 - Transversal X-basis measurement
        # block 1 - Transversal Z-basis measurement
        # block 2 - Teleported data
        def accept(key):
            # TODO also check additional normalizers.
            return not(key[0] & rejectMask) and not(key[1] & rejectMask)
        
        accepted = []
        rejected = []
        for count in counts:
            acceptedK = {}
            rejectedK = {}
            for key, c in count.iteritems():
                reducedKey = keyForBlock(key, 2, keyMeta)
                if accept(key):
                    acceptedK[reducedKey] = c
                else:
                    rejectedK[reducedKey] = c
                    
            accepted.append(acceptedK)
            rejected.append(rejectedK)
            
        return accepted, rejected
