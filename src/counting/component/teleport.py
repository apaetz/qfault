'''
Created on 2011-10-25

@author: adam
'''
from copy import copy
from counting import probability
from counting.component.base import Component, InputDependentComponent
from counting.countErrors import mapCounts, extendCounts
from counting.key import keySplitter, keyConcatenator, keyForBlock, keyExtender
from qec.error import Pauli
from util import bits
from util.cache import fetchable
import logging

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
        
    def outBlocks(self):
        #TODO
        raise NotImplementedError
    
    def keyPropagator(self, keyMeta, blockname):
        # Propagate the input through the CNOT of the Bell measurement and
        # then extend to the output block.
        bellMeas = self[self.bmName]
        bmPropagator = bellMeas.keyPropagator(bellMeas.measXName)
        extender = keyExtender(keyMeta, blocksAfter=1)
        
        def propagate(key):
            return extender(bmPropagator(key))
        
        return propagate
        
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
        splitter, meta0, meta1 = keySplitter(bpMeta, 1)
        bmPropagator, bmMeta = bm.keyPropagator(meta0, bm.measZName)
        concatenator, propMeta = keyConcatenator(bmMeta, meta1)
        
        def propagator(key):
            bpKey0, bpKey1 = splitter(key)
            key = concatenator(bmPropagator(bpKey0), bpKey1)
            return key

        bpRes.counts = mapCounts(bpRes.counts, propagator)
        bpRes.keyMeta = propMeta
        
        # Extend bell measurement results over all three blocks.
        bmRes.counts, bmRes.keyMeta = extendCounts(bmRes.counts, bmRes.keyMeta, blocksAfter=1)
                    
        bpRes.blocks = bmRes.blocks = bmRes.blocks + tuple([bpRes.blocks[1]])

        return super(UncorrectedTeleport, self)._convolve(results, noiseModels, pauli)

#    def _postCount(self, result):
#        # TODO: permute the output blocks to match outblockOrder
#        raise NotImplementedError
    
class TeleportED(InputDependentComponent):
    
    ucTeleportName = 'uct'
    inName = 'input'
    
    def __init__(self, kGood, bellPair, bellMeas):
        ucTeleport = UncorrectedTeleport(kGood, bellPair, bellMeas)
        super(TeleportED, self).__init__(kGood, subcomponents={self.ucTeleportName: ucTeleport})
    
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
    
    def _propagateInput(self, inputResult):
        ucTeleport = self[self.ucTeleportName]
        counts, keyMeta = ucTeleport.propagateCounts(inputResult.counts, inputResult.keyMeta, ucTeleport.inName)
                
        return counts, keyMeta
    
    def _postCount(self, result, noiseModels, pauli):
        '''
        Post-select for the trivial syndrome on the two measured blocks.
        '''
        result.counts, result.rejected = self._postselect(result)
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
