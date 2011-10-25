'''
Created on 2011-10-25

@author: adam
'''
from copy import copy
from counting import probability
from counting.component.base import Component
from counting.countErrors import mapCounts, extendCounts
from counting.key import keySplitter, keyConcatenator, keyForBlock
from qec.error import Pauli
from util import bits
from util.cache import fetchable
import logging

logger = logging.getLogger('component')

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
    @fetchable
    def _convolveBP(bellMeas, bellPairResult, bpMeasBlock):
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
        
    def _convolve(self, results, noiseModels, pauli):
        results[self.bpName] = self._convolveBP(self.subcomponents()[self.bmName],
                                                results[self.bpName], 
                                                self._bpMeasBlock)
        
        inRes = copy(self._data)
        # Propagate the input through the CNOT of the Bell measurement.
        bellMeas = self.subcomponents()[self.bmName]
        inRes.counts, inRes.keyMeta = bellMeas.propagateCounts(self._data.counts, 
                                                                         self._data.keyMeta,
                                                                         bellMeas.measXName)
        #self._data.blocks = results[self.bmName].blocks            
        results['input'] = inRes
        
        # Extend each of the results over all three blocks.
        bpRes = results[self.bpName]
        bpBlocksAfter = self._bpMeasBlock

        bmRes = results[self.bmName]
        bmRes.counts, bmRes.keyMeta = extendCounts(bmRes.counts, bmRes.keyMeta, not bpBlocksAfter, bpBlocksAfter)

        inRes.counts, inRes.keyMeta = extendCounts(inRes.counts, inRes.keyMeta, not bpBlocksAfter, bpBlocksAfter)
        
        if bpBlocksAfter:
            blocks = bpRes.blocks + inRes.blocks
        else:
            blocks = inRes.blocks + bpRes.blocks
            
        bpRes.blocks = bmRes.blocks = inRes.blocks = blocks

        return super(Teleport, self)._convolve(results, pauli)

    def _postCount(self, result):
        raise NotImplementedError
    
class TeleportED(Teleport):
    
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
    
    def _postCount(self, result, noiseModels, pauli):
        '''
        Post-select for the trivial syndrome on the two measured blocks.
        '''
        
        # TODO: compute and apply XZ corrections
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
        
        # The Teleport component blocks are setup as follows:
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
