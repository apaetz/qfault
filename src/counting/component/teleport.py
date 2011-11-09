'''
Created on 2011-10-25

@author: adam
'''
from copy import copy
from counting import probability, key
from counting.block import Block
from counting.component.base import Component, InputDependentComponent
from counting.countErrors import mapCounts
from counting.key import keyForBlock, KeyExtender, KeySplitter, KeyManipulator, \
    KeyConcatenator, IntegerKey
from counting.result import CountResult
from qec.error import Pauli
from util import bits
from util.cache import fetchable
import counting.key
import logging
import warnings

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
        
        bpRes.counts, bpRes.keyMeta = self._propagateBP(bpRes.counts, bpRes.keyMeta)
        bmRes.counts, bmRes.keyMeta = self._propagateBM(bmRes.counts, bmRes.keyMeta)

        bpRes.blocks = bmRes.blocks = bmRes.blocks + tuple([bpRes.blocks[1]])

        return super(UncorrectedTeleport, self)._convolve(results, noiseModels, pauli)
    
    def _propagateBP(self, bpCounts, bpMeta):
        # Propagate the first half of the bell pair through the bell measurement.
        extender = KeyExtender(bpMeta, blocksBefore=1)
        splitter = KeySplitter(extender, [2])
        bmPropagator = self.BMPropagator(splitter, self[self.bmName])
        concatenator = KeyConcatenator(bmPropagator)

        bpCounts = mapCounts(bpCounts, concatenator)
        bpMeta = concatenator.meta()
        
        return bpCounts, bpMeta
    
    def _propagateBM(self, bmCounts, bmMeta):
        # Extend bell measurement results over all three blocks.
        extender = KeyExtender(bmMeta, blocksAfter=1)
        bmCounts = mapCounts(bmCounts, extender)
        bmMeta = extender.meta()
        return bmCounts, bmMeta

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
        if None != kMax:
            kMaxY = kMax[Pauli.Y]
        else:
            kMaxY = None
            
        prBad = self.prBad(noiseModels[Pauli.Y], Pauli.Y, kMaxY)
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
                
        return CountResult(counts, keyMeta, self.outBlocks(), inputResult.rejected)
    
    def _postCount(self, result, noiseModels, pauli):
        '''
        Post-select for the trivial syndrome on the two measured blocks.
        '''
        result.counts, result.rejected = self._postselect(result, pauli)
        result.keyMeta = KeySplitter(result.keyMeta, [2]).meta()[1]
        result.blocks = (result.blocks[2],)
        return result
        
    def _postselect(self, result, pauli):
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
        
        warnings.warn('Detection of normalizer generators is not yet implemented')
        # The UncorrectedTeleport component blocks are setup as follows:
        # block 0 - Transversal X-basis measurement
        # block 1 - Transversal Z-basis measurement
        # block 2 - Teleported data
        def accept(key):
            # TODO also check additional normalizers.
            return not(key[0] & rejectMask) and not(key[1] & rejectMask)
        
        accepted = []
        rejected = []
        
        # Error detection involves looking at both the X and Z
        # syndromes.  We can upper bound the rejection probability
        # only when counting both types of errors.  When counting
        # only X, or only Z, we must assume (for the purpose of the
        # rejected counts) that all errors are rejected.
        rejectAll = (Pauli.Y != pauli)
        
        for count in counts:
            acceptedK = {}
            rejectedK = 0
            for key, c in count.iteritems():
                reducedKey = keyForBlock(key, 2, keyMeta)
                if accept(key):
                    acceptedK[reducedKey] = c
                    if rejectAll:
                        rejectedK += c
                else:
                    rejectedK += c
                    
            accepted.append(acceptedK)
            rejected.append({counting.key.rejectKey: rejectedK})
            
        return accepted, rejected
