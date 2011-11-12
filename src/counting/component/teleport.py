'''
Created on 2011-10-25

@author: adam
'''
from copy import copy
from counting import probability, key
from counting.block import Block
from counting.component.base import Component, InputDependentComponent,\
    PostselectingComponent, Empty, ConcatenatedComponent
from counting.countErrors import mapCounts
from counting.key import keyForBlock, KeyExtender, KeySplitter, KeyManipulator, \
    KeyConcatenator, IntegerKey, KeyCopier
from counting.result import CountResult
from qec.error import Pauli, xType, zType
from util import bits
from util.cache import fetchable
import counting.key
import logging
import warnings
from counting.component import transversal
from counting.component.transversal import TransRest

logger = logging.getLogger('component')

class UncorrectedTeleport(Component):
    
    inName = 'input'
    bpName = 'BP'
    bmName = 'BM'
    measOutName = 'measOut'
    
    def __init__(self, kGood, bellPair, bellMeas, enableRest=True):

        # Order of the Bell pair and Bell measurement output blocks
        # is hardcoded.
        measXBlock, measZBlock = bellMeas.outBlocks()
        outBlock = bellPair.outBlocks()[1]
        
        emptyMeasX = Empty(measXBlock.getCode(), measXBlock.name)
        emptyMeasZ = Empty(measZBlock.getCode(), measZBlock.name)
        emptyOut = Empty(outBlock.getCode(), outBlock.name)
        
        # Extend the Bell measurement and the Bell pair so that they
        # contain the same number of blocks, some of which are empty.
        bellPair = ConcatenatedComponent(kGood, emptyMeasX, bellPair)
        bellMeas = ConcatenatedComponent(kGood, bellMeas, emptyOut)
        
        subs = {self.bpName: bellPair,
                self.bmName: bellMeas}
        
        if enableRest:
            outBlock = bellPair.outBlocks()[1]
            measOut = TransRest(kGood, outBlock.getCode())
            
            # Extend the measurement to include the other two blocks, as well.
            measOut = ConcatenatedComponent(kGood, emptyMeasX, emptyMeasZ, measOut)
            
            subs[self.measOutName] = measOut
        
        super(UncorrectedTeleport, self).__init__(kGood, subcomponents=subs)
        #self._outblockOrder = outblockOrder
        
    def inBlocks(self):
        return (self[self.bmName].inBlocks()[0], )
    
    def outBlocks(self):
        # It is expected that UncorrectedTeleport will be wrapped by another
        # component which will be able to define the output blocks.
        raise NotImplementedError
    
    def keyPropagator(self, keyMeta):
        # Extend the input to all three blocks, then propagate
        # through the Bell measurement.
        extender = KeyExtender(keyMeta, blocksAfter=2)
        return self[self.bmName].keyPropagator(extender)


    def _convolve(self, results, noiseModels, pauli):
        
        # Propagate the Bell pair through the Bell measurement.
        bm = self[self.bmName]
        bpResult = results[self.bpName]
        bpResult.counts, bpResult.keyMeta = bm.propagateCounts(bpResult.counts, bpResult.keyMeta)

        # Now convolve normally.
        return super(UncorrectedTeleport, self)._convolve(results, noiseModels, pauli)
    
#    def _propagateBP(self, bpCounts, bpMeta):
#        # Propagate the first half of the bell pair through the bell measurement.
#        extender = KeyExtender(bpMeta, blocksBefore=1)
#        splitter = KeySplitter(extender, [2])
#        bmPropagator = self.BMPropagator(splitter, self[self.bmName])
#        concatenator = KeyConcatenator(bmPropagator)
#
#        bpCounts = mapCounts(bpCounts, concatenator)
#        bpMeta = concatenator.meta()
#        
#        return bpCounts, bpMeta
#
#    class BMPropagator(KeyManipulator):
#        
#        def __init__(self, bm, bpMeta):
#            splitter = KeySplitter(bpMeta, [3])
#            super(UncorrectedTeleport.BMPropagator, self).__init__(splitter)
#            bmMeta, self._outMeta = splitter.meta()
#            self._bmPropagator = bm.keyPropagator(bmMeta)
#            
#        def meta(self):
#            return (self._bmPropagator.meta(), self._outMeta)
#        
#        def _manipulate(self, key):
#            keyBM, keyOut = key
#            return (self._bmPropagator(keyBM), keyOut)
    
class TeleportED(PostselectingComponent):
    
    ucTeleportName = 'uct'
    inName = 'input'
    
    def __init__(self, kGood, bellPair, bellMeas):
        ucTeleport = UncorrectedTeleport(kGood, bellPair, bellMeas)
        super(TeleportED, self).__init__(kGood, subcomponents={self.ucTeleportName: ucTeleport})
        self._outBlock = bellPair.outBlocks()[1]
        
        # Delegate to ucTeleport
        self.inBlocks = ucTeleport.inBlocks
            
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
        Post-select for the trivial syndrome on the two measured blocks, and make
        logical corrections if necessary.
        '''
        counts = result.counts
        keyMeta = result.keyMeta
        outBlock = result.blocks[2]
        stabilizers = set(outBlock.getCode().stabilizerGenerators())
        logicals = outBlock.getCode().logicalOperators()
        if 1 != len(logicals):
            raise Exception('Incorrect number of logical qubits ({0})'.format(len(logicals)))
        
        logicals = logicals[0]
        
        parityChecks = keyMeta.parityChecks()
        
        syndromeBits = [(check in stabilizers) for check in parityChecks]
        logicalBitsX =  [check == logicals[xType] for check in parityChecks]
        logicalBitsZ =  [check == logicals[zType] for check in parityChecks]
        
        rejectMask = bits.listToBits(syndromeBits)
        logicalMaskX = bits.listToBits(logicalBitsX)
        logicalMaskZ = bits.listToBits(logicalBitsZ)
        
        warnings.warn('Detection of normalizer generators is not yet implemented')
        # The UncorrectedTeleport component blocks are setup as follows:
        # block 0 - Transversal X-basis measurement
        # block 1 - Transversal Z-basis measurement
        # block 2 - Teleported data
        def accept(key):
            # TODO also check additional normalizers.
            return not(key[0] & rejectMask) and not(key[1] & rejectMask)
        
        # Apply the logical operator parity checks (i.e. the logical corrections)
        # to the output key
        copier = KeyCopier(keyMeta, 0, 2, logicalMaskX)
        copier = KeyCopier(copier, 1, 2, logicalMaskZ)
        def outputKey(key):
            key = copier(key)
            return keyForBlock(key, 2, keyMeta)
        
        accepted = []
        rejected = []
        
        # Error detection involves looking at both the X and Z
        # syndromes.  We can upper bound the rejection probability
        # only when counting both types of errors.  When counting
        # only X, or only Z, we must assume (for the purpose of the
        # rejected counts) that all errors are rejected, unless
        # k = 0.
        rejectAll = (Pauli.Y != pauli)
        
        for k, count in enumerate(counts):
            acceptedK = {}
            rejectedK = 0
            for key, c in count.iteritems():
                if accept(key):
                    outKey = outputKey(key)
                    acceptedK[outKey] = acceptedK.get(outKey, 0) + c
                    
                    if rejectAll and (0 != k):
                        rejectedK += c
                else:
                    rejectedK += c
                    
            accepted.append(acceptedK)
            rejected.append({counting.key.rejectKey: rejectedK})
            
        return accepted, rejected
