'''
Created on 2011-10-25

@author: adam
'''
from copy import copy
from counting import probability, key
from counting.block import Block
from counting.component.base import Component, InputDependentComponent,\
    PostselectingComponent, Empty, Concatenator
from counting.countErrors import mapCounts
from counting.key import keyForBlock, KeyExtender, KeySplitter, KeyManipulator, \
    KeyConcatenator, IntegerKey, KeyCopier, KeyMasker, SyndromeKeyMeta
from counting.result import CountResult
from qec.error import Pauli, xType, zType
from util import bits
from util.cache import fetchable
import counting.key
import logging
import warnings
from counting.component import transversal
from counting.component.transversal import TransRest
from util.bits import listToBits

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
        bellPair = Concatenator(kGood, emptyMeasX, bellPair)
        bellMeas = Concatenator(kGood, bellMeas, emptyOut)
        
        subs = {self.bpName: bellPair,
                self.bmName: bellMeas}
        
        if enableRest:
            outBlock = bellPair.outBlocks()[1]
            measOut = TransRest(kGood, outBlock.getCode())
            
            # Extend the measurement to include the other two blocks, as well.
            measOut = Concatenator(kGood, emptyMeasX, emptyMeasZ, measOut)
            
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

class UCTSyndromeOut(UncorrectedTeleport):
    
    def outBlocks(self):
        return (self[self.bpName].outBlocks()[1], )
    
    def _postCount(self, result, noiseModels, pauli):
        code = self.outBlocks()[0].getCode()
        stabilizers = code.stabilizerGenerators()
        parityChecks = result.keyMeta.parityChecks()
        syndromeBits = [(check in stabilizers) for check in parityChecks]
        syndromeMask = listToBits(syndromeBits)
        
        reducer = self.KeyReducer(result.keyMeta, syndromeMask)
        result.counts = mapCounts(result.counts, reducer)
        result.keyMeta = reducer.meta()
        result.blocks = (result.blocks[2], )
        
        return result
        
    class KeyReducer(KeyManipulator):
        
        def __init__(self, meta, syndromeMask):
            super(UCTSyndromeOut.KeyReducer, self).__init__(meta)
            self.syndromeMask = syndromeMask
            self.masker = KeyMasker(self.meta(), syndromeMask)
        
        def meta(self):
            meta = super(UCTSyndromeOut.KeyReducer, self).meta()
            return SyndromeKeyMeta(meta.parityChecks(), 1)
        
        def _manipulate(self, key):
            key = keyForBlock(key, 2, self.meta())
            return self.masker(key)
            
    
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
    
    def _prepareInput(self, inputResult):
        ucTeleport = self[self.ucTeleportName]
        counts, keyMeta = ucTeleport.propagateCounts(inputResult.counts, inputResult.keyMeta)
                
        return CountResult(counts, keyMeta, None, inputResult.rejected)
    
    def _postCount(self, result, noiseModels, pauli):
        '''
        Post-select for the trivial syndrome on the two measured blocks.
        '''
        
        keyMeta = result.keyMeta
        outBlock = result.blocks[2]
        stabilizers = set(outBlock.getCode().stabilizerGenerators())
        logicals = outBlock.getCode().logicalOperators()
        if 1 != len(logicals):
            raise Exception('Incorrect number of logical qubits ({0})'.format(len(logicals)))
        
        logicals = logicals[0]
        
        accept = self._acceptor(keyMeta, stabilizers)
        outputKey = self._outputter(keyMeta, logicals)
        
        result.counts, result.rejected = self._postselect(result.counts, pauli, accept, outputKey)
        
        result.keyMeta = KeySplitter(result.keyMeta, [2]).meta()[1]
        result.blocks = (result.blocks[2],)
        return result
    
    def _acceptor(self, keyMeta, stabilizers):
        parityChecks = keyMeta.parityChecks()
        syndromeBits = [(check in stabilizers) for check in parityChecks]
        rejectMask = bits.listToBits(syndromeBits)
        
        warnings.warn('Detection of normalizer generators is not yet implemented')
        # The UncorrectedTeleport component blocks are setup as follows:
        # block 0 - Transversal X-basis measurement
        # block 1 - Transversal Z-basis measurement
        # block 2 - Teleported data
        def accept(key):
            # TODO also check additional normalizers.
            return not(key[0] & rejectMask) and not(key[1] & rejectMask)
    
        return accept
    
    def _outputter(self, keyMeta, logicalOperators):
        parityChecks = keyMeta.parityChecks()
        
        logicalBitsX =  [check == logicalOperators[xType] for check in parityChecks]
        logicalBitsZ =  [check == logicalOperators[zType] for check in parityChecks]
        
        logicalMaskX = bits.listToBits(logicalBitsX)
        logicalMaskZ = bits.listToBits(logicalBitsZ)
        
        # Apply the logical operator parity checks (i.e. the logical corrections)
        # to the output key
        copier = KeyCopier(keyMeta, 0, 2, logicalMaskX)
        copier = KeyCopier(copier, 1, 2, logicalMaskZ)
        def outputKey(key):
            key = copier(key)
            return keyForBlock(key, 2, keyMeta)
        
        return outputKey
        
    def _postselect(self, counts, pauli, accept, outputKey):
        '''
        Post-select for the trivial syndrome on the two measured blocks, and make
        logical corrections if necessary.
        '''
        
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
                    self._log(logging.DEBUG, 'accepted %s, output=%s, count=%d', key, outKey, c)
                    acceptedK[outKey] = acceptedK.get(outKey, 0) + c
                    
                    if rejectAll and (0 != k):
                        rejectedK += c
                else:
                    self._log(logging.DEBUG, 'rejected %s', key)
                    rejectedK += c
                    
            accepted.append(acceptedK)
            rejected.append({counting.key.rejectKey: rejectedK})
            
        return accepted, rejected
    
class TeleportEDFilter(TeleportED):
    '''
    Teleportation ED gadget that filters out rejected inputs.
    This ED doesn't actually propagate the input through the teleportation component.
    Rather, it simulates the action of the teleportation and outputs the original input
    only if it was accepted.
    '''
    
    def _convolve(self, results, noiseModels, pauli, inputResult):
        
        # We don't want any errors from the teleportation itself to go to the output.
        # Zero-out all of the syndromes on the output block of the UCT.
        uct = results[self.ucTeleportName]
        masker = KeyMasker(uct.keyMeta, 0, [2])
        uct.counts = mapCounts(uct.counts, masker)
        uct.keyMeta = masker.meta()
        
        return super(TeleportEDFilter, self)._convolve(results, noiseModels, pauli, inputResult)
    
    def _outputter(self, keyMeta, logicalOperators):
        
        def outputKey(key):
            return keyForBlock(key, 2, keyMeta)
        
        return outputKey
    
    def _prepareInput(self, inputResult):
        # 1. Extend the input to all three blocks, 
        # 2. copy from the original input block to the output block.
        # 3. Propagate through the Bell measurement  
        extender = KeyExtender(inputResult.keyMeta, blocksAfter=2)
        copier = KeyCopier(extender, 0, 2)
        uct = self[self.ucTeleportName]
        propagator = uct[uct.bmName].keyPropagator(copier)
        counts = mapCounts(inputResult.counts, propagator)
                
        return CountResult(counts, propagator.meta(), None, inputResult.rejected)