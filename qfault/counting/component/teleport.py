'''
Created on 2011-10-25

Teleportation-based components.

@author: adam
'''
from qfault.counting import key
from qfault.counting.block import Block
from qfault.counting.component.base import PostselectionFilter, Empty, \
    ParallelComponent, SequentialComponent
from qfault.counting.component.transversal import TransRest
from qfault.counting.count_errors import mapCounts
from qfault.counting.count_parallel import convolve
from qfault.counting.key import keyForBlock, KeyExtender, KeyManipulator, KeyCopier, \
    IdentityManipulator, SyndromeKeyGenerator, KeyRemover
from qfault.counting.result import TrivialResult
from qfault.qec.error import xType, zType, Pauli
from qfault.util import bits
from qfault.util.cache import memoize
import logging
from copy import copy
from qfault.counting.component.block import BlockDiscard, BlockInsert

logger = logging.getLogger('component')
       
        
class TeleportWithMeas(SequentialComponent):
    '''
    Teleportation component.  The output of this component is three blocks, rather than one.
    The first two blocks contain the residual errors resulting from the Bell measurement.
    The third block is the teleportation output (i.e., the teleported block).  This component
    *does* perform the logical Pauli corrections on the output block.
    '''

    
    def __init__(self, kGood, bellPair, bellMeas, enableRest=True):

        # Order of the Bell pair and Bell measurement output blocks
        # is hardcoded.
        measXBlock, measZBlock = bellMeas.outBlocks()
        outBlock = bellPair.outBlocks()[1]
        
        emptyMeasX = Empty(measXBlock.getCode(), measXBlock.name)
        emptyMeasZ = Empty(measZBlock.getCode(), measZBlock.name)
        emptyOut = Empty(outBlock.getCode(), outBlock.name)
        
        
#        # Extend the input
#        extender = BlockPrepend(bellMeas.inBlocks()[:1], bellPair.inBlocks())
#        self._extender = extender
        
        # Extend the Bell measurement and the Bell pair so that they
        # contain the same number of blocks, some of which are empty.
        bellPair = ParallelComponent(kGood, emptyMeasX, bellPair)
        bellMeas = ParallelComponent(kGood, bellMeas, emptyOut)
        
        subs = (bellPair, bellMeas)
        
        # TODO: Bell measurement is two time steps.  Need to include
        # both rests??
        if enableRest:
            outBlock = bellPair.outBlocks()[1]
            rest = TransRest(kGood, outBlock.getCode())
            
            # Extend the rest to include the other two blocks, as well.
            rest = ParallelComponent(kGood, emptyMeasX, emptyMeasZ, rest)
            
            subs += (rest,)
        
        super(TeleportWithMeas, self).__init__(kGood, subcomponents=subs)

    def inBlocks(self):
        return (self.subcomponents()[1].inBlocks()[0],)
    
    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
        # First, extend the input to all three blocks.
        # TODO: there are two input extension functions.  This is because there seems to be a problem
        # when trying to propagate the input through the Bell Pair (which is only necessary when using
        # a concatenated code at level-two or higher).

        extender = BlockInsert(self.inBlocks(), 1, self[0].outBlocks()[1:])
        extendedInput = extender.count(inputResult=inputResult)
#        
#        self._log(logging.DEBUG, "extended input before input propagation: {0}".format(extendedInput))
#        
#        
        # Now propagate the input through the Bell measurement and transversal rest.
        for sub in self.subcomponents()[1:]:
            extendedInput = sub.propagateCounts(extendedInput)
#  
#        self._log(logging.DEBUG, "extended input after input propagation: {0}".format(extendedInput))
          
        # Now count normally.
        result = self._countInternal(noiseModels, pauli, kMax)
        
        self._log(logging.DEBUG, "internal result: {0}".format(result))
  
        
        # TODO: this pattern of memoizing some internal result and then convolving is
        # likely to be present in other components.  Should generalize this behavior.
        # TODO: more robust way of getting key lengths?
        keyLengths = [len(SyndromeKeyGenerator(block.getCode(), None).parityChecks()) for block in result.blocks]
        inKeyLengths = [len(SyndromeKeyGenerator(block.getCode(), None).parityChecks()) for block in extendedInput.blocks]
        self._log(logging.DEBUG, "keyLengths={0}, inKeyLengths={1}".format(keyLengths, inKeyLengths))
        result.counts = convolve(extendedInput.counts, 
                                result.counts, 
                                kMax=kMax, 
                                convolveFcn=key.convolveKeyCounts, 
                                extraArgs=[inKeyLengths, keyLengths])
        
        self._log(logging.DEBUG, "result before corrections: {0}".format(result))
        
        # Finally, make the logical corrections necessary for teleportation.
        inBlock = self.inBlocks()[0]
        code = inBlock.getCode()
        corrector = self._corrector(code)
        result.counts = mapCounts(result.counts, corrector)
        result.blocks = extendedInput.blocks
        
        self._log(logging.DEBUG, "result after corrections: {0}".format(result))
        
        return result
    
    def prAccept(self, noiseModels, inputResult=None, kMax=None):
        # The Bell-pair preparation may be non-deterministic, but we don't
        # expect the Bell measurement to be non-deterministic.
        return self[0].prAccept(noiseModels)
    
#    def _extendInputA(self, inputResult):
#        bp = self[0]
#        nBP = len(bp.inBlocks())
#        nIn = len(self.inBlocks())
#        self._log(logging.DEBUG, 'Extending input by {0} blocks after block {1}'.format(2, 1))
#        extender = KeyExtender(IdentityManipulator(), 2, 1)
#        extendedInput = inputResult
#        extendedInput.counts = mapCounts(inputResult.counts, extender)
#        extendedInput.blocks = inputResult.blocks[:1]*(3) + inputResult.blocks[1:]
#        return extendedInput
#
#    def _extendInputB(self, inputResult):
#        bp = self[0]
#        nBP = len(bp.inBlocks())
#        nIn = len(self.inBlocks())
#        self._log(logging.DEBUG, 'Extending input by {0} blocks after block {1}'.format(nBP - nIn, nIn))
#        extender = KeyExtender(IdentityManipulator(), nBP - nIn, nIn)
#        extendedInput = inputResult
#        extendedInput.counts = mapCounts(inputResult.counts, extender)
#        extendedInput.blocks = inputResult.blocks[:1]*(nBP) + inputResult.blocks[1:]
#        return extendedInput

    
    @memoize
    def _countInternal(self, noiseModels, pauli, kMax):
        '''
        Counts the internal components with the trivial input.
        '''
        extender = BlockInsert(self.inBlocks(), 1, self[0].inBlocks()[1:])
        inputResult = extender.count(inputResult=TrivialResult(self.inBlocks()))
#        inputResult = self._extendInputB(inputResult)
        return super(TeleportWithMeas, self).count(noiseModels, pauli, inputResult, kMax)
        
    
#    def keyPropagator(self, keyMeta):
#        # Extend the input to all three blocks, then propagate
#        # through the Bell measurement.
#        extender = KeyExtender(keyMeta, blocksAfter=2)
#        return self[self.bmName].keyPropagator(extender)
        
    def _corrector(self, code):
        parityChecks = SyndromeKeyGenerator(code, None).parityChecks()
        
        # TODO: this assumes that the code contains only one logical qubit.
        logicalOperators = code.logicalOperators()[0]
        
        logicalBitsX =  [check == logicalOperators[xType] for check in parityChecks]
        logicalBitsZ =  [check == logicalOperators[zType] for check in parityChecks]
        
        logicalMaskX = bits.listToBits(logicalBitsX)
        logicalMaskZ = bits.listToBits(logicalBitsZ)
        
        # Apply the logical operator parity checks (i.e. the logical corrections)
        # to the output key
        copier = KeyCopier(IdentityManipulator(), 0, 2, logicalMaskX)
        copier = KeyCopier(copier, 1, 2, logicalMaskZ)
        def outputKey(key):
            key = copier(key)
            return key
        
        return outputKey
    
    
class Teleport(SequentialComponent):
    '''
    Teleportation component.  This is just like TeleportWithMeas, except that by default
    the measurement results are removed and only the teleported block remains.  The
    actual output block may be selected by setting `block_to_keep`.  See TeleportWithMeas
    for output block numbers.
    '''
    
    def __init__(self, kGood, bellPair, bellMeas, enableRest=True, output_block=2):
        teleportWM = TeleportWithMeas(kGood, bellPair, bellMeas, enableRest)
        remove_blocks = set(range(3))
        remove_blocks.remove(output_block)
        discard = BlockDiscard(teleportWM.outBlocks(), remove_blocks)
        super(Teleport, self).__init__(kGood, subcomponents=[teleportWM, discard])
        
#    def outBlocks(self):
#        return (self[0].outBlocks()[2],)
#    
#    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
#        
#        nblocks = len(inputResult.counts[0].keys()[0])
#        
#        result = super(Teleport, self).count(noiseModels, pauli, inputResult, kMax)
#        
#        # Now trace out the measurement results
#        tracer = self.KeyTracer()
#        result.counts = mapCounts(result.counts, tracer)
#        result.blocks = result.blocks[2:]
#        
#        return result
#    
#            
#    class KeyTracer(KeyManipulator):
#                
#        def _manipulate(self, key):
#            return keyForBlock(key, 2) + key[3:]
    
class TeleportED(SequentialComponent):
    '''
    Teleportation component with simultaneous error-detection and post-selection.
    The output of this component is a single teleported block.  Only those errors
    that pass error-detection are included in the output.
    '''
    
    
    def __init__(self, kGood, bellPair, bellMeas, enableRest=True, acceptFunction=None):
#        inBlock = bellMeas.inBlocks()[0]
        teleport = TeleportWithMeas(kGood, bellPair, bellMeas, enableRest)
        edFilter = TeleportEDFilter(teleport, acceptFunction)
        super(TeleportED, self).__init__(kGood, subcomponents=(teleport, edFilter))
#        self._outBlock = bellPair.outBlocks()[1]
        
        # Delegate to subcomponents
#        self.inBlocks = teleport.inBlocks
#        self.outBlocks = edFilter.outBlocks
        
    def prAccept(self, noiseModels, inputResult, kMax):
        inputResult = self[0].count(noiseModels, Pauli.Y, inputResult, kMax)
        return self[1].prAccept(noiseModels, inputResult, kMax)

    
class TeleportEDFilter(PostselectionFilter):
    '''
    Error-detection filter for the output of the TeleportWithMeas component.
    '''
    
    def __init__(self, teleport, acceptFunction=None):
        super(TeleportEDFilter, self).__init__()
        self.teleport = teleport
        
        # The QECC on the data block and the first ancilla block must match exactly.
        # This includes any gauge stabilizers.  The reason is that all of the stabilizer
        # generators are used in error detection, including gauge stabilizers.
        dataCode, anc1Code, _ = (block.getCode() for block in teleport.outBlocks())
        if not dataCode == anc1Code:
            raise Exception("QECC for data block ({0}) and ancilla block ({1}) are not equal.".format(dataCode, anc1Code))
        
        if None == acceptFunction:
            acceptFunction = self._defaultAcceptFunction()
            
        self.accept = self._keyAcceptor(acceptFunction, dataCode)
                 
    def inBlocks(self):
        return self.teleport.outBlocks()
        
    def outBlocks(self):
        return (self.inBlocks()[2],)
    
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        
        # Use the QECC of the input block for error detection.
        # This assumes that the code of the input block and the first
        # Bell-state ancilla are *exactly* the same, including all
        # gauge stabilizers.
#        code = self.inBlocks()[0].getCode()
#        accept = self._acceptor(code)
        keyAcceptor = self.KeyAcceptor(subPropagator, self.accept)
        return keyAcceptor
    
    def _keyAcceptor(self, acceptFunction, code):
        stabilizers = set(code.stabilizerGenerators())
        # TODO: more robust way to get parity checks?
        # TODO: clean up all of this.
        parityChecks = SyndromeKeyGenerator(code, None).parityChecks()
        syndromeBits = [(check in stabilizers) for check in parityChecks]
        l = len(syndromeBits)
        syndromeIndices = [l-i-1 for i in range(l) if syndromeBits[i]]
        syndromeIndices.reverse()
#        rejectMask = bits.listToBits(syndromeBits)
        
        # Expect all syndrome bits to be grouped together after the logical
        # operator bits.
        if len(syndromeIndices) == 0:
            mask = shift = 0
        else:
            firstBit = min(syndromeIndices)
            lastBit = max(syndromeIndices)
            if syndromeIndices != range(firstBit, lastBit+1):
                raise RuntimeError("Syndrome bits are not contiguous")
        
            shift = firstBit
            mask = bits.lsbMask(lastBit - firstBit + 1)
        
        # The UncorrectedTeleport component blocks are setup as follows:
        # block 0 - Transversal X-basis measurement
        # block 1 - Transversal Z-basis measurement
        # block 2 - Teleported data
        
        def accept(key):
            try:
                s0 = (key[0] >> shift) & mask
                s1 = (key[1] >> shift) & mask
            except:
                print key
            return acceptFunction(s0) and acceptFunction(s1)
        
        return accept
        
    def _defaultAcceptFunction(self):
        def accept(syndrome):
            return not syndrome

    
        return accept
        
    class KeyAcceptor(KeyManipulator):
        
        def __init__(self, manipulator, accept):
            super(TeleportEDFilter.KeyAcceptor, self).__init__(manipulator)
            self.accept = accept
        
        def _manipulate(self, key):
            if not self.accept(key):
                return None
            # Preserve/ignore keys for blocks outside of the input space
            return keyForBlock(key, 2) + key[3:]
        
        
class EDInputFilter(SequentialComponent):
    '''
    This component is similar to the TeleportED component, but it doesn't
    propagate the input through the teleporation gadget.  Rather, it filters
    the input by keeping only those errors that *would* have been accepted
    had they actually been propagated through the teleportation.  This is
    useful, for example, when counting errors in a Rectangle and then
    conditioning on acceptance of the trailing EC of the exRec.
    '''
    
    def __init__(self, ed):
        super(EDInputFilter, self).__init__(ed.kGood, subcomponents=(ed,))
        
    def outBlocks(self):
        return self.inBlocks()
        
    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
        teleport = self[0]
        numBlocks = len(teleport.inBlocks())
        
        # Copy the input block.  One of the copies will be propagated through the
        # teleportation ED, the other will be preserved (so long as the input
        # passes error detection).
        copier = KeyExtender(IdentityManipulator(), numBlocks, numBlocks)
        for block in range(numBlocks):
            copier = KeyCopier(copier, block, block+numBlocks)
        extendedInput = copy(inputResult)
        extendedInput.counts = mapCounts(inputResult.counts, copier)
        extendedInput.blocks = inputResult.blocks[:numBlocks] + inputResult.blocks
        
        result = super(EDInputFilter, self).count(noiseModels, pauli, extendedInput, kMax)
        
        # Now remove the output block of the teleportation.  We only want to keep the original
        # input.
        numBlocks = len(teleport.outBlocks())
        remover = KeyRemover(IdentityManipulator(), [1,2])
        result.counts = mapCounts(result.counts, remover)
        result.blocks = result.blocks[numBlocks:]
        
        return result