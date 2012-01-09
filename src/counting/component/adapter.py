'''
Created on 2011-11-14

@author: adam
'''
from counting.block import Block
from counting.component.base import Component, Filter
from counting.key import SyndromeKeyGenerator, IdentityManipulator, \
    KeyManipulator, SyndromeKeyDecoder
from counting.result import CountResult
from qec import qecc


class IdealDecoder(Filter):
    
    def __init__(self, code):
        super(IdealDecoder, self).__init__()
        self._code = code
        
    def inBlocks(self):
        return (Block('', self._code),)
    
    def outBlocks(self):
        return (Block('', qecc.TrivialStablizerCode()),)
    
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        return self.KeyDecoder(self._code, subPropagator)
    
    class KeyDecoder(KeyManipulator):
        
        def __init__(self, code, manipulator):
            super(IdealDecoder.KeyDecoder, self).__init__(manipulator)
            self.decoder = SyndromeKeyDecoder(code)
    
        def _manipulate(self, key):
            decoded = self.decoder.decode(key[0])
            return (decoded,) + key[1:]

class ComponentAdapter(Component):
    
    def __init__(self, adaptee):
        kGood = {}
        super(ComponentAdapter, self).__init__(kGood, subcomponents=(adaptee,))
        
        self.locations = adaptee.locations
        self.inBlocks = adaptee.inBlocks
        self.outBlocks = adaptee.outBlocks
        self.logicalStabilizers = adaptee.logicalStabilizers
        self.propagateCounts = adaptee.propagateCounts
        self.keyPropagator = adaptee.keyPropagator
        self.prBad = adaptee.prBad
        
    
class InputAdapter(ComponentAdapter):
    
    def __init__(self, component, inputCounts):
        # TODO: kludgey. need to clean up.
        if type(inputCounts) == tuple:
            inputCounts = [{inputCounts: 1}]
        
        blocks = component.inBlocks()
        code = blocks[0].getCode()
        parityChecks = SyndromeKeyGenerator.ParityChecks(code)
        keyMeta = SyndromeKeyMeta(parityChecks, len(blocks))
        inputResult = CountResult(inputCounts, blocks)
        
        self._component = component
        self._inResult = inputResult
        super(InputAdapter, self).__init__(component)
        
 
        
    def count(self, noiseModels, pauli):
        return self._component.count(noiseModels, pauli, self._inResult)
    
    def prAccept(self, noiseModels, kMax=None):
        return self._component.prAccept(noiseModels, self._inResult, kMax=kMax)
    
    def _hashStr(self):
        return super(InputAdapter, self)._hashStr() + str(self._inResult)
    
class NonzeroSyndromeOutputAdapter(ComponentAdapter):
    
    def _postCount(self, result, noiseModels, pauli):
        # hack!
        zeroKey = tuple([0]*result.keyMeta.nblocks)
        
        # Remove the trivial syndrome counts
        for count in result.counts:
            count.pop(zeroKey, 0)
        
        return self.adaptee._postCount(result, noiseModels, pauli)