'''
Created on 2011-11-14

@author: adam
'''
from counting.component.base import Component
from counting.key import SyndromeKeyGenerator, SyndromeKeyMeta
from counting.result import CountResult


class ComponentAdapter(Component):
    
    def __init__(self, adaptee):
        kGood = {}
        super(ComponentAdapter, self).__init__(kGood, subcomponents={'adaptee': adaptee})
        
        self.locations = adaptee.locations
        self.inBlocks = adaptee.inBlocks
        self.outBlocks = adaptee.outBlocks
        self.logicalStabilizers = adaptee.logicalStabilizers
        self.propagateCounts = adaptee.propagateCounts
        self.keyPropagator = adaptee.keyPropagator
        self.prBad = adaptee.prBad
        
    
class InputAdapter(ComponentAdapter):
    
    def __init__(self, component, inputKey):
        blocks = component.inBlocks()
        code = blocks[0].getCode()
        parityChecks = SyndromeKeyGenerator.ParityChecks(code)
        keyMeta = SyndromeKeyMeta(parityChecks, len(blocks))
        inputResult = CountResult([{inputKey: 1}], keyMeta, blocks)
        
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