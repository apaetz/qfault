'''
Created on 2011-11-14

This file contains components that adapt the output from other components in some way.

@author: adam
'''
from counting.block import Block
from counting.component.base import Filter, SequentialComponent,\
    ParallelComponent
from counting.key import IdentityManipulator, \
    KeyManipulator, SyndromeKeyDecoder, SyndromeKeyFilter
from qec import qecc


class IdealDecoder(Filter):
    '''
    An ideal decoder.
    Input syndrome keys are decoded to syndrome keys for the trivial code.
    TODO: should the keys be decoded directly to Paulis?
    '''
    
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
        
class FlaggingDecoder(IdealDecoder):
    '''
    An ideal decoder that also flags detectable errors.  The syndrome key
    for each count is transformed into a (flag, error) tuple.
    '''
    
    def __init__(self, code, decode_as_pauli=False):
        self._code = code
        self._decode_as_pauli = decode_as_pauli
    
    class KeyDecoder(KeyManipulator):
        
        def __init__(self, code, manipulator):
            super(FlaggingDecoder.KeyDecoder, self).__init__(manipulator)
            self._code = code
            self.decoder = SyndromeKeyDecoder(code)
    
        def _manipulate(self, key):
            syndrome = self.decoder.syndrome(key[0])
            flag = self._code.detectSyndrome(syndrome)
            decoded = self.decoder.decode(key[0])
            return ((flag, decoded),) + key[1:]
        
class SyndromeFilter(Filter):
    
    def __init__(self, code):
        super(SyndromeFilter, self).__init__()
        self._code = code
        
    def inBlocks(self):
        return (Block('', self._code),)
    
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        return SyndromeKeyFilter(self._code, subPropagator)
    
        
class DecodeAdapter(SequentialComponent):
    '''
    Applies ideal decoders to all output blocks of the given component.
    '''
    
    def __init__(self, component):
        outBlocks = component.outBlocks()
        idealDecoders = [IdealDecoder(block.getCode()) for block in outBlocks]
        decoder = ParallelComponent({}, *idealDecoders)
        super(DecodeAdapter, self).__init__(component.kGood, (component, decoder))
        
class SyndromeAdapter(SequentialComponent):
    
    def __init__(self, component):
        outBlocks = component.outBlocks()
        idealDecoders = [SyndromeFilter(block.getCode()) for block in outBlocks]
        decoder = ParallelComponent({}, *idealDecoders)
        super(SyndromeAdapter, self).__init__(component.kGood, (component, decoder))