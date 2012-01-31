'''
Created on 2011-11-14

@author: adam
'''
from counting.block import Block
from counting.component.base import Filter
from counting.key import IdentityManipulator, \
    KeyManipulator, SyndromeKeyDecoder
from qec import qecc


class IdealDecoder(Filter):
    '''
    An ideal decoder.
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