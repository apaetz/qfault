'''
Created on 2012-04-03

@author: adam
'''
from counting.component.base import Filter
from counting.key import IdentityManipulator, KeyPermuter, KeyRemover
from util import listutils

class BlockPermutation(Filter):
    '''
    Noiselessly permutes the logical blocks according to the given permutation.
    Permutation is given in one-line notation.  For example, the permutation
    [0,2,1,3] on (a,b,c,d) yields (a,c,b,d).
    '''
    
    def __init__(self, inBlocks, permutation):
        self.perm = permutation
        self._inBlocks = inBlocks
        super(BlockPermutation, self).__init__()
        
    def inBlocks(self):
        return self._inBlocks
        
    def outBlocks(self):
        return tuple(listutils.permute(self.inBlocks(), self.perm))
        
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        return KeyPermuter(subPropagator, self.perm)
    
    
class BlockDiscard(Filter):
    '''
    Traces out the last N blocks.
    '''
    
    def __init__(self, inBlocks, n):
        super(BlockDiscard, self).__init__()
        self._inBlocks = inBlocks
        self.n = n
        
    def inBlocks(self):
        return self._inBlocks
    
    def outBlocks(self):
        return tuple(self.inBlocks()[:-self.n])
        
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        l = len(self.inBlocks())
        return KeyRemover(subPropagator, l-self.n, l)