'''
Created on 2012-04-03

@author: adam
'''
from qfault.counting.component.base import Filter
from qfault.counting.key import IdentityManipulator, KeyPermuter, KeyRemover,\
    KeyExtender
from qfault.util import listutils
from qfault.counting.result import CountResult

class BlockPermutation(Filter):
    '''
    Noiselessly permutes the logical blocks according to the given permutation.
    Permutation is given in one-line notation.  For example, the permutation
    [0,2,1,3] on (a,b,c,d) yields (a,c,b,d).
    '''
    
    def __init__(self, inBlocks, permutation):
        self.perm = permutation
        self._in_blocks = inBlocks
        super(BlockPermutation, self).__init__()
        
    def inBlocks(self):
        return self._in_blocks
        
    def outBlocks(self):
        return tuple(listutils.permute(self.inBlocks(), self.perm))
        
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        return KeyPermuter(subPropagator, self.perm)
    
    
class BlockDiscard(Filter):
    '''
    Traces out the given blocks.
    '''
    
    def __init__(self, in_blocks, blocks_to_remove):
        super(BlockDiscard, self).__init__()
        self._in_blocks = in_blocks
        self._blocks_to_remove = blocks_to_remove
        
    def inBlocks(self):
        return self._in_blocks
    
    def outBlocks(self):
        return tuple(listutils.remove_subsequence(self.inBlocks(), self._blocks_to_remove))
        
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        return KeyRemover(subPropagator, self._blocks_to_remove)
    
class BlockInsert(Filter):
    '''
    Prepends the given blocks to the existing blocks.
    '''
    
    def __init__(self, in_blocks, insert_before, blocks_to_insert):
        super(BlockInsert, self).__init__()
        self._in_blocks = in_blocks
        self._out_blocks = in_blocks[:insert_before] + blocks_to_insert + in_blocks[insert_before:]
        self._num_new_blocks = len(blocks_to_insert)
        self._insert_before = insert_before
        
    def inBlocks(self):
        return self._in_blocks

    def outBlocks(self):
        return self._out_blocks
    
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        return KeyExtender(subPropagator, self._num_new_blocks, self._insert_before)
    
    
class BlockCombine(Filter):
    '''
    Combines the input blocks into a single block.  The maximum count for each
    syndrome and fault order is chosen for output.
    '''
    
    def __init__(self, in_blocks):
        super(BlockCombine, self).__init__()
        self._in_blocks = in_blocks
        
    def inBlocks(self):
        return self._in_blocks
    
    def outBlocks(self):
        return self.inBlocks()[:1]
    
    def propagateCounts(self, inputResult):
        # Separate each of the input blocks
        inblock_indices = set(range(len(self.inBlocks())))
        block_counts = []
        for i, block in enumerate(self.inBlocks()):
            discard = BlockDiscard((block,), inblock_indices - set([i]))
            block_counts.append(discard.count(inputResult=inputResult).counts)
        
        # Take the maximum count for each fault order k and each syndrome.
        counts = count_errors.maxCount(*block_counts)
        
        return CountResult(counts, self.outBlocks() + inputResult.blocks[len(block_counts):])