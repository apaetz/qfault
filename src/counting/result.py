'''
Created on 2011-10-22

@author: adam
'''    

from counting import key

class CountResult(object):
    ''' 
    Container for counting results.
    '''
    
    def __init__(self, counts, blocks, rejectedCounts=None, name=None):
                    
        self.blocks = blocks
        self.counts = counts
        if None == rejectedCounts:
            # TODO hack
            rejectedCounts = [{key.rejectKey:0}]
            
        self.rejected = rejectedCounts
        
    def __get__(self, etype):
        return self.counts[etype]

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self.counts) + str(self.blocks)


def TrivialResult(blocks):
    inputs = tuple([0]*len(blocks))
    counts = [{inputs: 1}]
    return CountResult(counts, blocks)