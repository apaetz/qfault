'''
Created on 2011-10-22

@author: adam
'''    

from qfault.counting import key
import logging

logger = logging.getLogger('counting.result')

class CountResult(object):
    ''' 
    Container for counting results.
    '''
    
    def __init__(self, counts, blocks):                    
        self.blocks = blocks
        self.counts = counts
        
    def is_valid(self, expNumBlocks=None):
        nblocks = len(self.blocks)
        
        logger.debug("nblocks={0}, expected={1}".format(nblocks, expNumBlocks))
        if None != expNumBlocks and nblocks != expNumBlocks:
            logger.error('nblocks={0}, expNumBlocks={1}'.format(nblocks, expNumBlocks))
            return False
        for count in self.counts:
            if any(nblocks - len(key) for key in count.keys()):
                logger.error('nblocks={0}, key lengths={1}'.format(nblocks, 
                                                                   [len(key) for key in count.keys()]))
                logger.debug('count={0}'.format(count))
                return False
            
        return True

        
    def __get__(self, k):
        return self.counts[k]

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self.counts) + str(self.blocks)


def TrivialResult(blocks):
    inputs = tuple([0]*len(blocks))
    counts = [{inputs: 1}]
    return CountResult(counts, blocks)