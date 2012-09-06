'''
Created on 2011-10-23

@author: adam
'''
from qfault.counting import count_locations
from qfault.util import bits
from qfault.util.concurrency import mapreduce_concurrent
from qfault.util.iteration import PartitionIterator
import logging
import operator

logger = logging.getLogger('counting.convolve')

def convolve_dict(counts1, 
                  counts2, 
                  keyOp=operator.xor, 
                  countMul=operator.mul, 
                  countAdd=operator.add, 
                  nullCount=0):
    '''
    Convolve counts from two dictionaries.
    Convolution takes every possible combination of a single
    element from counts1 and a single element from counts2.
    For each pair, the counts are multiplied using the
    countMul operator, and the keys
    are combined using keyOp operator.  If two pairs
    yield the same key, the counts from each pair is added
    using the countAdd operator.
    
    >>> a = {0: 1, 1: 2}
    >>> b = {0: 1, 2: 3}
    >>> convolve_dict(a, b)
    {0: 1, 1: 2, 2: 3, 3: 6}
    '''
    counts = {}
    logger.debug('convolving dictionaries {0}x{1}'.format(len(counts2), len(counts1)))
    for key2, count2 in counts2.iteritems():
        for key1, count1 in counts1.iteritems():
            key = keyOp(key1, key2)
            
            #TODO: not sure if it would be faster to use a try-except block here, instead
            # of counts.get()
            counts[key] = countAdd(counts.get(key, nullCount), countMul(count1, count2))
    
    return counts


class ConvolveCaller(object):
    
    def __init__(self, counts1, counts2, convolve_fcn):
        self._counts1 = counts1
        self._counts2 = counts2
        self._convolve_fcn = convolve_fcn
        
    def __call__(self, k1k2):
        k1, k2 = k1k2
        results = self._convolve_fcn(self._counts1[k1], self._counts2[k2])
        return self._convolve_fcn(self._counts1[k1], self._counts2[k2])

def convolve_counts(counts1, 
                    counts2, 
                    k0_max=None, 
                    k1_max=None, 
                    k_max=None, 
                    convolve_fcn=convolve_dict):
    '''
    >>> counts1 = [{0: 1, 1: 2}, {0: 3, 2: 1}]
    >>> counts2 = [{0: 1, 2: 3}]
    >>> convolve_counts(counts1, counts2)
    [{0: 1, 1: 2, 2: 3, 3: 6}, {0: 6, 2: 10}]
    '''
    
    if None == k0_max:
        k0_max = len(counts1)-1
    if None == k1_max:
        k1_max = len(counts2)-1
    if None == k_max:
        k_max = k0_max + k1_max
        
    k_max = min(k_max, len(counts1) + len(counts2) - 2)

    # Distribute the work
    convolved = []
    map_func = ConvolveCaller(counts1, counts2, convolve_fcn)
    for k in range(k_max+1):        
        convolved.append(mapreduce_concurrent(map_func, 
                                         count_locations.merge_counts, 
                                         PartitionIterator(k, 2, [k0_max, k1_max])))
    
    return convolved
    
def convolve_dict_tuples(bitlengths1, bitlengths2, counts1, counts2):
    '''
    Convolve dictionaries for which the keys are stored as tuples 
    of integers.
    Key operations (XOR) are performed element-wise on each tuple.
    
    The tuples in each of the counts need not be of the same length.
    However, the integer bit lengths of the tuples must match, up
    to the length of the shorter tuple.  More precisely,
    bitlengths1[i] == bitlengths2[i] for all 
    i < min(len(bitlengths1, bitlengths2)).
    
    If len(bitlengths1) != len(bitlengths2), then the resulting
    counts will contain tuples of length 
    min(len(bitlengths1, bitlengths2)).
    
    >>> counts1 = {(0,1): 1, (1,0): 1}
    >>> counts2 = {(1,1): 1, (2,0): 2}
    >>> convolve_dict_tuples(counts1, counts2, [2, 1], [2, 1])
    {(3, 0): 2, (0, 1): 1, (1, 0): 1, (2, 1): 2}
    '''
    
    # Return counts according to the larger of the two tuples.
    if len(bitlengths2) > len(bitlengths1):
        bitlengths = bitlengths2
        min_tuple_len = len(bitlengths1)
    else:
        bitlengths = bitlengths1
        min_tuple_len = len(bitlengths2)
    
    if bitlengths1[:min_tuple_len] != bitlengths2[:min_tuple_len]:
        raise Exception('Incompatible key lengths {0}, {1}'.format(bitlengths1, bitlengths2))

    counts1 = {bits.concatenate(key, bitlengths1, reverse=True): count 
               for key, count in counts1.iteritems()}
    counts2 = {bits.concatenate(key, bitlengths2, reverse=True): count 
               for key, count in counts2.iteritems()}
    
    counts = convolve_dict(counts1, counts2)
    
    return {bits.split(key, bitlengths, reverse=True): count 
            for key, count in counts.iteritems()}
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()