'''
Created on 2011-10-23

@author: adam
'''
import logging
import operator

logger = logging.getLogger('counting.convolve')

def convolveDict(counts1, counts2, keyOp=operator.xor, countMul=operator.mul, countAdd=operator.add, nullCount=0):
    '''
    Convolve counts from two dictionaries.
    
    >>> a = {0: 1, 1: 2}
    >>> b = {0: 1, 2: 3}
    >>> convolveDict(a, b)
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

if __name__ == '__main__':
    import doctest
    doctest.testmod()