'''
Created on 2011-10-23

@author: adam
'''
import logging

logging.basicConfig(level=logging.INFO,
                    format='[%(asctime)s] %(process)d %(levelname)-8s %(name)s: %(message)s')

logger = logging.getLogger('counting.convolve')

def convolveDict(counts1, counts2):
    '''
    Convolve counts from two dictionaries.
    
    >>> a = {0: 1, 1: 2}
    >>> b = {0: 1, 2: 3}
    >>> convolveDict(a, b)
    {0: 1, 1: 2, 2: 3, 3: 6}
    '''
    counts = {}
    logger.info('convolving dictionaries {0}x{1}'.format(len(counts2), len(counts1)))
    for key2, count2 in counts2.iteritems():
        for key1, count1 in counts1.iteritems():
            key = key1 ^ key2
            
            #TODO: not sure if it would be faster to use a try-except block here, instead.
            counts[key] = counts.get(key, 0) + (count2 * count1)
    
    return counts

if __name__ == '__main__':
    import doctest
    doctest.testmod()