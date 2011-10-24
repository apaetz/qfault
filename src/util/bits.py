'''
Some simple functions for counting and manipulating bits.
Some code was copied and adapted from http://wiki.python.org/moin/BitManipulation.

TODO: speed up with Cython or pyrex?

Created on 2011-08-29
@author: adam
'''
import warnings
#import util.cython.bits as cbits

def weight(x, n=0):
    '''
    Returns the Hamming weight x.
    
    Note: bitlength argument (n) is deprecated.
    
    >>> weight(0b10101)
    3
    '''
    if __debug__ and 0 != n:
        warnings.warn('Bitlength argument is deprecated', category=DeprecationWarning)
        
    s = 0
    while(x):
        x &= x - 1  # Clears the LSB
        s += 1
    return(s)


def parity(e, n = 0):
    '''
    Returns the parity (sum modulo 2) of x.
    
    >>> parity(0b1101)
    1
    >>> parity(0b110101)
    0
    '''
    if __debug__ and 0 != n:
        warnings.warn('Bitlength argument is deprecated', category=DeprecationWarning)

    return weight(e) & 1

def getBits(e, n = 32):
    bits = [0] * weight(e, n)
    j = 0
    for i in range(n+1):
        if e & 1:
            bits[j] = i
            j += 1
        e >>= 1
    return bits

def numbits(x):
    '''
    Returns the (base-one) index of the most significant 1-bit of x.
    
    >>> numbits(0)
    0
    >>> numbits(1<<42)
    43
    '''
    
    warnings.warn('Deprecated.  Use x.bit_length(), instead.', category=DeprecationWarning)
    
    i = 0
    while x > 0:
        x >>= 1
        i+=1
    return i

def bitsToStr(e, n = 32):
    """Converts the (32-bit) integer e into a string of bits.
    
    >>> bitsToStr(11, 6)
    '001011'
    """
    bits = ['0'] * n
    for i in range(n):
        if (e>>i)&1: 
            bits[n-1-i] = '1'
    return ''.join(bits)

def listToBits(values):
    '''
    >>> listToBits([1,0,1,1])
    11
    '''
    return reduce(lambda s,b: (s << 1) + bool(b), values, 0)
#    s = 0
#    for b in values:
#        s = (s<<1) + bool(b)
#    return s
    #return cbits.cython_listToBits(values)
    
def split(bits, lengths):
    '''
    Split the given value into a tuple according to the given bit lengths.
    The tuple ordering is big endian (the most significant bits are mapped
    to the first item in the tuple).
    >>> split(0b110100, [1,2,3])
    (1, 2, 4)
    '''
    n = len(lengths)
    items = [0] * n
    for i in reversed(range(n)):
        l = lengths[i]
        items[i] = bits & ((1 << l) - 1)
        bits >>= l
        
    return tuple(items)

def concatenate(seq, lengths):
    '''
    Concatenate the sequence of values according to the given bit lengths.
    The resulting value is big endian (the first item of the sequence maps
    to the most significant bits of the result).
    >>> cat = concatenate([1, 2, 3], [1, 3, 2])
    >>> '{0:b}'.format(cat)
    '100110'
    '''
    bits = seq[0]
    for i,l in enumerate(lengths[1:]):
        bits = (bits << l) + seq[i]
    return bits

if __name__ == '__main__':
    import doctest
    doctest.testmod()