'''
Some simple functions for counting and manipulating bits.
Some code was copied and adapted from 
http://wiki.python.org/moin/BitManipulation.

TODO: speed up with Cython or pyrex?

Created on 2011-08-29
@author: adam
'''
import warnings
#import util.cython.bits as cbits

def weight(integer, n=0):
    '''
    Returns the Hamming weight of 'integer'.
    
    Note: bitlength argument (n) is deprecated.
    
    >>> weight(0b10101)
    3
    '''
    if __debug__ and 0 != n:
        warnings.warn('Bitlength argument is deprecated', 
                      category=DeprecationWarning)
        
    s = 0
    while(integer):
        integer &= integer - 1  # Clears the LSB
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
        warnings.warn('Bitlength argument is deprecated', 
                      category=DeprecationWarning)

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
    
    warnings.warn('Deprecated.  Use x.bit_length(), instead.', 
                  category=DeprecationWarning)
    
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
    
def bitsToList(bits, n, bigEndian=True):
    '''
    >>> bitsToList(0b01101, 5)
    [0, 1, 1, 0, 1]
    >>> bitsToList(0b01101, 5, bigEndian=False)
    [1, 0, 1, 1, 0]
    '''
    blist = [0]*n
    indices = range(n)
    if bigEndian:
        indices = reversed(indices)
        
    for i in indices:
        blist[i] = bits & 1
        bits >>= 1
        
    return blist

def endianSwap(bits, n):
    '''
    >>> endianSwap(0b1101, 4) # 1101 -> 1011
    11
    >>> endianSwap(0b1011, 4) # 1011 -> 1101
    13
    '''
        
    return listToBits(bitsToList(bits, n, bigEndian=False))
    
def split(bits, lengths, reverse=False):
    '''
    Split the given value into a tuple according to the given bit lengths.
    The tuple ordering is big endian (the most significant bits are mapped
    to the first item in the tuple).
    >>> split(0b110100, [1,2,3])
    (1, 2, 4)
    >>> split(0b110100, [3,2,1], reverse=True)
    (4, 2, 1)
    '''
    
    n = len(lengths)
    items = [0] * n
    indices = range(n)
    if not reverse:
        indices = reversed(indices)
    
    for i in indices:
        l = lengths[i]
        items[i] = bits & ((1 << l) - 1)
        bits >>= l
        
    return tuple(items)

def concatenate(seq, lengths, reverse=False):
    '''
    Concatenate the sequence of values according to the given bit lengths.
    By default the resulting value is big endian (the first item of the sequence maps
    to the most significant bits of the result).  If reverse=True, then the
    resulting value is little endian (first item -> lsb)
    >>> cat = concatenate([1, 2, 3], [1, 3, 2])
    >>> '{0:b}'.format(cat)
    '101011'
    >>> cat = concatenate([1, 2, 3], [1, 3, 2], reverse=True)
    >>> '{0:b}'.format(cat)
    '110101'
    '''
    indices = range(len(seq))
    if reverse:
        indices = reversed(indices)
        
    bits = 0
    for i in indices:
        bits = (bits << lengths[i]) + seq[i]
    
    return bits

def lsbMask(n):
    '''
    Returns a mask that selects the n least significant bits.
    >>> mask = lsbMask(5)
    >>> '{0:b}'.format(mask)
    '11111'
    '''
    return (1 << n) - 1

def setbit(bit_index, val, new_bit_val):
    '''
    Returns the integer value in which the 'bit_index' bit
    of 'val' has been set to 'new_bit_val'.  The
    'bit_index' parameter is zero-based, i.e., zero
    corresponds to the least-significant bit.
    
    >>> val = 0b1101
    >>> setbit(2, val, 0)
    9
    >>> setbit(1, val, 1)
    15
    '''
    assert (1 == new_bit_val) or (0 == new_bit_val)
    
    # Save the bits below the bit to be set
    mask = lsbMask(bit_index)
    lower_bits = val & mask
    
    # Clear out the lower bits
    mask <<= 1
    mask += 1
    val ^= val & mask
    
    return val + (new_bit_val << bit_index) + lower_bits
    
    
    
    

if __name__ == '__main__':
    import doctest
    doctest.testmod()