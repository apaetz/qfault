'''
Created on 2011-08-29

@author: adam
'''

def weight(e, n = 32):
    """Returns the Hamming weight of the 32-bit integer e.
    
    >>> weight(3)
    2
    """
    s = 0
    for _ in range(n):
        if (e & 1):
            s = s+1
        e = e>>1
    return s

def parity(e, n = 32):
    return weight(e,n) & 1

def getBits(e, n = 32):
    bits = [0] * weight(e, n)
    j = 0
    for i in range(n+1):
        if e & 1:
            bits[j] = i
            j += 1
        e >>= 1
    return bits

def bitsToStr(e, n = 32):
    """Converts the (32-bit) integer e into a string of bits.
    
    >>> etostr(11, 6)
    '001011'
    """
    bits = ['0'] * n
    for i in range(n):
        if (e>>i)&1: 
            bits[n-1-i] = '1'
    return ''.join(bits)