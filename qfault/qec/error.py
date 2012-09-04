'''
Created on Jan 20, 2010

@author: adam
'''
from qfault.util import bits, listutils
from qfault.util.bits import listToBits

__all__ = ['Pauli', 'PauliError']

_strLookup = ('I', 'Z', 'X', 'Y')

xType = 'X'
zType = 'Z'
_dualTable = {xType: zType,
             zType: xType}

def dualType(etype):
    '''
    Returns the dual of the given error type.
    >>> dualType(xType)
    'Z'
    >>> dualType(zType)
    'X'
    '''
    return _dualTable[etype]

class PauliError(object):
    '''
    An n-qubit Pauli error, for arbitrary n.
    Does not account for overall phase factors.
    
    This class does not impose a particular ordering/endianness
    (i.e, ascending or descending) on the qubits.
    However, to be consistent with class methods
    fromstring() and fromlist(), the caller should
    specify qubit zero as the MSB, and qubit (n-1) as the LSB.  
    This convention also allows the tensor product to be read-off 
    directly from left to right.
    
    >>> PauliError(3, 0b110, 0b011)
    XYZ
    '''
    
    __slots__ = ['ebits', 'length']
    
    # TODO: generalize for a variable number of arguments.
    @classmethod
    def Tensor(cls, a, b):
        '''
        Returns the tensor product of 'a' and 'b'.
        >>> str(PauliError.Tensor(Pauli.Z, Pauli.Y))
        'ZY'
        >>> str(PauliError.Tensor(Pauli.X, Pauli.I))
        'XI'
        '''
        shift = b.length
        return cls(a.length + b.length, 
                   (a.ebits[xType]<<shift) + b.ebits[xType],
                   (a.ebits[zType]<<shift) + b.ebits[zType])
    
    def __init__(self, length, xbits=0, zbits=0):
        self.ebits = { xType: xbits, zType: zbits }
        self.length = length
        
    @classmethod
    def fromstring(cls, string):
        '''
        Returns an instance of PauliError corresponding
        to the given string of Pauli operators.
        
        >>> PauliError.fromstring('XYZ')
        XYZ
        '''
        xlist = [pauli.lower() in ('x', 'y') for pauli in string]
        zlist = [pauli.lower() in ('z', 'y') for pauli in string]
        return cls.fromlist(xlist, zlist)
    
    @classmethod
    def fromlist(cls, xlist, zlist):
        '''
        Returns an instance of PauliError corresponding to
        the symplectic form implied by the two lists.
        
        >>> PauliError.fromlist([1,1,0], [0,1,1])
        XYZ
        '''
        assert len(xlist) == len(zlist)
        return cls(len(xlist), listToBits(xlist), listToBits(zlist))
        
    def weight(self):
        return bits.weight(self.ebits[xType] | self.ebits[zType], self.length)
        
    def types(self):
        return [etype for etype in [xType, zType] if self.ebits[etype]]
    
    def asList(self):
        '''Returns the error as a list of single-qubit Pauli errors.'''
        xlist = bits.bitsToList(self.ebits[xType], self.length)
        zlist = bits.bitsToList(self.ebits[zType], self.length)
        
        return [PauliError(1, xbits=xlist[i], zbits=zlist[i]) for i in range(self.length)]
    
    def support(self):
        '''
        Returns a tuple of qubit on which the error has support.
        >>> PauliError.fromstring('XIYIZ').support()
        (0, 2, 4)
        '''
        return tuple(i for i, pauli in enumerate(self.asList())
                     if Pauli.I != pauli)
    
    def commutesWith(self, other):
        '''
        Returns True if the commutator with 'other' is zero, and False otherwise.
        
        >>> YY = PauliError(2, 0b11,0b11)
        >>> YZ = PauliError(2, 0b10,0b11)
        >>> XY = PauliError(2, 0b11,0b01)
        >>> YY.commutesWith(XY)
        False
        >>> YZ.commutesWith(XY)
        True
        '''
        return not bits.parity(self.ebits[xType] & other.ebits[zType]) ^ bits.parity(self.ebits[zType] & other.ebits[xType])
        
    def partial(self, pauli):
        '''
        Returns just the part of the error that corresponds
        to the given Pauli operator.
        
        >>> xyz = PauliError.fromstring('XYZ')
        >>> xyz.partial(Pauli.X)
        XXI
        >>> xyz.partial(Pauli.Z)
        IZZ
        >>> xyz.partial(Pauli.Y)
        IYI
        '''
        if Pauli.X == pauli:
            return self.__class__(self.length, xbits=self.ebits[xType])
        elif Pauli.Z == pauli:
            return self.__class__(self.length, zbits=self.ebits[zType])
        elif Pauli.Y == pauli:
            bitand = self.ebits[xType] & self.ebits[zType]
            return self.__class__(self.length, bitand, bitand)
    
    def __getitem__(self, index):
        '''
        >>> PauliError.fromstring('XYZZ')[1]
        Y
        '''
        shift = self.length - index - 1
        return self.__class__(1, 
                              (self.ebits[xType] >> shift) & 1,
                              (self.ebits[zType] >> shift) & 1)
        
    def __setitem__(self, index, value):
        '''
        >>> xyz = PauliError.fromstring('XYZ')
        >>> xyz[2] = Pauli.X
        >>> xyz
        XYX
        '''
        assert 1 == len(value)
        
        bit_index = self.length - index - 1
        self.ebits[xType] = bits.setbit(bit_index, 
                                        self.ebits[xType],
                                        value.ebits[xType])
        self.ebits[zType] = bits.setbit(bit_index, 
                                        self.ebits[zType],
                                        value.ebits[zType])
    
    def __hash__(self):
        return (self.ebits[xType], self.ebits[zType]).__hash__()
        
    def __eq__(self, other):
        return ((self.ebits[xType] == other.ebits[xType]) and
                (self.ebits[zType] == other.ebits[zType]))
    
    def __ne__(self, other):
        return not self == other
    
    def __xor__(self, other):
        '''
        >>> str(PauliError(2, 0b11, 0b01) ^ PauliError(2, 0b01, 0b11))
        'YI'
        '''
        return PauliError(self.length,
                          self.ebits[xType] ^ other.ebits[xType], 
                          self.ebits[zType] ^ other.ebits[zType])

    def __ixor__(self, other):
        return self ^ other

    def __mul__(self, other):
        return self ^ other
    
    def __imul__(self, other):
        return self ^ other

    def __iadd__(self, other):
        return self + other
    
    # TODO: is there a good reason to overload '+' as the tensor operation.
    # A more clear alternative would be to define a static method that accepts
    # two arguments.
    def __add__(self, other):
        '''
        Concatenate with another error.
        This is equivalent of taking the tensor product.
        >>> str(PauliError(1,0,1) + PauliError(1,1,1))
        'ZY'
        >>> str(PauliError(1,1,0) + PauliError(1,0,0))
        'XI'
        '''
        return self.Tensor(self, other)
        
    def __pow__(self, exp):
        '''
        Returns the tensor product of the error with itself 'exp' times.
        >>> str(PauliError(1,xbits=1) ** 3)
        'XXX'
        '''
        return sum([self]*(exp-1), self)
        
    def __rshift__(self, i):
        return PauliError(self.length - i, xbits=self[xType] >> i, zbits=self.ebits[zType] >> i)

    def __lshift__(self, i):
        return PauliError(self.length + i, xbits=self[xType] << i, zbits=self.ebits[zType] << i)

        
    def __len__(self):
        return self.length
    
    def __str__(self):
        X = self.ebits[xType]
        Z = self.ebits[zType]
        s = []
        for _ in range(self.length):
            s += _strLookup[2*(X&1) + int(Z&1)]
            X >>= 1
            Z >>= 1
            
        return ''.join(reversed(s))
    
    def __repr__(self):
        return str(self)
        #return 'PauliError({0:b},{1:b})'.format(self[xType], self[zType])

class _Pauli(object):
    '''
    Single-qubit Pauli operator.
    '''

    @property
    def I(self):
        return PauliError(1)
    @property
    def X(self):
        return PauliError(1, 1, 0)
    @property
    def Z(self):
        return PauliError(1, 0, 1)
    @property
    def Y(self):
        return PauliError(1, 1, 1)
    
    def __init__(self):
        self._dualTable = {
                      self.I: self.I,
                      self.X: self.Z,
                      self.Z: self.X,
                      self.Y: self.Y
                     }
        
        self._bitTable = [self.I, self.X, self.Z, self.Y]
        
#    @staticmethod
    def Dual(self, e):
        '''
        Returns the dual of the given Pauli operator
        (ignoring possible phase factors).
        >>> str(Pauli.Dual(Pauli.X))
        'Z'
        '''
        return self._dualTable[e]
    
#    @classmethod
    def fromint(self, bits):
        '''
        Returns the Pauli operator corresponding to the given
        2-bit value.  The MSB represents Pauli.Z, the LSB
        represents Pauli.X.
        >>> Pauli.fromint(0)
        I
        >>> Pauli.fromint(1)
        X
        >>> Pauli.fromint(2)
        Z
        >>> Pauli.fromint(3)
        Y
        '''
        return self._bitTable[bits]
    
# Public alias for _Pauli.  We need an actual
# object for the @property decorators to work correctly.
Pauli = _Pauli()        
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
        