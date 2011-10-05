'''
Created on Jan 20, 2010

@author: adam
'''
from util import bits

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
	>>> str(PauliError(0b110, 0b011))
	'XYZ'
	'''
	
	__slots__ = ['ebits']
	
	def __init__(self, xbits=0, zbits=0):
		self.ebits = { xType: xbits, zType: zbits }
		
	def types(self):
		return [etype for etype in [xType, zType] if self[etype]]
		
	def __getitem__(self, pauli):
		'''
		>>> PauliError(0b10, 0b01)[xType]
		2
		'''
		return self.ebits[pauli]
	
	def __hash__(self):
		return (self.ebits[xType], self.ebits[zType]).__hash__()
		
	def __eq__(self, other):
		return (self[xType] == other[xType]) and (self[zType] == other[zType])
	
	def __ne__(self, other):
		return not self == other
	
	def __xor__(self, other):
		'''
		>>> str(PauliError(0b11, 0b01) ^ PauliError(0b01, 0b11))
		'YI'
		'''
		return PauliError(self[xType] ^ other[xType], 
						  self[zType] ^ other[zType])

	def __ixor__(self, other):
		return self ^ other

	def __mul__(self, other):
		return self ^ other
	
	def __imul__(self, other):
		return self ^ other

	def __iadd__(self, other):
		return self + other
	
	def __add__(self, other):
		'''
		Concatenate with another error.
		This is equivalent of taking the tensor product.
		>>> str(PauliError(0,1) + PauliError(1,1))
		'ZY'
		>>> str(PauliError(1,0) + PauliError(0,0))
		'XI'
		'''
		shift = max(bits.numbits(other[xType]), bits.numbits(other[zType]), 1)
		return PauliError((self[xType]<<shift) + other[xType],
						  (self[zType]<<shift) + other[zType])
		
	def __str__(self):
		X = self[xType]
		Z = self[zType]
		s = []
		while (X > 0) or (Z > 0):
			s += _strLookup[2*(X&1) + int(Z&1)]
			X >>= 1
			Z >>= 1
			
		if 0 == len(s):
			s = ['I']
		return ''.join(reversed(s))
	
	def __repr__(self):
		return str(self)
		#return 'PauliError({0:b},{1:b})'.format(self[xType], self[zType])

class Pauli(object):
	'''
	'''

	I = PauliError()
	X = PauliError(1,0)
	Z = PauliError(0,1)
	Y = PauliError(1,1)
	
	_dualTable = {
				  I: I,
				  X: Z,
				  Z: X,
				  Y: Y
				 }
	
	_bitTable = [I, Z, X, Y]
	
	@staticmethod
	def Dual(e):
		'''
		Returns the dual of the given Pauli operator
		(ignoring possible phase factors).
		>>> str(Pauli.Dual(Pauli.X))
		'Z'
		'''
		return Pauli._dualTable[e]
	
	@staticmethod
	def FromBits(bits):
		'''
		Returns the Pauli operator corresponding to the given
		2-bit value.  The MSB represents Pauli.Z, the LSB
		represents Pauli.X.
		>>> str(Pauli.FromBits(3))
		'Y'
		'''
		return Pauli._bitTable[bits]
	

		
#class CompoundError(object):
#	
#	def __init__(self, e1, e2):
#		self.error = (e1, e2)
#	
#	def __xor__(self, other):
#		e1, e2 = np.bitwise_xor(self.error, other.error)
#		return CompoundError(e1, e2)
#	
#	def __add__(self, other):
#		return self.__xor__(other)
#	
#	def __getitem__(self, index):
#		return self.error[index]
#	
#	def __str__(self):
#		return str(self.error)
	
	
if __name__ == "__main__":
	import doctest
	doctest.testmod()
		