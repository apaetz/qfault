'''
Created on Jan 20, 2010

@author: adam
'''
from util import bits, listutils, iteration

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
	However, a good convention is to specify qubit
	zero as the MSB, and qubit (n-1) as the LSB.  This
	allows the tensor product to be read-off directly
	from left to right.
	
	>>> str(PauliError(3,0b110, 0b011))
	'XYZ'
	'''
	
	__slots__ = ['ebits', 'length']
	
	# TODO: generalize for a variable number of arguments.
	@staticmethod
	def Tensor(a, b):
		'''
		Returns the tensor product of 'a' and 'b'.
		>>> str(PauliError.Tensor(Pauli.Z, Pauli.Y))
		'ZY'
		>>> str(PauliError.Tensor(Pauli.X, Pauli.I))
		'XI'
		'''
		shift = b.length
		return PauliError(a.length + b.length, 
						  (a[xType]<<shift) + b[xType],
						  (a[zType]<<shift) + b[zType])
	
	def __init__(self, length, xbits=0, zbits=0):
		self.ebits = { xType: xbits, zType: zbits }
		self.length = length
		
	def types(self):
		return [etype for etype in [xType, zType] if self[etype]]
	
	def asList(self):
		'''Returns the error as a list of single-qubit Pauli errors.'''
		xlist = bits.bitsToList(self[xType], self.length)
		zlist = bits.bitsToList(self[zType], self.length)
		
		return [PauliError(1, xbits=xlist[i], zbits=zlist[i]) for i in range(self.length)]
	
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
		return not bits.parity(self[xType] & other[zType]) ^ bits.parity(self[zType] & other[xType])
		
	def __getitem__(self, pauli):
		'''
		>>> PauliError(2, 0b10, 0b01)[xType]
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
		>>> str(PauliError(2, 0b11, 0b01) ^ PauliError(2, 0b01, 0b11))
		'YI'
		'''
		return PauliError(self.length,
						  self[xType] ^ other[xType], 
						  self[zType] ^ other[zType])

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
		return PauliError(self.length - i, xbits=self[xType] >> i, zbits=self[zType] >> i)

	def __lshift__(self, i):
		return PauliError(self.length + i, xbits=self[xType] << i, zbits=self[zType] << i)

		
	def __str__(self):
		X = self[xType]
		Z = self[zType]
		s = []
		for _ in range(self.length):
			s += _strLookup[2*(X&1) + int(Z&1)]
			X >>= 1
			Z >>= 1
			
		return ''.join(reversed(s))
	
	def __repr__(self):
		return str(self)
		#return 'PauliError({0:b},{1:b})'.format(self[xType], self[zType])

class Pauli(object):
	'''
	'''

	I = PauliError(1)
	X = PauliError(1,1,0)
	Z = PauliError(1,0,1)
	Y = PauliError(1,1,1)
	
	_dualTable = {
				  I: I,
				  X: Z,
				  Z: X,
				  Y: Y
				 }
	
	_bitTable = [I, X, Z, Y]
	
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
		>>> str(Pauli.FromBits(0))
		'I'
		>>> str(Pauli.FromBits(1))
		'X'
		>>> str(Pauli.FromBits(2))
		'Z'
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
		