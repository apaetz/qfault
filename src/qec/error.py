'''
Created on Jan 20, 2010

@author: adam
'''

class Pauli(object):
	'''
	classdocs
	'''

	I = 'I'
	X = 'X'
	Z = 'Z'
	Y = 'Y'
	
class Error(object):
	
	def __init__(self, xbits=0, zbits=0):
		self.xbits = xbits
		self.zbits = zbits
		
	def getBits(self, pauli):
		
		
	def __str__(self):
		return ''.join(self.paulis)
		
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
	
		