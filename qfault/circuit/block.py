'''
Created on May 1, 2011

@author: Adam
'''

class Block(object):
	'''
	Represents a single encoded block.
	'''

	def __init__(self, name, code):
		'''
		Constructor
		'''
		self.code = code
		self.name = name
		
	def length(self):
		return self.code.n
	
	def get_code(self):
		return self.code
		
	def __str__(self):
		return self.name + '.' + str(self.code)
	
	def __repr__(self):
		return 'Block(' + self.name + ', ' + str(self.code) + ')'
	
	def __eq__(self, other):
		if self.length == other.length and self.code == other.code:
			return True
		
		return False