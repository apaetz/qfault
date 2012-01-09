'''
Created on May 1, 2011

@author: Adam
'''
from warnings import warn

class Block(object):
	'''
	classdocs
	'''


	def __init__(self, name, code, keyMeta=None):
		'''
		Constructor
		'''
		self.code = code
		self.name = name
		
	def length(self):
		return self.code.n
	
	def getCode(self):
		return self.code
	
	def __str__(self):
		return self.name + '.' + str(self.code)
	
	def __repr__(self):
		return 'Block(' + self.name + ', ' + str(self.code) + ')'
	
	def __eq__(self, other):
		if self.code == other.code:
			return True
		
		return False