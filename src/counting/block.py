'''
Created on May 1, 2011

@author: Adam
'''

class Block(object):
	'''
	classdocs
	'''


	def __init__(self, name, code):
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