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
	
class CountedBlock(Block):
	
	def __init__(self, name, code, counts):
		super(CountedBlock, self).__init__(name, code)
		self._counts = counts
		
	def __get__(self, etype):
		return self._counts[etype]
	
	def counts(self):
		return self._counts
		