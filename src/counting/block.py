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
	'''
	An encoded block, or possibly multiple encoded blocks, along with its
	associated error counts.
	There are two ways to initialize: as a single block and an
	error correcting code, or from multiple sub-blocks.
	
	Example (single counted block):
	CountedBlock(counts, code=someCode, name='some block')
	
	Example (from multiple sub-blocks):
	CountedBlock(counts, subblocks=(block1, block2))
	
	Order is important when initializing from multiple sub-blocks.  In particular,
	the order given by subblocks should match the order in each of the counts.  In
	the example above, count errors should be of the format (error on block1, error on block2).
	
	Also note that CountedBlock(counts, subblocks=block1) is equivalent to
	CountedBlock(counts, code=block1.getCode(), name=block1.name) 
	'''
	
	def __init__(self, counts, subblocks=(), code=None, name=None):
		
		if 0 < len(subblocks):
			name = subblocks[0].name.join('.' + block.name for block in subblocks[1:])
			
			if 1 == len(subblocks):
				code = subblocks[0].getCode()
			else:
				if None != code:
					raise Exception('Code {0} is given, but the number of subblocks is nonzero.')

				code = tuple(block.getCode() for block in subblocks)
			
		super(CountedBlock, self).__init__(name, code)
			
		self._subblocks = subblocks
		self._counts = counts
		
	def __get__(self, etype):
		return self._counts[etype]
	
	def length(self):
		return sum(block.getCode() for block in self._subblocks)
	
	def counts(self):
		return self._counts
	
	def subblocks(self):
		return self._subblocks
		