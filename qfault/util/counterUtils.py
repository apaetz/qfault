#!/usr/bin/python
# by Ben Reichardt, created 12/22/2009
# 
# Some utility functions for manipulating bit strings.  
#
# It might be useful to have code that prints the different errors using I,X,Y,Z.
# 
from fractions import Fraction
from qfault.util.listutils import addLists, nonZeroIndices
import fractions
import logging

logger = logging.getLogger('counterUtils')

def weight(e, n = 32):
	"""Returns the Hamming weight of the 32-bit integer e.
	
	>>> weight(3)
	2
	"""
	sum = 0
	for _ in range(n):
		if (e & 1):
			sum = sum+1
		e = e>>1
	return sum

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

def lcm(i1, i2):
	'''
	Returns the least common multiple of integers i1, i2
	>>> lcm(6,8)
	24
	'''
	return i1*i2/fractions.gcd(i1, i2)

def findFactorMultiple(a,b):
	'''
	Finds a common factor delta of rational numbers a and b such that 
	a=delta*alpha and b=delta*beta where delta is a rational number and alpha, beta are integers.
	
	Returns [alpha, beta, delta]
	
	This is useful, for example, when adding two sets of normalized counts for which the 
	normalization factors (a and b) differ.  Then, use 'a*c1 + b*c1 = delta * (alpha*c1 + beta*c2)'
	to preserve the integer counts.

	>>> from fractions import Fraction
	>>> findFactorMultiple(Fraction(1,2), Fraction(3,4))
	[2, 3, Fraction(1, 4)]
	'''
	
	# No need to go through all of the calculations if a and b are equal.
	if a == b:
		return [a, b, 1]

	if type(a) is not Fraction:
		a = Fraction(a, 1)
		
	if type(b) is not Fraction:
		b = Fraction(b, 1)
	
	gcdNum = fractions.gcd(a.numerator, b.numerator)
	lcmDen = lcm(a.denominator, b.denominator)
	delta = Fraction(gcdNum, lcmDen)
	
	alpha = (a.numerator/gcdNum) * (lcmDen/a.denominator)
	beta = (b.numerator/gcdNum) * (lcmDen/b.denominator)
	
	return [alpha, beta, delta]

def etostr(e, n = 32):
	"""Converts the 32-bit integer e into a string of bits.
	
	>>> etostr(11, 6)
	'001011'
	"""
	bits = ['0'] * n
	for i in range(n):
		if (e>>i)&1: 
			bits[n-1-i] = '1'
	return ''.join(bits)

def readBinaryString(string): 
	"""Converts a string to binary, with '1', 'X', 'Y' and 'Z' mapped to 1 and all other characters mapped to 0.
	
	The bits are indexed from right to left.
	>>> readBinaryString('.....11.1')
	13
	"""
	err = 0
	for i in range(len(string)): 
		err = err<<1
		if string[i] in ('1', 'X', 'Y', 'Z'):
			err += 1
	return err

def exztostr(ex, ez, n = 32):
	"""Converts the bit strings (integers) ex and ez into a string of I,X,Y,Z Pauli errors.
	
	Note that qubit 0 is at the *right* end of the string.
	>>> exztostr(6, 3, 4)
	'IXYZ'
	"""
	bits = ['0'] * n
	for i in range(n):
		bits[n-1-i] = ('I', 'X', 'Z', 'Y')[(ez>>i&1)*2 + (ex>>i&1)]
	return ''.join(bits)

def errorsToStr(errors, n = 32, onlyNontrivial = False):
	"""Converts a combined bunch of errors to a readable IXYZ string per block.
	
	If the third argument is True, then the output string only includes blocks including nontrivial entries.  
	(Note that these are not reduced modulo the stabilizer, though.)
	>>> errorsToStr({'X':{'a':0, 'b':1, 'c':0}, 'Z':{'a':2, 'b':3, 'c':0}}, 2, True)
	'a:ZI, b:ZY'
	"""
	blocks = sorted(errors['X'].keys())
	if onlyNontrivial:
		blocks = [b for b in blocks if (errors['X'][b] or errors['Z'][b])]
	results = [b + ":" + exztostr(errors['X'][b], errors['Z'][b], n)
				for b in blocks]
	return ', '.join(results)

def errorsMerge(e1, e2):
	"""Merges the two sets of errors with disjoint keys.
	
	Returns a new dictionary, i.e., the merging is not done in place.  
	>>> errorsMerge({'X':{'a':0, 'b':1, 'c':0}, 'Z':{'a':2, 'b':3, 'c':0}}, {'X':{'d':4},'Z':{'d':5}})
	{'X': {'a': 0, 'c': 0, 'b': 1, 'd': 4}, 'Z': {'a': 2, 'c': 0, 'b': 3, 'd': 5}}
	"""
	import copy
	errors = copy.deepcopy(e1)
	for k in e2['X'].keys():
		errors['X'][k] = e2['X'][k]
		errors['Z'][k] = e2['Z'][k]
	return errors

def hashError(e, stabilizers):
	"""Returns a minimum-weight representative for the error e.
	
	The error is reduced modulo the stabilizers.  
	The stabilizers should be given as a list, like [[0,1,2],[0,3],...].  
	To speed this up, I might add a third argument shortCircuitWeight=0, and stop the error reduction 
	as soon as the weight reaches that threshold.
	
	>>> hashError(7, [[0,1],[2]])
	0
	>>> hashError(15, [[0,1,2],[0,1]])
	8
	"""
	if isinstance(stabilizers[0], list):
		stabilizers = [reduce(lambda x, y: x+y, 
							map(lambda x: 1<<x, z)) 
						for z in stabilizers]
	# the following two lines are equivalent to the above "list comprehension"
	#stabilizers = map(lambda z: map(lambda x: 1<<x, z), stabilizers)		# convert each stabilizer binary
	#stabilizers = map(lambda z: reduce(lambda x, y: x+y, z), stabilizers)	# sum each stabilizer up
	#print stabilizers
	bestErrorSoFar = e
	bestWeightSoFar = weight(e)
	for k in range(1, 1<<(len(stabilizers))):
		error = e
		for digit in range(len(stabilizers)):
			if ((k>>digit)&1): 
				error = error ^ stabilizers[digit]
		if weight(error) < bestWeightSoFar:	# or (trialweight == bestweightsofar and trialerror < besterrorsofar): 
			bestErrorSoFar = error
			bestWeightSoFar = weight(error)
	return bestErrorSoFar

# A location is like a struct.  I have implemented it as a dictionary, but could also have used a blank class, e.g., 
#class Employee:
#	pass
#john = Employee() # Create an empty employee record
#john.name = 'John Doe'
#john.dept = 'computer lab' #...

# I considered adding a loccopy(block1, block2).  This would copy X and Z errors from the first block 
# into the second, basically creating a snapshot of the computation that could later be referred back to.  
# But I didn't.  Instead of taking a snapshot of, say, a data block, I give the same data block two names, 
# for before and after.  Later I will manually exor the first part of it into the second.  

def locrest(block, bit):
	return {'type': 'rest', 'block1': block, 'bit1': bit}

def locXprep(block, bit):
	return {'type': 'prepX', 'block1': block, 'bit1': bit}

def locZprep(block, bit):
	return {'type': 'prepZ', 'block1': block, 'bit1': bit}

def loccnot(block1, bit1, block2, bit2):
	return {'type': 'cnot', 'block1': block1, 'bit1': bit1, 'block2': block2, 'bit2': bit2}

def locXmeas(block, bit):
	return {'type': 'measX', 'block1': block, 'bit1': bit}

def locZmeas(block, bit):
	return {'type': 'measZ', 'block1': block, 'bit1': bit}

def loctostr(loc, verbose = False, n = 32):
	"""Meant to convert a location to a string in a short, useful format.
	
	The verbose version is incomplete.  It needs to filter, so that it only shows nonzero 
	entries of the error dictionaries, and only dictionaries that have nonzero entries.
	"""
	#l = loc.copy()
	#del l['X1'], l['Z1'], 
	#if 'X2' in l:
	#	del l['X2'], l['Z2']
	#print l
	t = loc['type']
	if t == 'cnot':
		locstr = '%s(%s,%d,%s,%d)' % (t, loc['block1'], loc['bit1'], loc['block2'], loc['bit2'])
		if verbose:
			pass # DEBUG
	elif t in ['rest', 'prepX', 'prepZ', 'measX', 'measZ']:
		locstr = '%s(%s,%d)' % (t, loc['block1'], loc['bit1'])
		if verbose:
			locstr += " " + errorsToStr(loc['X1'], n, True)
			# DEBUG
	return locstr

def sumLocs(locList):
	'''
	Sum the number of locations by type.  Returns a dictionary index by type.
	'''
	sums = {}
	for l in locList:
		type = l['type']
		sums[type] = sums.get(type, 0) + 1
		
	return sums

def swapXZloc(loc):
	"""Returns a *copy* of the input location, conjugated by Hadamards."""
	l = loc.copy()
	#transposes = [0, 3, 6, 1, 4, 7, 2, 5, 8]
	#l['bit1'] = transposes[l['bit1']]
	#if 'bit2' in l:
	#	l['bit2'] = transposes[l['bit2']]
	t = l['type']
	t = t.replace('Z', 'ZZ').replace('X', 'XX')
	t = t.replace('ZZ', 'X').replace('XX', 'Z')
	l['type'] = t
	if l['type'] == 'cnot':
		l['block1'], l['block2'] = l['block2'], l['block1']
		l['bit1'], l['bit2'] = l['bit2'], l['bit1']
	return l

def swapXZlocations(locations):
	"""Conjugates the circuit by Hadamards."""
	return map(swapXZloc, locations)
	# I am following a convention of always copying a location before editing it, 
	# as though the location were immutable.  Thus I don't need to copy.deepcopy(locations).

def allBlocks(locations):
	"""Returns a set containing all the blocks used by locations in the input list."""
	blocks = set([])
	for loc in locations:
		blocks.add(loc['block1'])
		if 'block2' in loc: 
			blocks.add(loc['block2'])
	return blocks

def blockLengths(locations):
	"""Returns a dictionary containing length of each block used by locations in the input list."""
	blocks = {}
	for loc in locations:
		blockname = loc['block1']
		bit = loc['bit1']
		blocks[blockname] = max(blocks.get(blockname, 0), bit+1)
		if 'block2' in loc: 
			blockname = loc['block2']
			bit = loc['bit2']
			blocks[blockname] = max(blocks.get(blockname, 0), bit+1)
	return blocks


propagateNoOpTypes = set(['rest', 'prepX', 'prepZ', 'measX', 'measZ'])

def propagateErrorsThroughLocation(errors, loc):
	# currently only cnot location types affect the errors, but a natural addition might be Hadamards
	if loc['type'] in propagateNoOpTypes:
		return errors
	if loc['type'] == 'cnot':
		b1, b2 = errors['X'][loc['block1']], errors['X'][loc['block2']]
		b2 = b2 ^ ((b1>>loc['bit1']&1)<<loc['bit2'])
		errors['X'][loc['block2']] = b2
		b1, b2 = errors['Z'][loc['block1']], errors['Z'][loc['block2']]
		b1 = b1 ^ ((b2>>loc['bit2']&1)<<loc['bit1'])
		errors['Z'][loc['block1']] = b1
	return errors

def propagateErrors(errors, locations):
	"""Starting with the dictionary of errors, updates them according to the locations.
	
	errors should be arranged as {'X':{'blockID1':xerr1,'blockID2':xerr2,...}, 'Z':...}.
	Note that I am not taking advantage of the CSS structure.  Since X errors can only 
	transform to more X errors in a circuit lacking Hadamards, we could track X errors only 
	while propagating an initial X error.  Since this is all happening in precomputation, 
	though, the running time isn't constraining yet.  
	"""
	for loc in locations:
		errors = propagateErrorsThroughLocation(errors, loc)
	#blocks = errors['X'].keys()
	return errors

def propagateAllErrors(locations):
	"""Propagates all bit errors.
	Errors are constructed using a descending (little endian) bit ordering.
	A fault that propagates to a single X error on qubit k will results in a bit string
	in which all bits are zero except bit k.  For example, a single X error on qubit
	2 is represented by 0b100 (=4). 
	"""
	# First extract a list of all the block IDs used in the computation
	blocks = allBlocks(locations)
	#if len(blocks) != 14:	# used to have only 10 blocks, before adding data snapshots
	#	print "Error!  Too few blocks!"
	
	# now for each location, consider the types of errors that can occur there
	# 	for a 1-qubit location, it is sufficient to consider X and Z errors separately
	# 	for a 2-qubit location, it is sufficient to consider XI, IX, ZI, IZ
	# for each type of error, propagate it forward, and store a dictionary in the location
	def initializeErrors(blocks):
		"""This returns a properly formatted errors array, initialized to 0."""
		errors = {'X':{},'Z':{}}
		for blockID in blocks:
			errors['X'][blockID] = 0
			errors['Z'][blockID] = 0
		return errors
	for k in range(len(locations)):
		loc = locations[k]
		remainingLocations = locations[k+1:]
		
		# AEP: So each location has (at least) 2 'sub' locations X1 and Z1?
		# loc['X1'] represents the total error that results when an X-error occurs at loc (qubit 1) 
		# and is then propagated through the remaining locations in the circuit.  Similarly for loc['Z1'].
		loc['X1'] = initializeErrors(blocks)
		loc['X1']['X'][loc['block1']] = 1<<loc['bit1']
		loc['X1'] = propagateErrors(loc['X1'], remainingLocations)
		loc['Z1'] = initializeErrors(blocks)
		loc['Z1']['Z'][loc['block1']] = 1<<loc['bit1']
		loc['Z1'] = propagateErrors(loc['Z1'], remainingLocations)
		if loc['type'] == 'cnot':
			loc['X2'] = initializeErrors(blocks)
			loc['X2']['X'][loc['block2']] = 1<<loc['bit2']
			loc['X2'] = propagateErrors(loc['X2'], remainingLocations)
			loc['Z2'] = initializeErrors(blocks)
			loc['Z2']['Z'][loc['block2']] = 1<<loc['bit2']
			loc['Z2'] = propagateErrors(loc['Z2'], remainingLocations)

def addErrors(e1, e2):
	"""Adds the X and Z errors stored in e1 and in e2.
	
	This is useful for combining precomputed error consequences.  
	>>> addErrors({'X':{'block0':1,'block1':3},'Z':{'block0':2,'block1':5}}, \
					{'X':{'block0':6,'block1':2},'Z':{'block0':1,'block1':4}})
	{'X': {'block1': 1, 'block0': 7}, 'Z': {'block1': 1, 'block0': 3}}
	"""
	blocks = e1['X'].keys()
	e = {'X':{}, 'Z':{}}
	for key in blocks: 
		e['X'][key] = e1['X'][key] ^ e2['X'][key]
		e['Z'][key] = e1['Z'][key] ^ e2['Z'][key]
	return e

def duplicateErrorBlock(errors, block, newblockname):
	"""Takes a dictionary of errors, and creates a new block with a copy of the X and Z errors in block.
	
	Although it returns the modified error dictionary for convenience, modifications are made in place.  
	>>> duplicateErrorBlock({'X':{'a':0, 'b':1}, 'Z':{'a':2, 'b':3}}, 'b', 'c')
	{'X': {'a': 0, 'c': 1, 'b': 1}, 'Z': {'a': 2, 'c': 3, 'b': 3}}
	"""
	for etype in ('X', 'Z'):
		errors[etype][newblockname] = errors[etype][block]
	return errors


### NOTE: There are built-in itertools functions that do this kind of thing (http://docs.python.org/py3k/library/itertools.html#itertools.combinations_with_replacement)
class SubsetIterator:
	"""Class for iterating over all subsets of a given size of (a shallow copy of) the given list.
	
	The subsets are given in lexicographic order of their indices in range(len(list)).
	>>> [t for t in SubsetIterator(range(1,5),2)]
	[[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]
	"""
	def __init__(self, universe, subsetsize):
		self.universe = universe[:]
		self.subsetsize = subsetsize
	def __iter__(self):
		if self.subsetsize == 0:
			return [[]].__iter__()
		
		self.a = range(self.subsetsize)	# a is used to store the indices of the subset
		self.a[-1] -= 1		# decrement the last index so first increment will return first subset
		return self
	def next(self):
		n = len(self.universe)
		r = self.subsetsize
		a = self.a
		if a[0] == n-r:
			raise StopIteration		
		a[-1] = a[-1] + 1
		if a[-1] > n-1:
			j = r-2
			while a[j] == n-1 - (r-1-j):
				j = j-1
			a[j] = a[j] + 1
			for i in range(j+1, r):
				a[i] = a[j] + i - j
		return [self.universe[x] for x in a]	

class TupleIterator:
	"""Class for iterating over all lists a where a[j] in range(0,ranges[j]).
	
	>>> [t for t in TupleIterator([1,2,3])]
	[[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 1, 0], [0, 1, 1], [0, 1, 2]]
	"""
	def __init__(self, ranges):
		self.ranges = ranges
	def __iter__(self):
		self.a = [0] * len(self.ranges)
		self.a[-1] = -1		# set the last element to -1, so the first call to next returns [0,0,..,0]
		return self
	def next(self):
		a = self.a
		ranges = self.ranges
		r = len(a)
		j = r-1
		while a[j] == ranges[j] - 1:
			j = j - 1
			if j < 0:
				a = [0] * len(ranges)
				raise StopIteration
		a[j] = a[j] + 1
		a[j+1:] = [0] * (r-j-1)
		return a[:]			# we return a copy of the internal list (not necessary for this application)
	
	
class PartitionIterator:
	'''
	Iterates over all possible assignments of of 'whole' objects into 'numParts' parts.
	Returns a list of size numParts during each iteration.  Each element in this
	list specifies the portion of the whole allotted to that element.  The sum of
	all of the elements in the list will always be equal to whole.
	The optional maxInParts argument is a list that specifies the maximum portion
	that can be assigned to each part.  For example, if maxInParts[i]=k, the value of the
	i'th element of the list returned for each iteration will be at most k.
	
	>>> [partition for partition in PartitionIterator(3,3, [3,1,2], [1,0,0])]
	[[3, 0, 0], [2, 1, 0], [2, 0, 1], [1, 1, 1], [1, 0, 2]]
	'''
	
	def __init__(self, whole, numParts, maxInParts=None, minInParts=None):
		self.whole = whole
		self.numParts = numParts
		if maxInParts == None:
			self.maxInParts = [whole] * numParts
		else:
			self.maxInParts = maxInParts
		assert numParts <= len(self.maxInParts)
		
		if minInParts == None:
			self.minInParts = [0] * numParts
		else:
			self.minInParts = minInParts
		assert numParts <= len(self.minInParts)
	
	def __iter__(self):
		if self.numParts == 0:
			return [].__iter__()
		
		self.leadPart = min(self.whole, self.maxInParts[0])
		
		if self.leadPart < self.minInParts[0]:
			return [].__iter__()
		
		if self.numParts == 1:
			if (self.whole == self.leadPart) and (self.whole >= self.minInParts[0]):
				return [[self.leadPart]].__iter__()
			return [].__iter__()
		
		self.subIterator = PartitionIterator(self.whole - self.leadPart, self.numParts - 1, self.maxInParts[1:]).__iter__()		
		return self
	
	def next(self):

		while True:
			
			try:
				return [self.leadPart] + self.subIterator.next()
			except StopIteration:
				self.leadPart -= 1
				
				if self.minInParts[0] > self.leadPart:
					raise StopIteration
				
				self.subIterator = PartitionIterator(self.whole - self.leadPart, self.numParts - 1, self.maxInParts[1:]).__iter__()								
		
class SliceIterator:
	'''
	>>> l = range(7)
	>>> [s for s in SliceIterator(l, 2)]
	[[0, 1], [2, 3], [4, 5], [6]]
	'''
	
	def __init__(self, iterable, sliceLen=1):
		self._iterable = iterable
		self._sliceLen = sliceLen
		
	def __iter__(self):
		self._iterator = self._iterable.__iter__()
		return self
	
	def next(self):

		slice = []
		try:
			for _ in xrange(self._sliceLen):
				slice.append(self._iterator.next())
		except StopIteration:
			pass
		
		if 0 == len(slice):
			raise StopIteration
		
		return slice		


### The following functions define the number of possible errors for every location, and then convert indices into these errors.  
### The suffix says whether it is for X and Z errors, just X errors, or just Z errors.  
def errorRangeXZ(loc):
	if loc['type'] in ['prepX', 'prepZ', 'measX', 'measZ']:
		return 1
	if loc['type'] == 'rest':
		return 3
	if loc['type'] == 'cnot':
		return 15
def locErrorXZ(loc, index):
	if loc['type'] in ['prepX', 'measX'] and index == 0:
		return loc['Z1']
	if loc['type'] in ['prepZ', 'measZ'] and index == 0:
		return loc['X1']
	if loc['type'] == 'rest':
		return reduce(addErrors, [[loc['X1']], [loc['Z1']], [loc['X1'], loc['Z1']]][index])
	if loc['type'] == 'cnot':
		l = [[loc['X1']], [loc['Z1']], [loc['X1'], loc['Z1']]]
		l = l + [a+[loc['X2']] for a in l] + [a+[loc['Z2']] for a in l] + [a+[loc['X2'],loc['Z2']] for a in l]
		l = l + [[loc['X2']], [loc['Z2']], [loc['X2'], loc['Z2']]]
		return reduce(addErrors, l[index])
def errorRangeYZ(loc):
	if loc['type'] in ['prepX', 'prepZ', 'measX', 'measZ']:
		return 1
	if loc['type'] == 'rest':
		return 2
	if loc['type'] == 'cnot':
		return 12
def errorRangeX(loc):	# actually prepX and measX can't cause any X errors, so these locations can be stripped out
	if loc['type'] in ['prepX', 'prepZ', 'measX', 'measZ', 'rest']:
		return 1	
	if loc['type'] == 'cnot':
		return 3
def locErrorX(loc, index):
	if loc['type'] in ['prepX', 'prepZ', 'measX', 'measZ', 'rest']:
		return loc['X1']
	if loc['type'] == 'cnot':
		return reduce(addErrors, [[loc['X1']], [loc['X2']], [loc['X1'], loc['X2']]][index])
def errorRangeXIIX(loc):	# consider only X1 and X2 errors for a CNOT, not XX errors
	if loc['type'] in ['prepX', 'prepZ', 'measX', 'measZ', 'rest']:
		return 1	
	if loc['type'] == 'cnot':
		return 2
def locErrorXIIX(loc, index):
	if loc['type'] in ['prepX', 'prepZ', 'measX', 'measZ', 'rest']:
		return loc['X1']
	if loc['type'] == 'cnot':
		return [loc['X1'], loc['X2']][index]
def errorRangeXI(loc):	# consider only X1 errors for a CNOT, not IX or XX errors
	if loc['type'] in ['prepX', 'prepZ', 'measX', 'measZ', 'rest']:
		return 1	
	if loc['type'] == 'cnot':
		return 1
def locErrorXI(loc, index):
	if loc['type'] in ['prepX', 'prepZ', 'measX', 'measZ', 'rest']:
		return loc['X1']
	if loc['type'] == 'cnot':
		return loc['X1']
def errorRangeZ(loc):	# actually prepZ and measZ can't cause any Z errors, so these locations can be stripped out
	if loc['type'] in ['prepX', 'prepZ', 'measX', 'measZ', 'rest']:
		return 1	
	if loc['type'] == 'cnot':
		return 3
def locErrorZ(loc, index):
	if loc['type'] in ['prepX', 'prepZ', 'measX', 'measZ', 'rest']:
		return loc['Z1']
	if loc['type'] == 'cnot':
		return reduce(addErrors, [[loc['Z1']], [loc['Z2']], [loc['Z1'], loc['Z2']]][index])
	

def convolveCounts(counts1, counts2):
	"""Convolves the input two lists under the exor operation.
	
	That is, returns a list with counts[s] = sum_{s1,s2:s1+s2=s} counts1[s1]*counts2[s2].
	"""
	
	length = max(len(counts1), len(counts2))
	#assert length == len(counts2)
	counts = [0] * length
	
	s1List = nonZeroIndices(counts1)
	s2List = nonZeroIndices(counts2)
	
	logger.info('Convolving {0}x{1}'.format(len(s1List), len(s2List)))
	for s1 in s1List:
		c1 = counts1[s1]  
		for s2 in s2List:
			c2 = counts2[s2]
			s = s1 ^ s2
			counts[s] += c1 * c2
	return counts

def convolvedSum(counts1, counts2):
	'''
	Returns the sum of the convolution of the two counts.
	This is equivalent to sum(convolveCounts(counts1, counts2)), but is
	much faster, in general.
	'''
	s1Sum = sum(counts1)
	s2Sum = sum(counts2)
	return s1Sum * s2Sum

def countLocationTypes(locations):
	'''
	Counts the number of locations according to location type (prepX, prepZ, cnot, ...).
	Returns a dictionary indexed by type.
	'''
	counts = {}
	for l in locations:
		type = l['type']
		counts[type] = counts.get(type, 0) + 1
		
	return counts

def convolveKCounts(kCounts1, kCounts2, kMax=None):
	if None == kMax:
		kMax = len(kCounts1) + len(kCounts2)
		
	kCounts = [0] * (kMax+1)
	
	nonZero1 = nonZeroIndices(kCounts1)
	nonZero2 = nonZeroIndices(kCounts2)
	for k1 in nonZero1:
		c1 = kCounts1[k1]
		for k2 in nonZero2:
			k = k1 + k2
			if kMax >= k:
				kCounts[k1 + k2] += c1 * kCounts2[k2]
			
	return kCounts
	
	
def convolve(counts0, counts1, partMax0=None, partMax1=None, kMax=None, convolveFcn=convolveCounts, extraArgs=[]):
	if None == partMax0:
		partMax0 = len(counts0)-1
	if None == partMax1:
		partMax1 = len(counts1)-1
	if None == kMax:
		kMax = partMax0 + partMax1
		
	convolved = [0] * (kMax+1)
	for k in range(kMax+1):
		convolvedK  = []
		logger.info('Convolving for k={0}'.format(k))
		for k0, k1 in PartitionIterator(k, 2, [partMax0, partMax1]):
			args = [counts0[k0], counts1[k1]] + extraArgs
			convolvedK.append(convolveFcn(*args))
			
		convolved[k] = addLists(*convolvedK)
			
	return convolved


def isFaultTolerant(counts, t, corrector):
	'''
	Checks the given set of counts to make sure that they are fault tolerant.
	i.e. Any error that can occur for k <= t has weight <= k.
	'''
	
	kMax = min(t, len(counts)-1)
	
	logical = (len(counts[0]) == (1<<11))
	success = True
	for k in range(kMax+1):
		sNonZero = nonZeroIndices(counts[k])
		for s in sNonZero:
			e = corrector.getError(s)
			e = corrector.hashError(e, logical)
			if not weight(e) <= k:
				logger.warning('s={0}, e={1}, w(e)={2}, k={3}'.format(s, e, weight(e), k))
				success = False		 
				
	return success	

# to run the doctests, run python or python -v directly on this script
if __name__ == "__main__":
	import doctest
	doctest.testmod()