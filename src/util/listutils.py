'''
Utilities for lists and other iterables.
Created on 2011-04-28

@author: adam
'''
import operator

def imax(list):
	'''
	Returns the index of the largest element in the list.
	'''
	maxIndex = 0
	max = list[0]
	for i, val in enumerate(list):
		if val > max:
			maxIndex = i
			max = val
			
	return maxIndex

def imin(list):
	'''
	Returns the index of the smallest element in the list.
	'''
	minIndex = 0
	min = list[0]
	for i, val in enumerate(list):
		if val < max:
			minIndex = i
			min = val
			
	return minIndex

def addLists(*lists):
	'''
	Returns the entry-wise sum of the given lists.  All of the lists
	must be the same length.
	'''
#	return [sum(l[x] for l in lists) for x in range(len(lists[0]))]
#	return [sum(x) for x in zip(*lists)]
	return [sum(l[x] for l in lists) for x in xrange(len(lists[0]))]

def mulLists(*lists):
	'''
	Returns the entry-wise product of the given lists.
	'''
	return [reduce(operator.mul, (l[i] for l in lists), 1) for i in xrange(len(lists[0]))]

def mul(list):
	'''
	Returns the product of all elements in the list.
	'''
	return reduce(operator.mul, list, 1)

def addDicts(*dicts):
	'''
	Returns the key-wise sum of the given dictionaries.
	'''
	
	result = {}
	for dict in dicts:
		for key, val in dict.iteritems():
			result[key] = result.get(key, 0) + val
		
	return result


def numNonzeroEntries(li):
	"""Returns the number of nonzero entries in the input list."""
	numEntries = len(li)
	for el in li:
		if el == 0:
			numEntries -= 1
	return numEntries

def nonZeroIndices(list):
	'''
	Returns a list indicies for which list[index] != 0
	'''
	return [s for s in xrange(len(list)) if list[s]]

def uniqify(seq):
	'''
	Returns a list identical to seq, except with all duplicate
	elements removed.  Relative ordering is preserved.
	
	Code copied from http://www.peterbe.com/plog/uniqifiers-benchmark.
	
	>>> uniqify([1,1,3,2,3])
	[1, 3, 2]
	'''

	seen = {} 
	result = [] 
	for item in seq: 
		if item in seen: continue 
		seen[item] = 1 
		result.append(item) 
	return result

if __name__ == "__main__":
	import doctest
	doctest.testmod()
