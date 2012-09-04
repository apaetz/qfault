'''
Utilities for lists and other iterables.
Created on 2011-04-28

@author: adam
'''
import operator
import math

def imax(seq):
	'''
	Returns the index of the largest element in the list.
	>>> imax([1,2,3,5,4])
	3
	'''
	return seq.index(max(seq))

def imin(seq):
	'''
	Returns the index of the smallest element in the list.
	>>> imin([5,4,1,2,3])
	2
	'''
	return seq.index(min(seq))

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

def mul(seq):
	'''
	Returns the product of all elements in the list.
	>>> mul([1,2,3])
	6
	'''
	seq = tuple(seq)
	return reduce(operator.mul, seq[1:], seq[0])

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
	
def cycle(seq):
	pass

def permute(seq, permutation):
	"""Permutes the input sequence, like Mathematica's Permute by with indices starting at 0.
	>>> permute(['a','b','c','d'], [0, 2, 3, 1])
	['a', 'c', 'd', 'b']
	"""
	return [seq[p] for p in permutation]

def chop(seq, lengths):
	'''
	Slice the sequence into smaller sequences of the given lengths.
	>>> chop(['a','b','c','d','e'], [1,2,2])
	[['a'], ['b', 'c'], ['d', 'e']]
	'''
	
	slices = []
	index = 0
	for l in lengths:
		slices.append(seq[index:index+l])
		index += l
	
	return slices

def equal_chop(seq, num_subsequences):
	'''
	Chops the sequence into equal size subsequences.
	The last subsequence may contain fewer elements.
	>>> equal_chop(range(10),3)
	[[0, 1, 2, 3], [4, 5, 6, 7], [8, 9]]
	'''
	step = long(math.ceil(len(seq) / float(num_subsequences)))
	chopped = [seq[i * step : i * step + step]
			   for i in range(num_subsequences)]
	return chopped

def remove_subsequence(seq, remove_indices):
	'''
	Returns a sequence in which the given subsequence has been removed.
	>>> remove_subsequence(['a', 'b', 'c', 'd'], [1,3])
	['a', 'c']
	'''
	return [seq[i] for i in range(len(seq)) if i not in set(remove_indices)]
	


if __name__ == "__main__":
	import doctest
	doctest.testmod()