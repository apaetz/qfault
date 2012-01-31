'''

Functions counting error propagation within a given set of (physical) locations. 

@author: adam
'''


from collections import namedtuple
from counting.key import DefaultErrorKeyGenerator, \
	MultiBlockSyndromeKeyGenerator
from qec.error import Pauli, PauliError, xType, zType
from util import iteration, counterUtils
import logging
import qec.error as error


logging.basicConfig(level=logging.INFO,
					format='[%(asctime)s] %(process)d %(levelname)-8s %(name)s: %(message)s')

logger = logging.getLogger('Counting')
	
CountResult = namedtuple('CountResult', ['counts', 'countsRejected', 'locTotals', 'prAccept', 'prBad'])


def countErrorsParallel(k, locations, lsCountingFcn, extraArgs=[]):

	counts = {}
	
#	if 0 == len(locations):
#		# No locations (and at least 1 failure), so all counts must be zero.
#		return counts
		
	import countParallel
	pool = countParallel.getPool()
	
	lSets = [ls for ls in iteration.SubsetIterator(locations, k)]
	lSetSlices = countParallel.packTasks(100, lSets, [1] * len(lSets))
	lsResults = []
	
	for lSetSlice in lSetSlices:
		lsResults.append(pool.apply_async(lsCountingFcn, [lSetSlice] + extraArgs))
		
	for result in lsResults:
		lsCounts = result.get()
		for e, count in lsCounts.iteritems():
			counts[e] = counts.get(e,0) + count
			
	return counts


def countBlocksBySyndrome(locations, blocks, noise, kMax):
	counterUtils.propagateAllErrors(locations)
	
	keyGenerator = MultiBlockSyndromeKeyGenerator(blocks)
	counts = [countErrors(k, locations, noise, keyGenerator) for k in range(kMax+1)]
	
	return counts

def mapCounts(counts, keymap):

	newCounts = []
	for countsK in counts:
		newCountsK = {}
		for key,count in countsK.iteritems():
			mappedKey = keymap(key)
			newCountsK[mappedKey] = newCountsK.get(mappedKey, 0) + count
			
		newCounts.append(newCountsK)
		
	return newCounts



#@fetchable
def countErrors(k, locations, noise, keyGenerator=DefaultErrorKeyGenerator):
	""" Identical to countErrors1BlockYZ(), except that the locations may now span two blocks.
	When considering both Y and Z errors, there can be (for encoded |0>) as many as 2^44
	syndromes across two blocks.  In order to manage this, counts are kept in dictionaries
	indexed by syndrome.  Only syndromes with non-zero counts are entered into the dictionary.
	
	Returned counts are indexed by syndromes of the following form.  The first zBits2 are the
	Z syndrome of the second block, the next xBits2 are the X syndrome of the second block,
	the next zBits1 are the Z syndrome of the first block, and the final xBits1 are the X
	syndrome of the first block.
	
	TODO: This documentation is old. Update documentation.
	
	"""
	
	blocknames = locations.blocknames()
	
	if 0 == k:
		# Only possible error with 0 failures is the trivial one
		error = {name: Pauli.I for name in blocknames}
		return {keyGenerator.getKey(error): 1}
		
	extraArgs = [keyGenerator, blocknames, noise]
	return countErrorsParallel(k, locations, countLocationSets, extraArgs)


def countLocationSets(lSets, keyGenerator, blocknames, noise):

	counts = {}
	
	for ls in lSets:
		# The set of possible errors can be different for each location.
		errorLookup = map(noise.errorList, ls)
		
		k = len(ls)
		for eIndexes in iteration.TupleIterator(map(len, errorLookup)):
			# ls is the set of locations under consideration.
			# eIndexes is a list of indexes, one for each location, that represent error types.
			
			X = [0 for _ in blocknames]
			Z = [0 for _ in blocknames]
			totalWeight = 1
			
			for j in range(k):
				l = ls[j]				
				e = errorLookup[j][eIndexes[j]]
				
				totalWeight *= noise.getWeight(l,e)
				#print 'l=', l
				#print 'k=', k, 'j=', j, 'eIndex=', eIndexes[j]
				#print 'e=', e, 'weight=', totalWeight
				
				# X-error at location l on qubit 1 (i.e., IX, or IY, or XX, etc.)
				if e[xType] & 1:
					X = [X[i] ^ l['X1']['X'][name] for i,name in enumerate(blocknames)]
					
				# X-error at location l on qubit 2 (i.e., XI, or YI, or XX, etc.)
				if e[xType] & 2:
					X = [X[i] ^ l['X2']['X'][name] for i,name in enumerate(blocknames)]
					
				# (IZ)
				if e[zType] & 1:
					Z = [Z[i] ^ l['Z1']['Z'][name] for i,name in enumerate(blocknames)]					
								
				# (ZI)
				if e[zType] & 2:
					Z = [Z[i] ^ l['Z2']['Z'][name] for i,name in enumerate(blocknames)]

			blockErrors = {name: PauliError(X[i], Z[i]) for i,name in enumerate(blocknames)}
			errorKey = keyGenerator.getKey(blockErrors)
			#print 'blockErrors=', blockErrors, 'key=', errorKey
			counts[errorKey] = counts.get(errorKey, 0) + totalWeight

	return counts

def pauliFilter(locations, pauli):
	# TODO: this function doesn't belong here.  Where should it go?
	# |+> preparations and X-basis measurements cannot cause X errors
	# (and similarly for Z). So they need not be counted.
	if pauli != Pauli.Y:
		locations = locations.filterAgainst('meas' + str(pauli))
		locations = locations.filterAgainst('prep' + str(pauli))
		
	return locations



