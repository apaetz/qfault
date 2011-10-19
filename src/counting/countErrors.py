'''

This file contains functions counting error propagation within a given set of locations. 

@author: adam
'''


from collections import namedtuple
from counting.key import SyndromeKeyGenerator, DefaultErrorKeyGenerator, \
	MultiBlockSyndromeKeyGenerator
from qec.error import Pauli, PauliError, xType, zType
from qec.qecc import StabilizerCode
from util import iteration
from util.cache import fetchable
from util.listutils import nonZeroIndices
import logging
import math
import operator
import qec.error as error
import util.counterUtils


logging.basicConfig(level=logging.INFO,
					format='[%(asctime)s] %(process)d %(levelname)-8s %(name)s: %(message)s')

logger = logging.getLogger('Counting')
	
CountResult = namedtuple('CountResult', ['counts', 'countsRejected', 'locTotals', 'prAccept', 'prBad'])


def countErrorsParallel(k, locations, lsCountingFcn, extraArgs=[]):

	counts = {}
	
	if 0 == k:
		# Only possible syndrome with 0 failures is the trivial one
		counts[0] = 1
		return counts
		
	if 0 == len(locations):
		# No locations (and at least 1 failure), so all counts must be zero.
		return counts
		
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


def countBlocksBySyndrome(locations, blocks, pauli, noise, kMax):
	filtered = filterAndPropagate(locations, pauli)
	
	keyGenerator = MultiBlockSyndromeKeyGenerator(blocks)
	counts = [countErrors(k, filtered, noise, keyGenerator) for k in range(kMax+1)]
	
	return counts, keyGenerator


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
	
	extraArgs = [locations.blocknames(), noise, keyGenerator]
	return countErrorsParallel(k, locations, countLocationSets, extraArgs)


def countLocationSets(lSets, blocknames, noise, keyGenerator):

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
			counts[errorKey] = counts.get(errorKey, 0) + totalWeight

	return counts


def convolveCountsPostselect(counts1, counts2, syndromeBits):
	'''
	Convolve counts from the input ancilla (counts1) and the error
	detection and postselection network (counts2).  The resulting
	counts are only those which pass error detection.
	'''
	counts = [0] * (1<<syndromeBits)
	mask = (1<<syndromeBits)-1
	s2List = nonZeroIndices(counts2)
	
	logger.info('Postselect counting {0}'.format(len(s2List)))
	
	for s2 in s2List:
		s2a = s2 >> syndromeBits
		s2b = s2 & mask
		counts[s2a ^ s2b] += counts1[s2b] * counts2[s2]
	return counts


def mapKeys(counts, keymap):
	mapped = []
	for kcounts in counts:
		mapped.append({keymap(e): count for e,count in kcounts.items()})
		
	return mapped

def convolveDict(counts1, counts2):
	'''
	Convolve counts from two dictionaries.
	'''
	counts = {}
	logger.info('Convolving {0}x{1}'.format(len(counts2), len(counts1)))
	for key2, count2 in counts2.iteritems():
		for key1, count1 in counts1.iteritems():
			key = key1 ^ key2
			
			#TODO: not sure if it would be faster to use a try-except block here, instead.
			counts[key] = counts.get(key, 0) + (count2 * count1)
	
	return counts


###################
# TODO: most of the convolution functions have the number of syndrome bits hardcoded.
# Need to parameterize.
###################

def convolveXZRejected(countsA1, countsVX):
	
	mask11 = (1<<11) - 1
	mask12 = (1<<12) - 1
	counts = [0] * (1<<11)
	logger.info('Convolving {0}x{1}'.format(len(countsA1), len(countsVX)))
	for sVX, countVX in countsVX.iteritems():
		sXb = (sVX >> 11) & mask12
		sZa = (sVX >> 23) & mask11
		for sA1, countA1 in countsA1.iteritems():
			sX = sA1 >> 11
			if sX != sXb:
				sZ = (sA1 & mask11) ^ sZa
				counts[sZ] += countA1 * countVX

	return counts



def convolveXZPostselectX(countsA1, countsVX):
	
	mask12 = (1<<12) - 1
	counts = {}
	logger.info('Convolving {0}x{1}'.format(len(countsA1), len(countsVX)))
	for sVX, countVX in countsVX.iteritems():
		sXb = (sVX >> 11) & mask12
		sAVX = sVX >> 23
		for sA1, countA1 in countsA1.iteritems():
			sX = sA1 >> 11
			if sX == sXb:
				sA = sA1 ^ sAVX
				counts[sA] = counts.get(sA,0) + countA1 * countVX

	return counts

def convolveXZRejectedX(countsA1, countsVX):
	
	mask12 = (1<<12) - 1
	counts = {}
	logger.info('Convolving {0}x{1}'.format(len(countsA1), len(countsVX)))
	for sVX, countVX in countsVX.iteritems():
		sXb = (sVX >> 11) & mask12
		sAVX = sVX >> 23
		for sA1, countA1 in countsA1.iteritems():
			sX = sA1 >> 11
			if sX != sXb:
				sA = sA1 ^ sAVX
				counts[sA] = counts.get(sA,0) + countA1 * countVX

	return counts


def convolveXZPostselectZ(countsA1, countsVZ):
	
	mask11 = (1<<11) - 1
	counts = {}
	logger.info('Convolving {0}x{1}'.format(len(countsA1), len(countsVZ)))
	for sVZ, countVZ in countsVZ.iteritems():
		sZb = sVZ & mask11
		sAVZ = sVZ >> 23
		for sA1, countA1 in countsA1.iteritems():
			sZ = sA1 & mask11
			if sZ == sZb:
				sA = sA1 ^ sAVZ
				counts[sA] = counts.get(sA,0) + countA1 * countVZ

	return counts

def convolveXZRejectedZ(countsA1, countsVZ):
	
	mask11 = (1<<11) - 1
	counts = {}
	logger.info('Convolving {0}x{1}'.format(len(countsA1), len(countsVZ)))
	for sVZ, countVZ in countsVZ.iteritems():
		sZb = sVZ & mask11
		sAVZ = sVZ >> 23
		for sA1, countA1 in countsA1.iteritems():
			sZ = sA1 & mask11
			if sZ != sZb:
				sA = sA1 ^ sAVZ
				counts[sA] = counts.get(sA,0) + countA1 * countVZ

	return counts

def convolveCountsXZPostselectX(counts1, counts2):
	'''
	Convolve counts (of X and Z errors) from the input ancilla (counts1) and the error
	detection and postselection network (counts2).  The resulting
	counts are only those which pass X-error detection.
	'''
	
	print 'Convolving {0}x{1}'.format(len(counts1), len(counts2))
	counts = [0] * (1<<23)
	mask12 = ((1<<12)-1)
	mask11 = ((1<<11)-1)
	for s2 in counts2.keys():
		# Syndrome s2 is bit string composed of (msb first): block1 X, block1 Z, block2 X, block2 Z syndromes
		s2Xa = s2 >> 34
		s2Za = (s2 >> 23) & mask11
		s2Xb = (s2 >> 11) & mask12
		sXOut = (s2Xa ^ s2Xb) << 11
		count2 = counts2[s2]
		s1List = [s for s in counts1.keys() if s2Xb == (s >> 11)]
		for s1 in s1List:
			s1Z = s1 & mask11
			sZOut = s1Z ^ s2Za
			sOut = sXOut + sZOut
			counts[sOut] += counts1[s1] * count2
						
	return counts


def convolveABB(countsAB, countsB):
	'''
	Convolves counts over two blocks A, B with counts on a single block B
	'''
	counts = [0] * len(countsAB)
	sABList = nonZeroIndices(countsAB)
	sBList = nonZeroIndices(countsB)
	logger.info('Convolving ABB {0}x{1}'.format(len(sABList), len(sBList)))
	for sB in sBList:
		countB = countsB[sB]
		for sAB in sABList:
			counts[sAB ^ sB] += countB * countsAB[sAB]
						
	return counts




def convolveABA(countsAB, countsA, bitsA=12, bitsB=12):
	'''
	Convolves counts over two blocks A, B with counts on a single block A
	'''
	counts = [0] * (1<<(bitsA + bitsB))
	sABList = nonZeroIndices(countsAB)
	sAList = nonZeroIndices(countsA)
	logger.info('Convolving ABA {0}x{1}'.format(len(sABList), len(sAList)))
	for sA in sAList:
		countA = countsA[sA]
		sAShift = sA << bitsB
		for sAB in sABList:
			counts[sAB ^ sAShift] += countA * countsAB[sAB]
	return counts


#def reduceAllErrorSyndromes(locations, code, eigenstate=None):
#	"""Reduces errors to error syndromes for locations on an encoded block."""
#	for l in locations:
#		for name in l['X1']['X'].keys(): 
#			reduceError(l, name, code, eigenstate)	
#				
#
#def reduceAllErrorSyndromesZero(locations, code):
#	return reduceAllErrorSyndromes(locations, code, Pauli.Z)
#
#def reduceAllErrorSyndromesPlus(locations, code):
#	return reduceAllErrorSyndromes(locations, code, Pauli.X)
#		
#def reduceError(l, name, code, eigenstate):
#	for pauli in [str(Pauli.X), str(Pauli.Z)]:
#		# All locations have a node on block 1.  Only some have a node
#		# on block 2.
#		l[pauli+'1'][pauli][name] = code.getSyndrome(l[pauli+'1'][pauli][name])
#		if pauli+'2' in l:
#			l[pauli+'2'][pauli][name] = code.getSyndrome(l[pauli+'2'][pauli][name])


def filterAndPropagate(locations, pauli):
	# |+> preparations and X-basis measurements cannot cause X errors
	# (and similarly for Z). So they need not be counted.
	if pauli != Pauli.Y:
		locations = locations.filterAgainst('meas' + str(pauli))
		locations = locations.filterAgainst('prep' + str(pauli))
	util.counterUtils.propagateAllErrors(locations)
		
	return locations



