'''

This file contains functions counting error propagation within a given set of locations. 

@author: adam
'''


from util.cache import fetchable
from util.listutils import nonZeroIndices
from collections import namedtuple
import countParallel
import util.counterUtils
from golay import golayCode
import copy
import logging


logging.basicConfig(level=logging.INFO,
					format='[%(asctime)s] %(process)d %(levelname)-8s %(name)s: %(message)s')

logger = logging.getLogger('Counting')

corrector = golayCode.Corrector()	# initialize the ideal decoder
	
CountResult = namedtuple('CountResult', ['counts', 'countsRejected', 'locTotals', 'prAccept', 'prBad'])


zeroBitsForType = {'X': 12, 'Z': 11}
		
def countXerrorsZero(k, locations, name, noise):
	'''
	Takes a set of locations on a single Golay encoded |0> block and returns a list indexed by syndrome with the
	(weighted) number of X errors leading to each syndrome when exactly k locations fail.
	'''
	return countErrors1Block(k, locations, name, 12, 'X', noise['X'])

def countErrorsZero(k, locations, name, noise, type):
	return countErrors1Block(k, locations, name, zeroBitsForType[type], type, noise[type])

def countXErrorsGolay(k, locations, name, noise):
	'''
	Takes a set of locations on a single Golay encoded block and returns a list indexed by syndrome with the
	(weighted) number of X errors leading to each syndrome when exactly k locations fail.
	'''
	return countXerrorsZero(k, locations, name, noise)

def countZerrorsZero(k, locations, name, noise):
	'''
	Takes a set of locations on a single Golay encoded |0> block and returns a list indexed by syndrome with the
	(weighted) number of Z errors leading to each syndrome when exactly k locations fail.
	'''
	return countErrors1Block(k, locations, name, 11, 'Z', noise['Z'])

def countZerrorsPlus(k, locations, name, noise):
	'''
	Takes a set of locations on a single Golay encoded |+> block and returns a list indexed by syndrome with the
	(weighted) number of Z errors leading to each syndrome when exactly k locations fail.
	'''
	return countErrors1Block(k, locations, name, 12, 'Z', noise['Z'])

def countZerrorsGolay(k, locations, name, noise):
	'''
	Takes a set of locations on a single Golay encoded block and returns a list indexed by syndrome with the
	(weighted) number of Z errors leading to each syndrome when exactly k locations fail.
	'''
	return countZerrorsPlus(k, locations, name, noise)

def countXZErrorsZero(k, locations, name, noise):
	'''
	'''
	return countErrors1BlockXZ(k, locations, name, 11, 12, noise['XZ'])	

def countXerrorsZeroZero(k, locations, name1, name2, noise):
	'''
	Takes a set of locations on two encoded |0> blocks and returns a list indexed by syndrome with the
	(weighted) number of X errors leading to each syndrome when exactly k locations fail.
	'''
	return countErrors2blocks(k, locations, name1, 12, name2, 12, 'X', noise['X'])

def countErrorsZeroZero(k, locations, name1, name2, noise, type):
	'''
	Takes a set of locations on two encoded |0> blocks and returns a list indexed by syndrome with the
	(weighted) number of X errors leading to each syndrome when exactly k locations fail.
	'''
	bits = zeroBitsForType[type]
	return countErrors2blocks(k, locations, name1, bits, name2, bits, type, noise[type])

def countZerrorsZeroZero(k, locations, name1, name2, noise):
	'''
	Takes a set of locations on two encoded |0> blocks and returns a list indexed by syndrome with the
	(weighted) number of Z errors leading to each syndrome when exactly k locations fail.
	'''
	return countErrors2blocks(k, locations, name1, 11, name2, 11, 'Z', noise['Z'])

def countXZErrorsZeroZero(k, locations, name1, name2, noise):
	return countErrors2blocksXZ(k, locations, name1, 11, 12, name2, 11, 12, noise['XZ'])


@fetchable
def countErrors1Block(k, locations, blockname, errorBits, type, noise):
	"""This takes a set of locations on one block, and returns a list indexed by syndrome, 
	with the number of X (or Z) error possibilities leading to each syndrome when exactly k locations fail.
	A single gate can fail as XI, IX or XX (or ZI, IZ, ZZ).  For a logical |0>, the X error syndrome is 12 bits, the Z syndrome is 11 bits.
	
	When considering subsets of different kinds of locations, reweighting is necessary.  
	This is because, working with likelihoods, the CNOT parameter is 4g/(1-12g), where g=p/15 and p is the total error 
	rate (e.g., p=1/1000); the rest likelihood parameter is 8g/(1-8g); and the prep/meas parameter is 4g/(1-4g).  
	"""
	counts = [0] * (1<<errorBits)
	
	dictCounts = countErrorsParallel(k, locations, countLocationSets1Block, [type, blockname, noise])
	
	# Convert from a dictionary to a list.
	for s, count in dictCounts.iteritems():
		counts[s] = count
			
	return counts

def countErrorsParallel(k, locations, lsCountingFcn, extraArgs=[]):

	counts = {}
	
	if 0 == k:
		# Only possible syndrome with 0 failures is the trivial one
		counts[0] = 1
		return counts
		
	if 0 == len(locations):
		# No locations (and at least 1 failure), so all counts must be zero.
		return counts
		
	lSets = [ls for ls in util.counterUtils.SubsetIterator(locations, k)]
	lSetSlices = countParallel.packTasks(100, lSets, [1] * len(lSets))
	lsResults = []
	pool = countParallel.getPool()
	for lSetSlice in lSetSlices:
		lsResults.append(pool.apply_async(lsCountingFcn, [lSetSlice] + extraArgs))
		
	for result in lsResults:
		lsCounts = result.get()
		for e, count in lsCounts.iteritems():
			counts[e] = counts.get(e,0) + count
			
	return counts

def countLocationSets1Block(lSets, type, blockname, noise):
	type1 = type + '1'
	type2 = type + '2'
	counts = {}
	for ls in lSets:
		k = len(ls)
		for es in util.counterUtils.TupleIterator(map(util.counterUtils.errorRangeX, ls)):
			# ls is the set of locations under consideration, es the error types: XI, IX or XX
			totalError = 0
			totalWeight = 1
			for j in range(k):
				l = ls[j]
				e = es[j] + 1	# after adding 1, e is now in {1,2,3}, corresponding to XI, IX and XX
				totalWeight *= noise.getWeight(l,e-1)
				if e&1:
					totalError ^= l[type1][type][blockname]
				if e&2:
					totalError ^= l[type2][type][blockname]
			counts[totalError] = counts.get(totalError,0) + totalWeight

	return counts


@fetchable
def countErrors1BlockXZ(k, locations, blockname, zBits, xBits, noise):
	"""
	Counts Y and Z errors for a single (23-qubit) encoded block.  X errors (i.e. errors that
	contain no Z Pauli) are not counted.
	
	Returned counts are a dictionary indexed by a 23-bit syndrome.  The least significant zBits bits are
	the Z syndrome, the remaining xBits bits are the X syndrome.  Syndromes for which the count is 0 are
	not included in the dictionary.
	"""
	return countErrorsParallel(k, locations, countLocationSets1BlockXZ, [blockname, zBits, noise])



def countLocationSets1BlockXZ(lSets, blockname, zBits, noise):
	
	eXI = 1
	eZI = 4
	errorSets = {
				 'prepX': [eZI],
				 'prepZ': [eXI],
				 'measX': [eZI],
				 'measZ': [eXI],
				 'rest': [eXI, eZI, eXI + eZI],
				 'cnot': range(16)
				 }

	counts = {}
	for ls in lSets:
		k = len(ls)
		for es in util.counterUtils.TupleIterator(map(util.counterUtils.errorRangeXZ, ls)):
			# ls is the set of locations under consideration, es the error types: XI, IX, ...
			totalError = 0
			totalWeight = 1
			for j in range(k):
				l = ls[j]
				
				# after adding 1, e is now in {1,2,...,15}, corresponding to XI, IX, ..., YY
				# i.e. bit0 ~ XI, bit1 ~ IX, bit2 ~ ZI, bit3 ~ IZ
				eIndex = es[j]
				e = errorSets[l['type']][eIndex]
				
				totalWeight *= noise.getWeight(l, e-1)
				
				if e&1:
					totalError ^= (l['X1']['X'][blockname] << zBits)
				if e&2:
					totalError ^= (l['X2']['X'][blockname] << zBits)
				if e&4:
					totalError ^= l['Z1']['Z'][blockname]
				if e&8:
					totalError ^= l['Z2']['Z'][blockname]
					
			counts[totalError] = counts.get(totalError, 0) + totalWeight

	return counts




@fetchable
def countErrors2blocks(k, locations, name1, errorBits1, name2, errorBits2, type, noise):
	""" Identical to countErrors1Block(), except that the locations may now span two blocks."""
	counts = [0] * (1 << (errorBits1 + errorBits2))
	
	extraArgs = [type, name1, errorBits1, name2, errorBits2, noise]
	
	dictCounts = countErrorsParallel(k, locations, countLocationSets2Blocks, extraArgs)
	
	for s, count in dictCounts.iteritems():
		counts[s] = count
	
	return counts


def countLocationSets2Blocks(lSets, type, name1, errorBits1, name2, errorBits2, noise):
	type1 = type + '1'
	type2 = type + '2'
	
	counts = {}
	for ls in lSets:
		k = len(ls)
		for es in util.counterUtils.TupleIterator(map(util.counterUtils.errorRangeX, ls)):
			# ls is the set of CNOT locations under consideration, es the error types: XI, IX or XX
			totalError = 0
			totalWeight = 1
			for j in range(k):
				l = ls[j]
				e = es[j] + 1	# after adding 1, e is now in {1,2,3}, corresponding to XI, IX and XX
				totalWeight *= noise.getWeight(l,e-1)
				if e&1:
					totalError ^= (l[type1][type][name1]<<errorBits2) + l[type1][type][name2]
				if e&2:
					totalError ^= (l[type2][type][name1]<<errorBits2) + l[type2][type][name2]
			counts[totalError] = counts.get(totalError,0) + totalWeight
	
	return counts

@fetchable
def countErrors2blocksXZ(k, locations, name1, zBits1, xBits1, name2, zBits2, xBits2, noise):
	""" Identical to countErrors1BlockYZ(), except that the locations may now span two blocks.
	When considering both Y and Z errors, there can be (for encoded |0>) as many as 2^44
	syndromes across two blocks.  In order to manage this, counts are kept in dictionaries
	indexed by syndrome.  Only syndromes with non-zero counts are entered into the dictionary.
	
	Returned counts are indexed by syndromes of the following form.  The first zBits2 are the
	Z syndrome of the second block, the next xBits2 are the X syndrome of the second block,
	the next zBits1 are the Z syndrome of the first block, and the final xBits1 are the X
	syndrome of the first block.
	
	TODO: It might make more sense to arrange the syndrome by X and Z instead of by block.
	
	"""
	extraArgs = [name1, xBits1, zBits1, name2, xBits2, zBits2, noise]
	return countErrorsParallel(k, locations, countLocationSets2BlocksXZ, extraArgs)


def countLocationSets2BlocksXZ(lSets, name1, xBits1, zBits1, name2, xBits2, zBits2, noise):
		
	block2Bits = zBits2 + xBits2
	block1XShift = zBits1 + block2Bits
	
	eXI = 1
	eZI = 4
	errorSets = {
				 'prepX': [eZI],
				 'prepZ': [eXI],
				 'measX': [eZI],
				 'measZ': [eXI],
				 'rest': [eXI, eZI, eXI + eZI],
				 'cnot': range(16)
				 }
	
	counts = {}
	
	for ls in lSets:
		k = len(ls)
		for es in util.counterUtils.TupleIterator(map(util.counterUtils.errorRangeXZ, ls)):
			# ls is the set of CNOT locations under consideration, es the error types: XI, IX, ...
			totalError = 0
			totalWeight = 1
			for j in range(k):
				l = ls[j]
				# after adding 1, e is now in {1,2,...,15}, corresponding to XI, IX, ..., YY
				# i.e. bit0 ~ XI, bit1 ~ IX, bit2 ~ ZI, bit3 ~ IZ
				eIndex = es[j]
				e = errorSets[l['type']][eIndex]
				totalWeight *= noise.getWeight(l,e-1)
				
				if e&1:
					totalError ^= (l['X1']['X'][name1] << block1XShift) + \
								  (l['X1']['X'][name2] << zBits2)
				if e&2:
					totalError ^= (l['X2']['X'][name1] << block1XShift) + \
								  (l['X2']['X'][name2] << zBits2)					 
				if e&4:
					totalError ^= (l['Z1']['Z'][name1]<< block2Bits) + \
								  l['Z1']['Z'][name2]
				if e&8:
					totalError ^= (l['Z2']['Z'][name1] << block2Bits) + \
								  l['Z2']['Z'][name2]

			counts[totalError] = counts.get(totalError, 0) + totalWeight

	return counts

def convolveCountsPostselectX(counts1, counts2):
	return convolveCountsPostselect(counts1, counts2, 12)

def convolveCountsPostselectZ(counts1, counts2):
	return convolveCountsPostselect(counts1, counts2, 11)

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


def convolveDict(counts1, counts2):
	'''
	Convolve counts from two dictionaries.
	'''
	counts = {}
	logger.info('Convolving {0}x{1}'.format(len(counts2), len(counts1)))
	for s2, count2 in counts2.iteritems():
		for s1, count1 in counts1.iteritems():
			s = s1 ^ s2
			counts[s] = counts.get(s, 0) + (count2 * count1)
	
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


def reduceAllErrorSyndromes(locations):
	"""Reduces errors to error syndromes for Golay encoded locations."""
	for l in locations:
		for name in l['X1']['X'].keys(): 
			reduceGolay(l, name)	
				

def reduceAllErrorSyndromesZero(locations):
	"""Reduces errors to error syndromes for encoded |0>."""
	for l in locations:
		for name in l['X1']['X'].keys(): 
			reduceZero(l, name)

def reduceAllErrorSyndromesPlus(locations):
	"""Reduces errors to error syndromes for encoded |+>."""
	for l in locations:
		for name in l['X1']['X'].keys(): 
			reducePlus(l, name)
			


def reducePlus(l, name):
	l['X1']['X'][name] = corrector.getSyndrome(l['X1']['X'][name])
	l['Z1']['Z'][name] = corrector.getLogicalSyndrome(l['Z1']['Z'][name])
	if 'X2' in l:
		l['X2']['X'][name] = corrector.getSyndrome(l['X2']['X'][name])
	if 'Z2' in l:
		l['Z2']['Z'][name] = corrector.getLogicalSyndrome(l['Z2']['Z'][name])
		
def reduceZero(l, name):
	l['X1']['X'][name] = corrector.getLogicalSyndrome(l['X1']['X'][name])
	l['Z1']['Z'][name] = corrector.getSyndrome(l['Z1']['Z'][name])
	if 'X2' in l:
		l['X2']['X'][name] = corrector.getLogicalSyndrome(l['X2']['X'][name])
	if 'Z2' in l:
		l['Z2']['Z'][name] = corrector.getSyndrome(l['Z2']['Z'][name])
		
def reduceGolay(l, name):
	l['X1']['X'][name] = corrector.getLogicalSyndrome(l['X1']['X'][name])
	l['Z1']['Z'][name] = corrector.getLogicalSyndrome(l['Z1']['Z'][name])
	if 'X2' in l:
		l['X2']['X'][name] = corrector.getLogicalSyndrome(l['X2']['X'][name])
	if 'Z2' in l:
		l['Z2']['Z'][name] = corrector.getLogicalSyndrome(l['Z2']['Z'][name])


#
# X-error verification
#

def propagateAndReduceZeroX(locations):
	return propagateAndReduceZero(locations, 'X')

def propagateAndReduceZeroZ(locations):
	return propagateAndReduceZero(locations, 'Z')

def propagateAndReduceZeroXZ(locations):
	locations = copy.copy(locations)
	util.counterUtils.propagateAllErrors(locations)
	reduceAllErrorSyndromesZero(locations)
	
	return locations


def propagateAndReduceZero(locations, type=None):
	assert type == 'X' or type == 'Z' or type == None
	# |+> preparations and X-basis measurements cannot cause X errors
	# (and similarly for Z). So they need not be counted.
	if type != None:
		locations = locations.filterAgainst('meas' + type)
		locations = locations.filterAgainst('prep' + type)
	util.counterUtils.propagateAllErrors(locations)
	reduceAllErrorSyndromesZero(locations)
	
	return locations

def propagateReduceAndCountZero(locations, type, blockname, noise, kMax):
	locations = propagateAndReduceZero(locations, type)
	counts = [countErrorsZero(k, locations, blockname, noise, type) for k in range(kMax + 1)]
	return locations, counts



