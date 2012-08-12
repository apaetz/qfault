'''
 ec.py

 This file contains functions for counting X and Z errors produced at the
 output of Steane error correction.

 Steane error correction works by measuring the X-error and Z-error syndromes of the
 encoded data block separately.  First, Z errors are copied from the data block to
 an encoded |0> ancilla by a transversal CNOT from the ancilla to the data.  Then,
 X errors are copied onto an encoded |+> ancilla by a transveral CNOT from the data
 to the ancilla.  The error syndrome is obtained by transversally measuring the |0>
 ancilla in the X-basis, and transversally measuring the |+> ancilla in the Z-basis
 and then classically processing the result.
 
 There are up to four error corrections in a single exRec.  However, for the
 purposes of exRec error counting, the (up to) two leading error corrections 
 are equivalent, as are the two trailing error corrections.
 
 There are two functions indended for external use:
 countLEC_xOnly
 countTEC_xOnly
 countLEC_zOnly
 countTEC_zOnly
'''

from counting.probability import calcPrBadECZ, calcPrBadECX
from component.ancilla import countVerifiedZeroXOnly, countVerifiedZeroZOnly, \
	countVerifiedPlusZOnly, countVerifiedPlusXOnly
from counting import countParallel
from counting.countErrors import countXErrorsGolay, CountResult, \
	countErrors2blocks, convolveABB, countZerrorsGolay, convolveABA, reduceZero, \
	reduceGolay, reducePlus
from counting.countParallel import convolve
from counting.location import Locations
from util.cache import fetchable
from util.listutils import nonZeroIndices
import util.counterUtils
from golay import golayCode
import logging
import time

logger = logging.getLogger('component.ec')
corrector = golayCode.Corrector()

def locationsRCM_ECZ(type):
	# Error Correction includes the transversal CNOT, transversal X-basis measurement,
	# and rests due to the Z-error verification measurement.  
	locationsEC = [util.counterUtils.locrest('B', i) for i in range(23)] + \
				  [util.counterUtils.loccnot('B', i, 'A', i) for i in range(23)] + \
				  [util.counterUtils.locXmeas('B', i) for i in range(23)]
	locationsEC = Locations(locationsEC, 'rcm.ecZ')
	
	locationsEC = propagateAndReduceECZ(locationsEC, type)
	
	return locationsEC

def locationsRCM_ECX(type):
	# Error Correction includes the transversal CNOT, transversal Z-basis measurement,
	# and rests due to the Z-error verification measurement.  
	locationsEC = [util.counterUtils.locrest('B', i) for i in range(23)] + \
				  [util.counterUtils.loccnot('A', i, 'B', i) for i in range(23)] + \
				  [util.counterUtils.locZmeas('B', i) for i in range(23)]
	locationsEC = Locations(locationsEC, 'rcm.ecX')
	
	locationsEC = propagateAndReduceECX(locationsEC, type)
	
	return locationsEC

def propagateAndReduceECZ(locations, type=None):
	assert type == 'X' or type == 'Z' or type == None
	# |+> preparations and X-basis measurements cannot cause X errors
	# (and similarly for Z). So they need not be counted.
	if type != None:
		locations = locations.filterAgainst('meas' + type)
		locations = locations.filterAgainst('prep' + type)
	util.counterUtils.propagateAllErrors(locations)
	reduceErrorSyndromesECZ(locations, 'A', 'B')
	
	return locations

def propagateAndReduceECX(locations, type=None):
	assert type == 'X' or type == 'Z' or type == None
	# |+> preparations and X-basis measurements cannot cause X errors
	# (and similarly for Z). So they need not be counted.
	if type != None:
		locations = locations.filterAgainst('meas' + type)
		locations = locations.filterAgainst('prep' + type)
	util.counterUtils.propagateAllErrors(locations)
	reduceErrorSyndromesECX(locations, 'A', 'B')
	
	return locations	

def reduceErrorSyndromesECZ(locations, dataName, zeroName):
	'''
	Reduces errors to syndromes for encoded |psi> (data) on the first block, and encoded |0>
	on the second block.
	'''
	for l in locations:
		reduceZero(l, zeroName)
		reduceGolay(l, dataName)

def reduceErrorSyndromesECX(locations, dataName, plusName):
	'''
	Reduces errors to syndromes for encoded |psi> (data) on the first block, and encoded |+>
	on the second block.
	'''
	for l in locations:
		reducePlus(l, plusName)
		reduceGolay(l, dataName)
	
	

@fetchable
def countECZ_xOnly(prepPairs, settings, noise):
	settingsVZ = settings.getSubcomponent('vz')
	settingsRCM = settings.getSubcomponent('rcm')
	kGoodRCM = settingsRCM['kGood']
	kGood = settings['kGood']
	
	resultVZ = countVerifiedZeroXOnly(prepPairs[0], prepPairs[1], settingsVZ, noise)
	
	logger.info('Computing Z correction counts for X-only and k<={0}'.format(kGoodRCM))
		
	locationsEC = locationsRCM_ECZ('X')
	
	# Count only the errors on block A, the data block.  X errors on block B have no effect.
	countsEC = [countXErrorsGolay(k, locationsEC, 'A', noise) for k in range(kGoodRCM+1)]
	
	kGood = min(len(resultVZ.counts) + len(countsEC) - 2, kGood)
	convolved = convolve(resultVZ.counts, countsEC, kMax=kGood)
	
	prBad = calcPrBadECZ(resultVZ.locTotals, locationsEC.getTotals(), resultVZ.prBad, resultVZ.prAccept, kGoodRCM, kGood, noise['X'])
	
	return CountResult(convolved, None, locationsEC.getTotals() + resultVZ.locTotals, resultVZ.prAccept, prBad)

@fetchable
def countECZ_zOnly(prepPairs, settings, noise):
	'''
	Counts the errors on the two-block Z-error correction, without actually making a correction.
	To complete the correction, convolve this result with the input data counts, then compute
	and apply the correction to data block.
	'''
	
	# Data is block A, correction is block B
	
	settingsVZ = settings.getSubcomponent('vz')
	settingsRCM = settings.getSubcomponent('rcm')
	kGoodRCM = settingsRCM['kGood']
	kGoodEC = settings['kGood']
	
	resultVZ = countVerifiedZeroZOnly(prepPairs[0], prepPairs[1], settingsVZ, noise)
		
	locationsEC = locationsRCM_ECZ('Z')
	
	# There are 2^11 inequivalent Z errors on the encoded |0> block, but 2^12 Z errors on the data block
	countsEC = [countErrors2blocks(k, locationsEC, 'A', 12, 'B', 11, 'Z', noise['Z']) for k in range(kGoodRCM+1)]
	
	kGood = min(len(resultVZ.counts) + len(countsEC) - 2, kGoodEC)
	convolved = convolve(countsEC, resultVZ.counts, kMax=kGood, convolveFcn=convolveABB)
	
	prBad = calcPrBadECZ(resultVZ.locTotals, locationsEC.getTotals(), resultVZ.prBad, resultVZ.prAccept, kGoodRCM, kGood, noise['Z'])
	
	return CountResult(convolved, None, locationsEC.getTotals() + resultVZ.locTotals, resultVZ.prAccept, prBad)


@fetchable
def countECX_zOnly(prepPairs, settings, noise):
	
	settingsVZ = settings.getSubcomponent('vz')
	settingsRCM = settings.getSubcomponent('rcm')
	kGoodRCM = settingsRCM['kGood']
	kGoodEC = settings['kGood']
	
	# Preparing |+> is the same as swapping X and Z errors
	resultPlus = countVerifiedPlusZOnly(prepPairs[0], prepPairs[1], settingsVZ, noise)
	
	logger.info('Computing X correction counts for Z-only and k<={0}'.format(kGoodRCM))
	
	locationsEC = locationsRCM_ECX('Z')
	
	# Count only the errors on block A, the data block.  Z errors on block B have no effect.
	countsEC = [countZerrorsGolay(k, locationsEC, 'A', noise) for k in range(kGoodRCM+1)]
	
	kGood = min(len(resultPlus.counts) + len(countsEC) - 2, kGoodEC)
	convolved = convolve(resultPlus.counts, countsEC, kMax=kGood)

	prBad = calcPrBadECX(resultPlus.locTotals, locationsEC.getTotals(), resultPlus.prBad, resultPlus.prAccept, kGoodRCM, kGood, noise['Z'])
	
	return CountResult(convolved, None, locationsEC.getTotals() + resultPlus.locTotals, resultPlus.prAccept, prBad)

	
@fetchable
def countECX_xOnly(prepPairs, settings, noise):
	'''
	Counts the errors on the two-block X-error correction, without actually making a correction.
	To complete the correction, convolve this result with the input data counts, then compute
	and apply the correction to data block.
	'''

	# TODO: this function recomputes all k whenver ecTotalMaxFailures changes.  Instead,
	# use a fetchable sub-function to compute for a fixed k so that previous results can
	# be reused.

	settingsVZ = settings.getSubcomponent('vz')
	settingsRCM = settings.getSubcomponent('rcm')
	kGoodRCM = settingsRCM['kGood']
	kGood = settings['kGood']

	# Preparing |+> is the same as swapping X and Z errors
	resultPlus = countVerifiedPlusXOnly(prepPairs[0], prepPairs[1], settingsVZ, noise)
	
	locationsEC = locationsRCM_ECX('X')
	countsEC = [countErrors2blocks(k, locationsEC, 'A', 12, 'B', 11, 'X', noise['X']) for k in range(kGoodRCM+1)]
		
	kGood = min(len(resultPlus.counts) + len(countsEC) - 2, kGood)
	convolved = convolve(countsEC, resultPlus.counts, kMax=kGood, convolveFcn=convolveABB)

	prBad = calcPrBadECX(resultPlus.locTotals, locationsEC.getTotals(), resultPlus.prBad, resultPlus.prAccept, kGoodRCM, kGood, noise['X'])

	return CountResult(convolved, None, locationsEC.getTotals() + resultPlus.locTotals, resultPlus.prAccept, prBad)



@fetchable
def countEC_xOnly(prepPairs, settings, noise):
	
	resultECZ = countECZ_xOnly(prepPairs, settings, noise)
	resultECX = countECX_xOnly(prepPairs, settings, noise)
	
	countsECZ = resultECZ.counts
	
	mask11 = (1<<11) - 1
	# Propagate the X errors of the Z correction though the CNOT of the X correction
	for k in range(len(countsECZ)):
		propagated = [0] * (1<<23)
		countList = countsECZ[k]
		for s in range(len(countList)):
			# The bottom block is encoded |+>.  So syndrome parity is ignored.
			sB = mask11 & s
			propagated[(s<<11) + sB] = countList[s]
			
		countsECZ[k] = propagated
	
	convolved = convolve(countsECZ, resultECX.counts, kMax=settings['kGood'], splitListsInto=[2,2])
	
	return CountResult(convolved, 
					   None, 
					   resultECZ.locTotals + resultECX.locTotals,
					   resultECZ.prAccept * resultECX.prAccept,
					   resultECZ.prBad + resultECX.prBad)


@fetchable
def countEC_zOnly(prepPairs, settings, noise):
	kGoodEC = settings['kGood']
	resultECZ = countECZ_zOnly(prepPairs, settings, noise)
	resultECX = countECX_zOnly(prepPairs, settings, noise)
			
	# Convolve the Z errors of the X correction into the top block of the Z correction.
	convolved = convolve(resultECZ.counts, resultECX.counts, kMax=kGoodEC, convolveFcn=convolveABA, extraArgs=[12,11], splitListsInto=[2,2])
	
	return CountResult(convolved, 
					   None, 
					   resultECZ.locTotals + resultECX.locTotals,
					   resultECZ.prAccept * resultECX.prAccept,
					   resultECZ.prBad + resultECX.prBad)
	
	
@fetchable
def countTEC_xOnly(prepPairs, settings, noise):
	'''
	Computes the weighted X-error counts for the trailing error correction.
	
	For the trailing error correction, it is sufficient to consider, for each input
	error equivalence class, the weighted count for logical errors and for the event
	in which no logical error occurs.  The counts returned by this function are indexed
	by [k][syndrome][logical], where "logical" is a boolean indicating that there is
	a logical error (True) or not (False).
		
	@param prepPairs: A list of two pairs of |0> encoding circuits that define ancilla verification.
	@type prepPairs   List of tuples  - [(Locations, Locations), (Locations, Locations)]
	@param settings:  Error correction settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''
	
	# Counts are for two blocks.  The top block has syndromes with 12 bits, the bottom has syndromes
	# with 11 bits.
	resultEC = countEC_xOnly(prepPairs, settings, noise)
	return CountResult(countTEC(resultEC.counts, settings['kGood']), 
					   resultEC.countsRejected,
					   resultEC.locTotals,
					   resultEC.prAccept,
					   resultEC.prBad)


@fetchable
def countTEC_zOnly(prepPairs, settings, noise):
	'''
	Computes the weighted Z-error counts for the trailing error correction.
	
	For the trailing error correction, it is sufficient to consider, for each input
	error equivalence class, the weighted count for logical errors and for the event
	in which no logical error occurs.  The counts returned by this function are indexed
	by [k][syndrome][logical], where "logical" is a boolean indicating that there is
	a logical error (True) or not (False).
		
	@param prepPairs: A list of two pairs of |0> encoding circuits that define ancilla verification.
	@type prepPairs   List of tuples  - [(Locations, Locations), (Locations, Locations)]
	@param settings:  Error correction settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''
	
	# ecCounts are 2-block counts. Block A has 2^12 possible syndromes, block B has 2^11.
	resultEC = countEC_zOnly(prepPairs, settings, noise)
	return CountResult(countTEC(resultEC.counts, settings['kGood']), 
					   resultEC.countsRejected,
					   resultEC.locTotals,
					   resultEC.prAccept,
					   resultEC.prBad)
	

def countTEC(ecCounts, maxFailuresEC):
	'''
	Computes the logical X-error outcome of error correction for each k \in good, and 
	each possible syndrome on the input data block. 
	'''
	
	# Counts are for two blocks.  The top block has syndromes with 12 bits, the bottom has syndromes
	# with 11 bits.
	
	# Precompute syndrome->correction, and syndrome->logical error conversions
	# This has a *very* dramatic affect on performance.
	corrections = [corrector.correctLogicalSyndrome(s) for s in xrange(1<<11)]
	decoded = [corrector.decodeError(corrector.getError(s), 0) for s in xrange(1<<12)]
	
	logicalCounts = [0] * len(ecCounts)	
	minSplitSize = (1<<12)
	
	logger.info('Counting logical effect of error correction for 0 <= k < {0}'.format(len(logicalCounts)))
	
	for k in range(len(logicalCounts)):
				
		sECList = nonZeroIndices(ecCounts[k])
		
		logger.info('Counting Logical EC for k={0}. Correction has {1} syndromes.'.format(k, len(sECList)))
		
		ecCountsK = ecCounts[k]
		
		if len(sECList) < minSplitSize:
			# Computation is small, just perform it directly
			sECSlices = [sECList]
		else:
			# Computation is large.  Delegate to workers
			nSlots = countParallel.numSlots()
			sliceLen = len(sECList) / nSlots
			sECSlices = [sECList[i*sliceLen : (i+1)*sliceLen] for i in range(nSlots)]
			# Assign leftovers to the first worker
			sECSlices[0] += sECList[nSlots*sliceLen:]

		args = [ecCountsK, corrections, decoded]
		results = []
		pool = countParallel.getPool()
		for slice in sECSlices[1:]:
			results.append(pool.apply_async(countTECForSyndromes, [slice] + args))

		# Compute the first slice in this thread
		lCountsK = countTECForSyndromes(*([sECSlices[0]] + args))
		
		# Sum together all of the worker results
		for result in results:
			rCountsK = result.get()
			lCountsK = [[lCountsK[s][l] + rCountsK[s][l] for l in [0,1]] for s in xrange(1<<12)]
		
		logicalCounts[k] = lCountsK
				
	return logicalCounts


def countTECForSyndromes(sECList, ecCounts, corrections, decoded):
	'''
	Computes the logical X/Z-error outcome of error correction for each k \in good, and 
	each possible syndrome on the input data block. 
	'''
	mask11 = ((1<<11)-1)
	
	lCounts = [[0,0] for _ in range(1<<12)]
	
	start = time.time()
	for i, sEC in enumerate(sECList):
		# sECa is 12 bits, sECb is 11 bits
		sECa = sEC>>11
		sECb = sEC & mask11
		countEC = ecCounts[sEC]

		for sData in xrange(1<<12):
			# Propagate the data syndrome through the CNOT of the X/Z correction, then
			# calculate and apply the correction to determine the presence of a logical X/Z-error.
			
			# When syndromes from the data are propagated to logical |+>/|0>, syndromes for which
			# the parities differ but the other 11 bits are the same become equivalent (via the
			# logical X/Z operator).
			sReduced = sData & mask11
			correction = corrections[sReduced ^ sECb]
			sOut = sData ^ sECa ^ correction
			logical = decoded[sOut]
			lCounts[sData][logical] += countEC

		if 0 == i % 10000:
			end = time.time()
			logger.info('i={0} ({1}%), time={2}'.format(i, float(i)/(len(sECList))*100, end-start))
			start = time.time()
			
	return lCounts
					

@fetchable
def countLEC_xOnly(prepPairs, settings, noise):
	'''
	Computes the weighted X-error counts for the leading error correction when
	the syndrome of the input data block is zero.
	
	For the leading error correction, it is sufficient to consider only the
	case in which the syndrome of the input data block is zero.  This is because
	the output syndrome of the EC is independent of the input syndrome, and for
	exRec analysis, the logical state of the LEC output is no important.
		
	@param prepPairs: A list of two pairs of |0> encoding circuits that define ancilla verification.
	@type prepPairs  list of tuples  - [(Locations, Locations), (Locations, Locations)]
	@param settings:  Error correction settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''

	resultEC = countEC_xOnly(prepPairs, settings, noise)
	return CountResult(countLEC(resultEC.counts), 
					   resultEC.countsRejected,
					   resultEC.locTotals,
					   resultEC.prAccept,
					   resultEC.prBad)

@fetchable
def countLEC_zOnly(prepPairs, settings, noise):
	'''
	Computes the weighted Z-error counts for the leading error correction when
	the syndrome of the input data block is zero.
	
	For the leading error correction, it is sufficient to consider only the
	case in which the syndrome of the input data block is zero.  This is because
	the output syndrome of the EC is independent of the input syndrome, and for
	exRec analysis, the logical state of the LEC output is no important.
		
	@param prepPairs: A list of two pairs of |0> encoding circuits that define ancilla verification.
	@type prepPairs  list of tuples  - [(Locations, Locations), (Locations, Locations)]
	@param settings:  Error correction settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''

	resultEC = countEC_zOnly(prepPairs, settings, noise)
	return CountResult(countLEC(resultEC.counts), 
					   resultEC.countsRejected,
					   resultEC.locTotals,
					   resultEC.prAccept,
					   resultEC.prBad)
	
def countLEC(countsEC):	
	# We have the trivial syndrome on the input data block.  So just compute the
	# correction on the second block and apply it to the first block.
	
	corrections = [corrector.correctLogicalSyndrome(s) for s in xrange(1<<11)]
	
	# Malignancy of the rectangle depends only on the syndrome of the output
	# of the LEC, not on the logical state. We still need to keep the 12-bit
	# logical syndrome, but we can take all of the counts for syndromes that
	# decode to a logical error, and fold them into the equivalent syndrome
	# that doesn't decode to a logical error.
	syndromeConversions = [0] * (1<<12)
	# The logical operator for the Golay code is the string of 23 1's
	logicalOp = (1<<23) - 1
	for s in xrange(1<<12):		
		if corrector.decodeSyndrome(s):
			e = corrector.getError(s)
			syndromeConversions[s] = corrector.getLogicalSyndrome(e ^ logicalOp)
		else:
			syndromeConversions[s] = s
	
	counts = [0] * len(countsEC)
	for k in range(len(countsEC)):
		counts[k] = [0] * (1<<12)
		countsECk = countsEC[k]
		sECList = nonZeroIndices(countsECk)
		
		mask11 = ((1<<11)-1)
		for sEC in sECList:
			sECa = sEC>>11
			sECb = sEC & mask11
			countEC = countsECk[sEC]
			
			correction = corrections[sECb]
			s = syndromeConversions[sECa ^ correction]
			counts[k][s] += countEC
	
	return counts