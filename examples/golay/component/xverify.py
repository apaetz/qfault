'''
 xverify.py

 This file contains functions for counting X and Z errors produced at the
 output of X-error verification.

 In X-error verification, an encoded |0> ancilla on block A is checked
 for X errors by performing a transversal CNOT from block A to another 
 encoded |0> ancilla on block B and then doing a transversal Z-basis
 measurement on block B.  The ancilla on block A is accepted only
 if no errors are detected by the measurement.

 There are two functions indended for external use:
 countXVerifyXOnly - Counts X errors
 countXVerifyZOnly - Counts Z errors

 Other functions that may be of use:
 countXVerifyZOnly_uncorrected
 countXVerifyXZ

'''

from counting.convert import zCountsFromXZCounts, rescaleXZCounts
from counting.countErrors import filterAndPropagate, propagateAndReduceZeroX, \
	countXerrorsZero, countXerrorsZeroZero, convolveABB, convolveCountsPostselectX, \
	CountResult, propagateAndReduceZeroXZ, countXZErrorsZero, countXZErrorsZeroZero, \
	convolveXZPostselectX, convolveXZRejectedX, countZerrorsZero, \
	propagateReduceAndCountZero
from counting.countParallel import convolve
from counting.location import Locations
from counting.probability import prAcceptXPoly, calcPrBadVX
from util.cache import fetchable
from util.counterUtils import PartitionIterator, convolvedSum
import logging
import util.counterUtils

logger = logging.getLogger('component.xverify')
	
def locationsCMVX(type):
	'''Returns the CNOT and measurement locations for X-error verification.'''
	locationsCM = [util.counterUtils.loccnot('A', i, 'B', i) for i in range(23)] + \
				  [util.counterUtils.locZmeas('B', i) for i in range(23)]
	locationsCM = Locations(locationsCM, 'verifyX-CM')
	locationsRCM = filterAndPropagate(locationsCM, type)

	return locationsRCM

def countXVerifyXOnly(zeroPrepA, zeroPrepB, settings, noise):	
	'''
	Computes the weighted X-error counts for X-error verification.
	
	Only X errors are considered by this function.  To count
	Z errors see countXVerifyZOnly_uncorrected.
	
	In X-error verification, an encoded |0> ancilla on block A is checked
	for X errors by performing a transversal CNOT from block A to another 
	encoded |0> ancilla on block B and then doing a transversal Z-basis
	measurement on block B.  The ancilla on block A is accepted only
	if no errors are detected by the measurement.
	
	@param zeroPrepA: The encoding circuit for |0> on block A.
	@type zeroPrepA:  Locations
	@param zeroPrepB: The encoding circuit for |0> on block B.
	@type zeroPrepB:  Locations
	@param settings:  X-error verification settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''
	zeroSettings = settings.getSubcomponent('zero')
	kGood0 = zeroSettings['kGood']
	cmSettings = settings.getSubcomponent('cm')
	kGoodCM = cmSettings['kGood']
	kGood = settings['kGood']
	kMax = settings.parent()['kGood']
						
	return countXVerifyXOnly_fetch(zeroPrepA, zeroPrepB, kGood0, kGoodCM, kGood, kMax, noise)

@fetchable
def countXVerifyXOnly_fetch(zeroPrepA, zeroPrepB, kGood0, kGoodCM, kGood, kMax, noise):	
	'''
	Computes the weighted X-error counts for X-error verification.
	
	Only X errors are considered by this function.  To count
	Z errors see countXVerifyZOnly_uncorrected.
	
	In X-error verification, an encoded |0> ancilla on block A is checked
	for X errors by performing a transversal CNOT from block A to another 
	encoded |0> ancilla on block B and then doing a transversal Z-basis
	measurement on block B.  The ancilla on block A is accepted only
	if no errors are detected by the measurement.
	
	@param zeroPrepA: The encoding circuit for |0> on block A.
	@type zeroPrepA:  Locations
	@param zeroPrepB: The encoding circuit for |0> on block B.
	@type zeroPrepB:  Locations
	@param kGood0:    kGood for encoded |0> preparation
	@type kGood0:	  int
	@param kGoodCM:	  kGood for the CNOT and measurement
	@type kGoodCM:	  int     
	@param kGood:     kGood for X-error verification
	@type kGood:      int
	@param kMax:	  The maximum k to consider (i.e., kGood for Z-error verification)
	@type kMax:       int
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''
	
	# |+> preparations and X-basis measurements cannot cause X errors
	# so there is no use in counting them.
	locsPrepA = propagateAndReduceZeroX(zeroPrepA)
	locsPrepB = propagateAndReduceZeroX(zeroPrepB)
	locsCM = locationsCMVX('X')
	
	countsA1 = [countXerrorsZero(k, locsPrepA, 'A', noise) for k in range(kGood0+1)]
	countsA2 = [countXerrorsZero(k, locsPrepB, 'B', noise) for k in range(kGood0+1)]
	countsVX = [countXerrorsZeroZero(k, locsCM, 'A', 'B', noise) for k in range(kGoodCM+1)]
	
	countsVX = convolve(countsVX, countsA2, convolveFcn=convolveABB, kMax=kGood)
	
	# Convolve and postselect simultaneously.
	countsVerified = convolve(countsA1, countsVX, convolveFcn=convolveCountsPostselectX, kMax=kGood)
	
	# In order to compute a lower bound on the probability of acceptance,
	# we need a sum of the counts for each k *before* postselection.
	countsUnverified = [0] * len(countsVerified)
	for k in range(len(countsVerified)):
		for kA, kVX in PartitionIterator(k, 2, [len(countsA1)-1, len(countsVX)-1]):
			countsUnverified[k] += convolvedSum(countsA1[kA], countsVX[kVX])
			
	countsRejected = [countsUnverified[k] - sum(countsVerified[k]) for k in range(len(countsUnverified))]
					
	# Now compute a polynomial upper bound on the for the events that weren't counted.
	totals1 = locsPrepA.getTotals()
	totals2 = locsPrepB.getTotals()
	totalsCM = locsCM.getTotals()
	locTotals = totals1 + totals2 + totalsCM
	
	noiseX = noise['X']
	prAccept = prAcceptXPoly(countsRejected, kGood0, kGoodCM, totals1, totals2, totalsCM, noiseX)
	
	# Upper bound on Pr[bad, X]
	prBad = calcPrBadVX(totals1, totals2, totalsCM, kGood0, kGoodCM, kGood, kMax, noiseX)
						
	return CountResult(countsVerified, countsRejected, locTotals, prAccept, prBad)





@fetchable
def countXVerifyZOnly_uncorrected(zeroPrepA, zeroPrepB, settings, noise):
	'''
	Computes the weighted Z-error counts for X-error verification, without 
	any correlation correction.
	
	Only Z errors are considered by this function.  To count
	X errors see countXVerifyXOnly.
	
	In X-error verification, an encoded |0> ancilla on block A is checked
	for X errors by performing a transversal CNOT from block A to another 
	encoded |0> ancilla on block B and then doing a transversal Z-basis
	measurement on block B.  The ancilla on block A is accepted only
	if no errors are detected by the measurement.
	
	Since this function considers only Z errors, the Z-basis measurement
	is effectively ignored.  The resulting weighted error counts are
	upper bounds on the Z-error counts that would result from considering
	correlations between X and Z errors.
	
	@param zeroPrepA: The encoding circuit for |0> on block A.
	@type zeroPrepA:  Locations
	@param zeroPrepB: The encoding circuit for |0> on block B.
	@type zeroPrepB:  Locations
	@param settings:  X-error verification settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''
	
	kGood0 = settings.getSubcomponent('zero')['kGood']
	settingsCM = settings.getSubcomponent('cm')
	kGoodCM = settingsCM['kGood']
	kGood = settings['kGood']
	kMax = settings.parent()['kGood']
	
	# X-error verification has three parts, two encoded |0> preparations
	# and a transversal CNOT and Z-basis measurement on block B.
	locsPrepA, countsPrepA = propagateReduceAndCountZero(zeroPrepA, 'Z', 'A', noise, kGood0)
	locsPrepB, countsPrepB = propagateReduceAndCountZero(zeroPrepB, 'Z', 'B', noise, kGood0)
	
	# CNOT and Z-basis measurements
	locsCM = locationsCMVX('Z')
		
	# For CNOT and measurement, count only the errors on block A.  Errors from
	# block B do not propagate (to block A).
	countsC = [countZerrorsZero(k, locsCM, 'A', noise) for k in range(kGoodCM + 1)]
		
	countsVerify = convolve(countsPrepB, countsC, kMax=kGood)
	countsVerified = convolve(countsPrepA, countsVerify, kMax=kGood)
	
	totals1 = locsPrepA.getTotals()
	totals2 = locsPrepB.getTotals()
	totalsCM = locsCM.getTotals()
	locTotals = totals1 + totals2 + totalsCM

	# Acceptance is based on X errors only.
	prAccept = countXVerifyXOnly(zeroPrepA, zeroPrepB, settings, noise).prAccept
	
	# Upper bound on Pr[bad, X]	
	prBad = calcPrBadVX(totals1, totals2, totalsCM, kGood0, kGoodCM, kGood, kMax, noise['Z'])

	return CountResult(countsVerified, None, locTotals, prAccept, prBad)


@fetchable
def countXVerifyXZ(zeroPrepA, zeroPrepB, settings, noise):
	'''
	Computes the weighted error counts for X-error verification.
	
	In this function, X and Z errors are counted together.  When X and
	Z errors are positively correlated, this has the effect of eliminating
	Z errors that would otherwise be counted when ignoring X errors.
		
	@param zeroPrepA: The encoding circuit for |0> on block A.
	@type zeroPrepA:  Locations
	@param zeroPrepB: The encoding circuit for |0> on block B.
	@type zeroPrepB:  Locations
	@param settings:  X-error verification settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''	
	
	kBest = settings['kBest']

	zeroPrepA = propagateAndReduceZeroXZ(zeroPrepA)
	
	locsCM = [util.counterUtils.loccnot('A', i, 'B', i) for i in range(23)]\
					+ [util.counterUtils.locZmeas('B', i) for i in range(23)]
	
	# Combine |0> on block B with X-verification.
	verifyX = propagateAndReduceZeroXZ(zeroPrepB + Locations(locsCM, 'verifyX-CM'))	
		
	countsA = [countXZErrorsZero(k, zeroPrepA, 'A', noise) for k in range(kBest+1)]
	countsVX = [countXZErrorsZeroZero(k, verifyX, 'A', 'B', noise) for k in range(kBest+1)]
	
	logger.debug('|0> A count sums: %s', [sum(c.values()) for c in countsA])
	logger.debug('VX count sums: %s', [sum(c.values()) for c in countsVX])
	logger.info('Syndrome counts A1: %s', [len(countsA[k]) for k in range(len(countsA))])
	logger.info('Syndrome counts VX: %s', [len(countsVX[k]) for k in range(len(countsVX))])
	
	# splitListsInto=[1,1] avoids splitting the large dictionaries of counts, which is very time consuming.
	acceptedCounts = convolve(countsA, countsVX, kMax=kBest, convolveFcn=convolveXZPostselectX, splitListsInto=[1,1])
	rejectedCounts = convolve(countsA, countsVX, kMax=kBest, convolveFcn=convolveXZRejectedX, splitListsInto=[1,1])

	return CountResult(acceptedCounts, rejectedCounts, zeroPrepA.getTotals() + verifyX.getTotals(), 0, 0)





def getXZCorrections(zeroPrepA, zeroPrepB, settings, noise, zTotals):
	xzResult = countXVerifyXZ(zeroPrepA, zeroPrepB, settings, noise)
	zRejected = zCountsFromXZCounts(xzResult.countsRejected)
	
	corrections, _ = rescaleXZCounts(zRejected, xzResult.locTotals, zTotals, 'Z', noise)
	return corrections

@fetchable
def countXVerifyZOnly(zeroPrepA, zeroPrepB, settings, noise):
	'''
	Computes the weighted Z-error counts for X-error verification, including
	low-order correlation corrections.
	
	Only Z errors are considered by this function.  To count
	X errors see countXVerifyXOnly.
	
	In X-error verification, an encoded |0> ancilla on block A is checked
	for X errors by performing a transversal CNOT from block A to another 
	encoded |0> ancilla on block B and then doing a transversal Z-basis
	measurement on block B.  The ancilla on block A is accepted only
	if no errors are detected by the measurement.

	To account for correlations between X and Z errors, this function
	computes both Z-only error counts and low-order XZ error counts. The
	XZ rejected error counts are converted and then subtracted from the
	Z-only (accepted) error counts.	
	
	@param zeroPrepA: The encoding circuit for |0> on block A.
	@type zeroPrepA:  Locations
	@param zeroPrepB: The encoding circuit for |0> on block B.
	@type zeroPrepB:  Locations
	@param settings:  X-error verification settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''

						
	goodResult = countXVerifyZOnly_uncorrected(zeroPrepA, zeroPrepB, settings, noise)
	
	kBest = settings['kBest']
	
	if kBest > 0:
		corrections = getXZCorrections(zeroPrepA, zeroPrepB, settings, noise, goodResult.locTotals)
		
		# Subtract the XZ corrections
		goodCounts = goodResult.counts
		for k in range(len(corrections)):
			goodCounts[k] = [goodCounts[k][s] - corrections[k][s] for s in range(len(goodCounts[k]))]
		
		logger.debug('Corrected count sums: %s', [sum(c) for c in goodCounts])
		
	return goodResult
	
	
