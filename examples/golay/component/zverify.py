'''
 zverify.py

 This file contains functions for counting X and Z errors produced at the
 output of Z-error verification.

 In Z-error verification, an X-verified encoded |0> ancilla on block A is checked
 for Z errors by performing a transversal CNOT to block A from another 
 X-verified encoded |0> ancilla on block B and then doing a transversal Z-basis
 measurement on block B.  The ancilla on block A is accepted only
 if no errors are detected by the measurement.

 There are two functions indended for external use:
 countZVerifyXOnly - Counts X errors
 countZVerifyZOnly - Counts Z errors

 Other functions that may be of use:
 countZVerifyXOnly_uncorrected
 countZVerifyXZ

'''

from util.counterUtils import PartitionIterator, convolvedSum
from counting.probability import calcPrBadVZ, prAcceptZPoly
from component.xverify import countXVerifyXOnly, countXVerifyZOnly, \
	countXVerifyXZ
from counting.countErrors import filterAndPropagate, countXerrorsZero, \
	CountResult, countZerrorsZeroZero, convolveCountsPostselectZ, \
	countXZErrorsZeroZero, convolve_dict, convolveXZPostselectZ, convolveXZRejectedZ
from counting.countParallel import convolve
from counting.location import Locations
from counting.convert import xCountsFromXZCounts, rescaleXZCounts
from util.cache import fetchable
import util.counterUtils
#
# Z-error verification
#

def locationsRCMVZ(type):
	'''
	Construct the Z verification circuit, including rests on the two input blocks
	to account for the postselection measurement during X verification.
	'''
	locationsR0 = [util.counterUtils.locrest('A', i) for i in range(23)]
	locationsR1 = [util.counterUtils.locrest('B', i) for i in range(23)]
	locationsCM = [util.counterUtils.loccnot('B', i, 'A', i) for i in range(23)] + \
				  [util.counterUtils.locXmeas('B', i) for i in range(23)]
	locationsRCM = Locations(locationsR0 + locationsR1 + locationsCM, 'verifyZ-rcm')
	locationsRCM = filterAndPropagate(locationsRCM, type)

	return locationsRCM

@fetchable
def countZVerifyXOnly_uncorrected(prepPairA, prepPairB, settings, noise):
	'''
	Computes the weighted X-error counts for Z-error verification, without 
	any correlation correction.
	
	Only X errors are considered by this function.  To count
	Z errors see countZVerifyZOnly.
	
	In Z-error verification, an encoded |0> ancilla on block A is checked
	for Z errors by performing a transversal CNOT to block A from another 
	encoded |0> ancilla on block B and then doing a transversal X-basis
	measurement on block B.  The ancilla on block A is accepted only
	if no errors are detected by the measurement.
	
	In this case, the encoded |0> ancillas are verified for X errors first,
	and so a pair of encoding circuits is required for each block, A and B.
	
	Since this function considers only X errors, the X-basis measurement
	is effectively ignored.  The resulting weighted error counts are
	upper bounds on the X-error counts that would result from considering
	correlations between X and Z errors.
	
	@param prepPairA: The encoding circuit pair for |0> on block A.
	@type prepPairA:  tuple  - (Locations, Locations)
	@param prepPairB: The encoding circuit pair for |0> on block B.
	@type prepPairB:  tuple  - (Locations, Locations)
	@param settings:  Z-error verification settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''

	settingsVX = settings.getSubcomponent('vx')
	settingsRCM = settings.getSubcomponent('rcm')
	kGoodRCM = settingsRCM['kGood']
	kGoodVZ = settings['kGood']
	kMax = settings.parent()['kGood']
	
	resultA0 = countXVerifyXOnly(prepPairA[0], prepPairA[1], settingsVX, noise)
	resultA1 = countXVerifyXOnly(prepPairB[0], prepPairB[1], settingsVX, noise)
		
	locationsRCM = locationsRCMVZ('X')
		
	# The RCM has locations on two blocks, but X errors on block B don't affect
	# the transversal X-basis measurement and so have no effect.  We only need to
	# count block A.
	countsZ = [countXerrorsZero(k, locationsRCM, 'A', noise) for k in range(kGoodRCM+1)]
		
	# First, convolve the two X verified ancillas.  Then convolve that with the Z verification.
	countsVerified = convolve(resultA0.counts, resultA1.counts, kMax=kGoodVZ)
	countsVerified = convolve(countsVerified, countsZ, kMax=kGoodVZ)
	
	totalsRCM = locationsRCM.getTotals()
	locTotals = resultA0.locTotals + resultA1.locTotals + totalsRCM
	
	# Acceptance is based on Z errors only.
	prAccept = countZVerifyZOnly(prepPairA, prepPairB, settings, noise).prAccept
	
	# Upper bound on Pr[badZ | X1, X2]
	prBad = calcPrBadVZ(resultA0.locTotals, 
					    resultA1.locTotals, 
					    totalsRCM, 
					    resultA0.prBad, 
					    resultA0.prAccept, 
					    resultA1.prBad, 
					    resultA1.prAccept, 
					    kGoodRCM, 
					    kGoodVZ,
					    kMax,
					    noise['X'])
	
	return CountResult(countsVerified, None, locTotals, prAccept, prBad)



@fetchable
def countZVerifyZOnly(prepPairA, prepPairB, settings, noise):
	'''
	Computes the weighted Z-error counts for Z-error verification.
	
	Only Z errors are considered by this function.  To count
	X errors see countZVerifyXOnly.
	
	In Z-error verification, an encoded |0> ancilla on block A is checked
	for Z errors by performing a transversal CNOT to block A from another 
	encoded |0> ancilla on block B and then doing a transversal X-basis
	measurement on block B.  The ancilla on block A is accepted only
	if no errors are detected by the measurement.
	
	In this case, the encoded |0> ancillas are verified for X errors first,
	and so a pair of encoding circuits is required for each block, A and B.
		
	@param prepPairA: The encoding circuit pair for |0> on block A.
	@type prepPairA:  tuple  - (Locations, Locations)
	@param prepPairB: The encoding circuit pair for |0> on block B.
	@type prepPairB:  tuple  - (Locations, Locations)
	@param settings:  Z-error verification settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''

	settingsVX = settings.getSubcomponent('vx')
	settingsRCM = settings.getSubcomponent('rcm')
	kGoodRCM = settingsRCM['kGood']
	kGood = settings['kGood']
	kMax = settings.parent()['kGood']
	
	resultA0 = countXVerifyZOnly(prepPairA[0], prepPairA[1], settingsVX, noise)
	resultA1 = countXVerifyZOnly(prepPairB[0], prepPairB[1], settingsVX, noise)
		
	locationsRCM = locationsRCMVZ('Z')
		
	countsRCM = [countZerrorsZeroZero(k, locationsRCM, 'A', 'B', noise) for k in range(kGoodRCM+1)]
	
	# Convolve together the X-verified ancilla on block B, with the rests, CNOTs and measurements.
	countsZ = convolve(resultA1.counts, countsRCM, kMax=kGood)
	
	countsVerified = convolve(resultA0.counts, countsZ, convolveFcn=convolveCountsPostselectZ, kMax=kGood)

	maxFailures = [len(resultA0.counts)-1, len(countsZ)-1]	
	countsUnverified = [0] * len(countsVerified)
	for k in range(len(countsUnverified)):
		for kA, kZ in PartitionIterator(k, 2, maxFailures):
			countsUnverified[k] += convolvedSum(resultA0.counts[kA], countsZ[kZ])
			
	countsRejected = [countsUnverified[k] - sum(countsVerified[k]) for k in range(len(countsUnverified))]
	
	rcmTotals = locationsRCM.getTotals()
	locTotals = rcmTotals + resultA0.locTotals + resultA1.locTotals
	
	noiseZ = noise['Z']
	# Lower bound on Pr[Z | X1, X2]
	prZ = prAcceptZPoly(countsRejected, kGoodRCM, kGood, rcmTotals, locTotals, resultA0.prAccept, resultA1.prAccept, resultA0.prBad, resultA1.prBad, noiseZ)
	
	# Upper bound on Pr[bad | X1, X2]
	prBad = calcPrBadVZ(resultA0.locTotals, resultA1.locTotals, rcmTotals, resultA0.prBad, resultA1.prAccept, resultA1.prBad, resultA1.prAccept, kGoodRCM, kGood, kMax, noiseZ)

	return CountResult(countsVerified, countsRejected, locTotals, prZ, prBad)

@fetchable
def countZVerifyXZ(prepPairA, prepPairB, settings, noise):
	'''
	Computes the weighted error counts for X-error verification.
	
	In this function, X and Z errors are counted together for low-order faults.
	When X and Z errors are positively correlated, this has the effect of eliminating
	X errors that would otherwise be counted when ignoring Z errors.
	
	@param prepPairA: The encoding circuit pair for |0> on block A.
	@type prepPairA:  tuple  - (Locations, Locations)
	@param prepPairB: The encoding circuit pair for |0> on block B.
	@type prepPairB:  tuple  - (Locations, Locations)
	@param settings:  Z-error verification settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''

	settingsVX = settings.getSubcomponent('vx')
	kBestVZ = settings['kBest']
	
	result0 = countXVerifyXZ(prepPairA[0], prepPairA[1], settingsVX, noise)
	result1 = countXVerifyXZ(prepPairB[0], prepPairB[1], settingsVX, noise)
		
	locationsRCM = locationsRCMVZ(None)
		
	countsRCM = [countXZErrorsZeroZero(k, locationsRCM, 'A', 'B', noise) for k in range(kBestVZ+1)]
	
	# Propagate the X errors on the block B ancilla through the transveral CNOT.
	counts1Prop = [0] * len(result1.counts)
	for k in range(len(counts1Prop)):
		countsProp = {}
		for sB, count in result1.counts[k].iteritems():
			sXa = (sB >> 11) << 34
			countsProp[sXa + sB] = count
			
		counts1Prop[k] = countsProp
	
	# Convolve together the X-verified ancilla on block B, with the rests, CNOTs and measurements.
	countsZ = convolve(counts1Prop, countsRCM, kMax=kBestVZ, convolveFcn=convolve_dict)
	
	acceptedCounts = convolve(result0.counts, countsZ, kMax=kBestVZ, convolveFcn=convolveXZPostselectZ)
	rejectedCounts = convolve(result0.counts, countsZ, kMax=kBestVZ, convolveFcn=convolveXZRejectedZ)
		
	locTotals = result0.locTotals + result1.locTotals + locationsRCM.getTotals()
	return CountResult(acceptedCounts, rejectedCounts, locTotals, 0, 0)

def getXZCorrections(prepPairA, prepPairB, settings, noise, xTotals):
	xzResult = countZVerifyXZ(prepPairA, prepPairB, settings, noise)
	xRejected = xCountsFromXZCounts(xzResult.countsRejected)
	corrections, _ = rescaleXZCounts(xRejected, xzResult.locTotals, xTotals, 'X', noise)
	return corrections


@fetchable
def countZVerifyXOnly(prepPairA, prepPairB, settings, noise):
	'''
	Computes the weighted X-error counts for Z-error verification, including
	low-order correlation corrections.
	
	Only X errors are considered by this function.  To count
	Z errors see countZVerifyZOnly.
	
	In Z-error verification, an encoded |0> ancilla on block A is checked
	for Z errors by performing a transversal CNOT to block A from another 
	encoded |0> ancilla on block B and then doing a transversal X-basis
	measurement on block B.  The ancilla on block A is accepted only
	if no errors are detected by the measurement.
	
	In this case, the encoded |0> ancillas are verified for X errors first,
	and so a pair of encoding circuits is required for each block, A and B.
	
	To account for correlations between X and Z errors, this function
	computes both Z-only error counts and low-order XZ error counts. The
	XZ rejected error counts are converted and then subtracted from the
	Z-only (accepted) error counts.	

		
	@param prepPairA: The encoding circuit pair for |0> on block A.
	@type prepPairA:  tuple  - (Locations, Locations)
	@param prepPairB: The encoding circuit pair for |0> on block B.
	@type prepPairB:  tuple  - (Locations, Locations)
	@param settings:  Z-error verification settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''
				
	goodResult = countZVerifyXOnly_uncorrected(prepPairA, prepPairB, settings, noise)
	
	kBest = settings['kBest']
	
	if kBest > 0:
		corrections = getXZCorrections(prepPairA, prepPairB, settings, noise, goodResult.locTotals)
		
		# Subtract the XZ corrections
		goodCounts = goodResult.counts
		for k in range(len(corrections)):
			goodCounts[k] = [goodCounts[k][s] - corrections[k][s] for s in range(len(goodCounts[k]))]
	
	return goodResult

	