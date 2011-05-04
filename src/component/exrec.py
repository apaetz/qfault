'''
 exRec.py

 This file contains functions for counting malignant sets of X and Z errors 
 produced for CNOT, prep |0>, prep |+>, Z-basis measurement and X-basis measurement
 exRecs.
 
 There are two functions indended for external use:O
 countLEC_xOnly
 countTEC_xOnly
 countLEC_zOnly
 countTEC_zOnly
'''



#===============================================================================
# TODO: this file contains a ton of cut & paste duplication.  It would be nice
# to parameterize and eliminate most of the duplication.
#===============================================================================

from component.ancilla import countVerifiedZeroXOnly, countVerifiedPlusZOnly
from component.ec import countLEC_xOnly, countTEC_xOnly, countEC_xOnly, \
	countLEC_zOnly, countTEC_zOnly, countEC_zOnly
from counting import countParallel
from counting.bounding import computeMax
from counting.countErrors import CountResult, reduceAllErrorSyndromes, \
	countErrors2blocks, reduceAllErrorSyndromesZero, countXerrorsZero, \
	countZerrorsPlus, reduceAllErrorSyndromesPlus
from counting.countParallel import packTasks
from counting.location import Locations
from counting.probability import exRecConfigs, prBadPoly, calcPrBadExRec, \
	prMinFailures, countResultAsPoly, \
	calcPrBadCnotExRec_LEC_CNOT_ignore
from golay import golayCode
from util.cache import fetchable
from util.counterUtils import PartitionIterator
from util.listutils import nonZeroIndices
from util.plotting import plotPolyList
import copy
import logging
import time
import util.counterUtils

logger = logging.getLogger('component.exrec')

corrector = golayCode.Corrector()
#
# exRec
# 

@fetchable
def exRecXOnly(prepPairs, settings):
	settingsEC = settings.getSubcomponent('ec')
	kGoodCnot = settings.getSubcomponent('cnot')['kGood']
	kGoodEC = settingsEC['kGood']
	kGood_LEC_CNOT = settings['kGood-LEC-CNOT']
	kGoodExRec = settings['kGood']
	noise = settings['noise']
	
	# Precompute counts for the leading EC and the logical outcomes of the trailing EC
	countsLEC = countLEC_xOnly(prepPairs, settingsEC, noise)
	countsECLogical = countTEC_xOnly(prepPairs, settingsEC, noise)
	resultEC = countEC_xOnly(prepPairs, settingsEC, noise)

	goodConfigs, _ = exRecConfigs(kGoodEC, kGoodCnot, kGoodExRec, kGood_LEC_CNOT)
	locsCnot = locationsExRecCnot()
	countsC = [countErrors2blocks(k, locsCnot, 'A', 12, 'B', 12, 'X', noise['X']) for k in range(kGoodCnot+1)]
	
	# malignantCounts is a list indexed by [ec][error][k]
	malignantCounts = countCnotExRecParallel(countsLEC.counts, countsECLogical.counts, countsC, copy.copy(goodConfigs))
	
	
	# The malignant counts are indexed as [ec][error][k].
	# Package them up as CountResult objects index as [error][ec]
	totalsEC = resultEC.locTotals
	prBadEC = resultEC.prBad
	prAcceptEC = resultEC.prAccept
	totalsCnot = locsCnot.getTotals()
	prBadCnot = prBadPoly(kGoodCnot, totalsCnot, noise['X'], kGoodExRec)
	
	
	prIgnored = calcPrBadCnotExRec_LEC_CNOT_ignore(totalsEC, totalsCnot, kGoodEC, kGoodCnot, kGood_LEC_CNOT, prAcceptEC, noise['X'])
	
	results = [[0]*4 for _ in range(3)]
	for error in range(3):
		for ec in range(4):
			ecPresent = [1, 1, (ec & 1), (ec >> 1)]
			
			prBad = calcPrBadExRec(totalsEC, totalsCnot, prBadEC, prBadCnot, kGoodExRec, noise['X'], ecPresent, prAcceptEC, prIgnored)
			numECs = sum(ecPresent)
			locTotals = numECs * totalsEC + totalsCnot
			prAccept = prAcceptEC ** numECs 
			
			results[error][ec] = CountResult(malignantCounts[ec][error], 0, locTotals, prAccept, prBad)
	
	return results


def locationsExRecCnot():
	# Transversal CNOT
	locationsC = [util.counterUtils.loccnot('A', i, 'B', i) for i in range(23)]		
	util.counterUtils.propagateAllErrors(locationsC)
	reduceAllErrorSyndromes(locationsC)
	locationsC = Locations(locationsC, 'exRec.cnot')

	return locationsC

@fetchable
def exRecZOnly(prepPairs, settings):
	# Precompute counts for the leading EC and the logical outcomes of the trailing EC
	ecSettings = settings.getSubcomponent('ec')
	kGoodEC = ecSettings['kGood']
	cnotSettings = settings.getSubcomponent('cnot')
	kGoodCnot = cnotSettings['kGood']
	kGood_LEC_CNOT = settings['kGood-LEC-CNOT']
	kGoodExRec = settings['kGood']
	noise = settings['noise']
	
	countsLEC = countLEC_zOnly(prepPairs, ecSettings, noise)
	countsECLogical = countTEC_zOnly(prepPairs, ecSettings, noise)
	goodConfigs, _ = exRecConfigs(kGoodEC, kGoodCnot, kGoodExRec, kGood_LEC_CNOT)
	resultEC = countEC_zOnly(prepPairs, ecSettings, noise)
	
	# exRec counting for Z errors is the same as exRec counting for X errors except that the
	# top (CNOT control) and bottom (CNOT target) blocks are swapped.  In reality this
	# won't matter because the set of configurations being counted should be invariant
	# under swaps of the top and bottom leading EC and top and bottom trailing EC.
	# The CNOT errors are also symmetric with respect to the control and target blocks.
	# But, if for some reason the exRec configurations aren't invariant under the swap,
	# then the the permutation below will allow the Z counting to be done in the same
	# way as X counting.
	
	# Swap the top and bottom leading ECs, and the top and bottom trailing ECs.
	goodConfigs = [[c[1], c[0], c[2], c[4], c[3]] for c in goodConfigs]
	locsCnot = locationsExRecCnot()
	countsC = [countErrors2blocks(k, locsCnot, 'A', 12, 'B', 12, 'Z', noise['Z']) for k in range(kGoodCnot+1)]
	
	malignantCounts = countCnotExRecParallel(countsLEC.counts, countsECLogical.counts, countsC, copy.copy(goodConfigs))
	
	# The malignant counts are indexed as [ec][error][k].
	# Package them up as CountResult objects index as [error][ec]
	totalsEC = resultEC.locTotals
	prBadEC = resultEC.prBad
	prAcceptEC = resultEC.prAccept
	totalsCnot = locsCnot.getTotals()
	prBadCnot = prBadPoly(kGoodCnot, totalsCnot, noise['Z'], kGoodExRec)
	
	# countCnotExRecParallel assumes X-error propagation.  This has the effect
	# of swapping IZ and ZI results, and the TEC configurations.  Swap them back.
	# IZ->ZI, ZI->IZ, ZZ->ZZ
	errorPerm = [1,0,2]
	
	# 00->00, 01->10, 10->01, 11->11
	ecPerm = [0,2,1,3]
	
	prIgnored = calcPrBadCnotExRec_LEC_CNOT_ignore(totalsEC, totalsCnot, kGoodEC, kGoodCnot, kGood_LEC_CNOT, prAcceptEC, noise['Z'])
	
	results = [[0]*4 for _ in range(3)]
	for error in range(3):
		permError = errorPerm[error]
		for ec in range(4):
			ecPresent = [1, 1, (ec & 1), (ec >> 1)]
			permEC = ecPerm[ec]
			
			prBad = calcPrBadExRec(totalsEC, totalsCnot, prBadEC, prBadCnot, kGoodExRec, noise['Z'], ecPresent, prAcceptEC, prIgnored)
			numECs = sum(ecPresent)
			locTotals = numECs * totalsEC + totalsCnot
			prAccept = prAcceptEC ** numECs 
			
			results[permError][permEC] = CountResult(malignantCounts[ec][error], 0, locTotals, prAccept, prBad)
			
	return results


def countCnotExRecParallel(countsLEC, countsECLogical, countsCnot, configs):
	'''
	Counts malignant sets of X errors in the exRec by delegating the work to the pool
	of available worker processes.  Z errors can be counted by appropriately permuting
	the specified error configurations.
	
	Arguments
	---------
	countsLEC		-- A table indexed by [k][syndrome] containing weighted counts for each syndrome at the output
	                   of the leading error correction for fault order k.
	countsECLogical -- A table indexed by [k][syndrome][logical] containing weighted counts for each input syndrome
	              	   to the trailing error correction.  logical=TRUE is the count for the event that ideal decoding
	                   of the output of the TEC yeilds a logical X-error.
	countsCnot		-- A table indexed by [k][syndrome] containing weighted counts for each syndrome at the output
                  	   of the transversal CNOT for fault order k.
	configs			-- The configurations (kLECa, kLECb, kCnot, kTECa, kTECb) to count.
	'''
	
	# Precompute the syndromes that may occur for each component and value k
	sLEClists = [nonZeroIndices(countsLEC[k]) for k in range(len(countsLEC))]
	sClists = [nonZeroIndices(countsCnot[k]) for k in range(len(countsCnot))]
	
	costs = [len(sLEClists[c[0]]) * len(sLEClists[c[1]]) * len(sClists[c[2]]) for c in configs]
	totalCost = sum(costs)
	kMaxExRec = max([sum(c) for c in configs])
	configSlices = packTasks(countParallel.numSlots() * 100, configs, costs)
	
	# Defined so that the updateProgress() closure can modify.
	# See http://www.saltycrane.com/blog/2008/01/python-variable-scope-notes/
	countCnotExRecParallel.setsCompleted = 0
	
	# Callback for worker processes.
	def updateProgress(returnValue):
		cost = returnValue[1]
		countCnotExRecParallel.setsCompleted += cost
								
	results = []
	args = [sLEClists, sClists, countsLEC, countsCnot, countsECLogical]
	logger.info('Delegating computation to {0} workers'.format(len(configSlices)))
	
	from component.cython import exrecCython
	start = time.time()
	pool = countParallel.getPool()
	for cSlice in configSlices:
		result = pool.apply_async(exrecCython.countExRecForConfigs_c, args + [cSlice], callback=updateProgress)
		results.append(result)
		
	while countCnotExRecParallel.setsCompleted < totalCost:
		time.sleep(60)
		completed = countCnotExRecParallel.setsCompleted
		if 0 != completed:
			percentComplete = 100 * float(completed) / totalCost
			ave = (time.time() - start) / (3600 * completed)
			eta = ave * (totalCost - completed)
			logger.info('{0}% complete. ETA {1} hrs'.format(percentComplete, eta))
			
	partialCountsList = [result.get()[0] for result in results]
	
	mCounts = [[[0 for _ in range(kMaxExRec+1)] for _ in range(3)] for _ in range(4)]
	
	for partialCounts in partialCountsList:
		for ec in range(4):
			for e in range(3):
				for k in range(len(partialCounts[ec][e])):
					mCounts[ec][e][k] += partialCounts[ec][e][k]
				

	logger.info('Final Malignant Counts: {0}'.format(mCounts))
			
	return mCounts
	
	
def plotCnotExRecDetails(ancillaPairs, settings, gMaxAlt=1.):
	noise = settings['noise']
	gMin, gMax = noise['X'].noiseRange()
	gMax = min(gMaxAlt, gMax)
	
	zOnly = exRecZOnly(ancillaPairs, settings)	
	xOnly = exRecXOnly(ancillaPairs, settings)
	xMaxes, zMaxes = countCnotExRec(ancillaPairs, settings)
	
	badChars = ['.', ',', ' ', ')', '(', '}', '{']
	pairStr = str(ancillaPairs)
	for char in badChars:
		pairStr = pairStr.replace(char, '')
	
	errorStrX = ['IX', 'XI', 'XX']
	errorStrZ = ['IZ', 'ZI', 'ZZ']
	ecLabels = ['--', '-B', 'A-', 'AB', 'max']
	
	def plotDetail(polys, filename, error, ylabel='', legend='upper left'):
		plotPolyList(polys, gMin, gMax, filename, labelList=ecLabels, numPoints=20, legendLoc=legend, xLabel=r'$\gamma$', yLabel=ylabel)
		
	for error in reversed(range(3)):
		esX = errorStrX[error]
		prBadXList = [result.prBad for result in xOnly[error]]
		plotDetail(prBadXList, 'plot-prBad-cnot-' + esX + pairStr, error, ylabel='Pr[bad]')
		polysX = [countResultAsPoly(r, noise['X']) for r in xOnly[error]] + [xMaxes[error]]
		plotDetail(polysX, 'plot-cnot-exrec-malig-' + esX + pairStr, error, ylabel=r'Pr[mal$_{' + esX + r'}$]')

		#print 'Pr[bad](0.00014)=', [p(0.00014) for p in prBadXList]
		#prAcceptXList = [result.prAccept for result in  xOnly[error]]
		#plotDetail(prAcceptXList, 'plot-prAccept-cnot-' + errorStrX[error] + pairStr, error, legend='upper right')
		
		esZ = errorStrZ[error]
		polysZ = [countResultAsPoly(r, noise['Z']) for r in zOnly[error]] + [zMaxes[error]]
		plotDetail(polysZ, 'plot-cnot-exrec-malig-' + esZ + pairStr, error, ylabel=r'Pr[mal$_{' + esZ + r'}$]')
		prBadZList = [result.prBad for result in zOnly[error]]
		plotDetail(prBadZList, 'plot-prBad-cnot-' + esZ + pairStr, error, ylabel='Pr[bad]')
		#prAcceptZList = [result.prAccept for result in  xOnly[error]]
		#plotDetail(prAcceptZList, 'plot-prAccept-cnot-' + errorStrZ[error] + pairStr, error, legend='upper right')
		
	

@fetchable
def countCnotExRec(ancillaPairs, settings):
	
	noise = settings['noise']
	_, gMax = noise['X'].noiseRange()
	gMin = gMax/100
	
	zOnly = exRecZOnly(ancillaPairs, settings)	
	xOnly = exRecXOnly(ancillaPairs, settings)
	
	# For each logical error, find the EC configuration that
	# gives the maximum error probability.
	xMaxes = []
	zMaxes = []
	errorStrX = ['IX', 'XI', 'XX']
	errorStrZ = ['IZ', 'ZI', 'ZZ']
	ecLabels = ['--', '-B', 'A-', 'AB']
	for error in range(3):
		polysX = [countResultAsPoly(r, noise['X']) for r in xOnly[error]]
		xMax = computeMax(polysX, gMin, gMax)
		xMaxes.append(xMax)
		plotPolyList(polysX + [xMax], gMin, gMax, 'plot-cnot-withMax-' + errorStrX[error], labelList=ecLabels + ['xMax'], numPoints=10)
		
		polysZ = [countResultAsPoly(r, noise['Z']) for r in zOnly[error]]
		zMax = computeMax(polysZ, gMin, gMax)
		zMaxes.append(zMax)
		plotPolyList(polysZ + [zMax], gMin, gMax, 'plot-cnot-withMax-' + errorStrZ[error], labelList=ecLabels + ['zMax'], numPoints=10)

	
	return xMaxes, zMaxes

@fetchable
def countPrepZeroExRec(ancillaPairs, settings):
	noise = settings['noise']
	settingsEC = settings.getSubcomponent('ec')
	settingsVZ = settingsEC.getSubcomponent('vz')	
	kGoodExRec = settings['kGood']
	
	# For |0>, we only need to consider logical X errors.
	prepResult = countVerifiedZeroXOnly(ancillaPairs[0], ancillaPairs[1], settingsVZ, noise)
	prepCounts = prepResult.counts
	ecCountsLogicalX = countTEC_xOnly(ancillaPairs, settingsEC, noise)
	resultEC = countEC_xOnly(ancillaPairs, settingsEC, noise)
	
	maligCountsZero = [0] * (kGoodExRec+1)
	nonZeroPrepSyndromes = [nonZeroIndices(c) for c in prepCounts]
	
	for k in range(len(maligCountsZero)):
		for kPrep, kEC in PartitionIterator(k, 2, [len(prepCounts)-1, len(ecCountsLogicalX.counts)-1]):
			kPrepCounts = prepCounts[kPrep]
			prepSyndromes = nonZeroPrepSyndromes[kPrep]
			kECCountsLogicalX = ecCountsLogicalX.counts[kEC]
			
			maligCountsZero[k] += sum(kPrepCounts[s] * kECCountsLogicalX[s][1] for s in prepSyndromes)
		
	return CountResult(maligCountsZero, 
					   0, 
					   resultEC.locTotals + prepResult.locTotals, 
					   prepResult.prAccept * resultEC.prAccept, 
					   prepResult.prBad + resultEC.prBad)
	
@fetchable
def countPrepPlusExRec(ancillaPairs, settings):
	noise = settings['noise']
	settingsEC = settings.getSubcomponent('ec')
	settingsVZ = settingsEC.getSubcomponent('vz')	
	kGoodExRec = settings['kGood']
	
	# For |+> we only need to consider logical Z errors.
	prepResult = countVerifiedPlusZOnly(ancillaPairs[0], ancillaPairs[1], settingsVZ, noise)
	prepCounts = prepResult.counts
	ecCountsLogical = countTEC_zOnly(ancillaPairs, settingsEC, noise)
	resultEC = countEC_zOnly(ancillaPairs, settingsEC, noise)
	
	maligCounts = [0] * (kGoodExRec+1)
	nonZeroPrepSyndromes = [nonZeroIndices(c) for c in prepCounts]
	
	for k in range(len(maligCounts)):
		for kPrep, kEC in PartitionIterator(k, 2, [len(prepCounts)-1, len(ecCountsLogical.counts)-1]):
			kPrepCounts = prepCounts[kPrep]
			prepSyndromes = nonZeroPrepSyndromes[kPrep]
			kECCountsLogicalX = ecCountsLogical.counts[kEC]
			
			maligCounts[k] += sum(kPrepCounts[s] * kECCountsLogicalX[s][1] for s in prepSyndromes)
		
	return CountResult(maligCounts, 
					   0, 
					   resultEC.locTotals + prepResult.locTotals, 
					   prepResult.prAccept * resultEC.prAccept, 
					   # TODO: Pr[bad] is underestimated here.  need to use calcPrBadExRec().
					   prepResult.prBad + resultEC.prBad)
@fetchable
def countMeasZExRec(ancillaPairs, settings):

	settingsEC = settings.getSubcomponent('ec')
	kGoodMeas = settings.getSubcomponent('cnot')['kGood']
	kGoodExRec = settings['kGood']
	noise = settings['noise']
	
	countsLEC = countLEC_xOnly(ancillaPairs, settingsEC, noise)
	resultEC = countEC_xOnly(ancillaPairs, settingsEC, noise)

	locationsMZ = [util.counterUtils.locZmeas('A', i) for i in range(23)]
	
	util.counterUtils.propagateAllErrors(locationsMZ)
	reduceAllErrorSyndromesZero(locationsMZ)
	locationsMZ = Locations(locationsMZ, 'transversal.measZ')
		
	countsMZ = [countXerrorsZero(k, locationsMZ, 'A', noise) for k in range(kGoodMeas+1)]
	
	maligCountsZbasis = countMeasExRecMalig(countsLEC.counts, countsMZ, kGoodExRec)
	measTotals = locationsMZ.getTotals()
	prBadMeas = prMinFailures(kGoodMeas+1, measTotals, noise['X'], kGoodExRec)
	ecPresent = [True] + [False]*3
	prBad = calcPrBadExRec(resultEC.locTotals, locationsMZ.getTotals(), resultEC.prBad, prBadMeas, kGoodExRec, noise['X'], ecPresent, resultEC.prAccept)
			
	return CountResult(maligCountsZbasis, 0, resultEC.locTotals + locationsMZ.getTotals(), resultEC.prAccept, prBad)

@fetchable
def countMeasXExRec(ancillaPairs, settings):

	settingsEC = settings.getSubcomponent('ec')
	kGoodMeas = settings.getSubcomponent('cnot')['kGood']
	kGoodExRec = settings['kGood']
	noise = settings['noise']
	
	countsLEC = countLEC_zOnly(ancillaPairs, settingsEC, noise)
	resultEC = countEC_zOnly(ancillaPairs, settingsEC, noise)

	locationsMX = [util.counterUtils.locXmeas('A', i) for i in range(23)]
	
	util.counterUtils.propagateAllErrors(locationsMX)
	reduceAllErrorSyndromesPlus(locationsMX)
	locationsMX = Locations(locationsMX, 'transversal.measX')
		
	countsMX = [countZerrorsPlus(k, locationsMX, 'A', noise) for k in range(kGoodMeas+1)]
	
	maligCountsXbasis = countMeasExRecMalig(countsLEC.counts, countsMX, kGoodExRec)
	measTotals = locationsMX.getTotals()
	prBadMeas = prMinFailures(kGoodMeas+1, measTotals, noise['Z'], kGoodExRec)
	ecPresent = [True] + [False]*3
	prBad = calcPrBadExRec(resultEC.locTotals, locationsMX.getTotals(), resultEC.prBad, prBadMeas, kGoodExRec, noise['Z'], ecPresent, resultEC.prAccept)
			
	return CountResult(maligCountsXbasis, 0, resultEC.locTotals + locationsMX.getTotals(), resultEC.prAccept, prBad)

@fetchable
def countRestExRecZero(ancillaPairs, settings):
	
	settingsEC = settings.getSubcomponent('ec')
	kGoodRest = settings.getSubcomponent('cnot')['kGood']
	kGoodExRec = settings['kGood']
	noise = settings['noise']
	
	countsLECx = countLEC_xOnly(ancillaPairs, settingsEC, noise)	
	lCountsTECx = countTEC_xOnly(ancillaPairs, settingsEC, noise)
	resultEC = countEC_xOnly(ancillaPairs, settingsEC, noise)
	
	locationsRest = [util.counterUtils.locrest('A', i) for i in range(23)]
	
	util.counterUtils.propagateAllErrors(locationsRest)
	reduceAllErrorSyndromesZero(locationsRest)
	locationsRest = Locations(locationsRest, 'transversal.rest.zero')
		
	countsX = [countXerrorsZero(k, locationsRest, 'A', noise) for k in range(kGoodRest+1)]
	
	maligCounts = countExRecMalig(countsLECx.counts, countsX, lCountsTECx.counts, kGoodExRec)
	
	restTotals = locationsRest.getTotals()
	prBadRest = prMinFailures(kGoodRest+1, restTotals, noise['X'], kGoodExRec)
	ecPresent = [True] + [False]*3
	prBad = calcPrBadExRec(resultEC.locTotals, restTotals, resultEC.prBad, prBadRest, kGoodExRec, noise['X'], ecPresent, resultEC.prAccept)
	
	return CountResult(maligCounts, resultEC.countsRejected, 2*resultEC.locTotals + restTotals, resultEC.prAccept, prBad)

@fetchable
def countRestExRecPlus(ancillaPairs, settings):
	
	settingsEC = settings.getSubcomponent('ec')
	kGoodRest = settings.getSubcomponent('cnot')['kGood']
	kGoodExRec = settings['kGood']
	noise = settings['noise']
	
	countsLECz = countLEC_zOnly(ancillaPairs, settingsEC, noise)	
	lCountsTECz = countTEC_zOnly(ancillaPairs, settingsEC, noise)
	resultEC = countEC_zOnly(ancillaPairs, settingsEC, noise)
	
	locationsRest = [util.counterUtils.locrest('A', i) for i in range(23)]
	
	util.counterUtils.propagateAllErrors(locationsRest)
	reduceAllErrorSyndromesPlus(locationsRest)
	locationsRest = Locations(locationsRest, 'transversal.rest.plus')
		
	countsZ = [countZerrorsPlus(k, locationsRest, 'A', noise) for k in range(kGoodRest+1)]
	
	maligCounts = countExRecMalig(countsLECz.counts, countsZ, lCountsTECz.counts, kGoodExRec)
	
	
	restTotals = locationsRest.getTotals()
	prBadRest = prMinFailures(kGoodRest+1, restTotals, noise['Z'], kGoodExRec)
	ecPresent = [True] + [False]*3
	prBad = calcPrBadExRec(resultEC.locTotals, restTotals, resultEC.prBad, prBadRest, kGoodExRec, noise['Z'], ecPresent, resultEC.prAccept)
	
	return CountResult(maligCounts, resultEC.countsRejected, 2*resultEC.locTotals + restTotals, resultEC.prAccept, prBad)


	
	

def countExRecMalig(countsLEC, countsGa, lCountsTEC, kMax):	
	nonZeroSyndromesLEC = [nonZeroIndices(c) for c in countsLEC]
	nonZeroSyndromesGa = [nonZeroIndices(c) for c in countsGa]

	decodedSyndromes = [corrector.decodeSyndrome(s) for s in xrange(1<<12)]
	
	maligCounts = [0] * (kMax+1)
	# The measurement is malignant if it produces an uncorrectable error
	# (i.e. a syndrome that decodes as a logical error)
	for k in range(len(maligCounts)):
		logger.info('Computing malignant counts for k={0}'.format(k))
		for kLEC, kGa, kTEC in PartitionIterator(k, 3, [len(countsLEC)-1, len(countsGa)-1, len(lCountsTEC)-1]):
			for sLEC in nonZeroSyndromesLEC[kLEC]:
				maligDecode = not decodedSyndromes[sLEC]
				maligCounts[k] += countsLEC[kLEC][sLEC] * \
								  sum(countsGa[kGa][sGa] * lCountsTEC[kTEC][sLEC ^ sGa][maligDecode] 
									  for sGa in nonZeroSyndromesGa[kGa])
		
	return maligCounts

def countExRecMaligOld(countsLEC, countsGa, lCountsTEC, kMax):	
	nonZeroSyndromesLEC = [nonZeroIndices(c) for c in countsLEC]
	nonZeroSyndromesGa = [nonZeroIndices(c) for c in countsGa]

	decodedSyndromes = [corrector.decodeSyndrome(s) for s in xrange(1<<12)]
	
	maligCounts = [0] * (kMax+1)
	# The measurement is malignant if it produces an uncorrectable error
	# (i.e. a syndrome that decodes as a logical error)
	for k in range(len(maligCounts)):
		logger.info('Computing malignant counts for k={0}'.format(k))
		for kLEC, kGa, kTEC in PartitionIterator(k, 3, [len(countsLEC)-1, len(countsGa)-1, len(lCountsTEC)-1]):
			for sLEC in nonZeroSyndromesLEC[kLEC]:
				maligDecode = not decodedSyndromes[sLEC]
				for sGa in nonZeroSyndromesGa[kGa]:
					print 'kLEC={5}, sLec={0}, decLEC={1}, maligDecode={2}, kGa={6}, sGa={3}, decGa={4}'.format(sLEC, decodedSyndromes[sLEC], maligDecode,
																							 sGa, decodedSyndromes[sGa], kLEC, kGa)
					print 'sLEC ^ sGa={0}, dec(sLEC^sGa)={1}'.format(sLEC^sGa, decodedSyndromes[sLEC^sGa])
					print 'countLEC={0}, countGa={1}, kTEC={2}, lCount={3}'.format(countsLEC[kLEC][sLEC], countsGa[kGa][sGa],
																				   kTEC, lCountsTEC[kTEC][sLEC ^ sGa][maligDecode])
					maligCounts[k] += countsLEC[kLEC][sLEC] * \
								      countsGa[kGa][sGa] * lCountsTEC[kTEC][sLEC ^ sGa][maligDecode]
		
	return maligCounts


def countMeasExRecMalig(countsLEC, countsMeas, kMax):
	
	# The measurement exRec has not TEC.  Or equivalently, its TEC
	# is simply a no-op.
	lCountsTEC = [[0,0] for _ in xrange(1<<12)]
	for s in xrange(1<<12):
		lCountsTEC[s][corrector.decodeSyndrome(s)] = 1
	
	return countExRecMalig(countsLEC, countsMeas, [lCountsTEC], kMax)