'''
Created on 2010-06-17

@author: adam
'''

from counting import countParallel
from counting.location import LocationCount
from fractions import Fraction
from settings.noise import NoiseModelXSympy, Bound
from util.counterUtils import PartitionIterator, loccnot, locrest, locXprep, \
	locXmeas, locZprep, locZmeas, SliceIterator
import gmpy
import logging
import operator
from util.polynomial import SymPolyWrapper, sympoly1d

logger = logging.getLogger('counting.probability')


def constructLocPoly(kList, nList, prIdealList, prFailList):
	def prKLoc(n, k, prIdeal, prFail):
		B = long(gmpy.comb(n, k))
		#print n, k
		p = B * pow(prIdeal, n - k) * pow(prFail, k)
		return p
			
	prList = [prKLoc(nList[i], kList[i], prIdealList[i], prFailList[i]) for i in range(len(kList))]
	pr = reduce(operator.mul, prList)
	return pr

def sumLocationPartitions(partitions, nList, prList):
	return sum(constructLocLikely(kList, nList, prList) for kList in partitions)

def constructLocLikely(kList, nList, prList):
	def prKLoc(n, k, pr):
		B = long(gmpy.comb(n, k))
		p = B * pow(pr, k)
		return p
			
	prList = [prKLoc(nList[i], kList[i], prList[i]) for i in range(len(kList))]
	pr = reduce(operator.mul, prList)
	return pr


def prMinFailuresOld(kMin, locTotals, noiseModel, kMax=None):
	'''
	Computes the an upper bound on the probability that at least kMin failures occur for the given
	numbers of locations.
	'''
	logger.debug('Computing poly for {0}, {1}, {2}, {3}'.format(kMin, locTotals, noiseModel, kMax))
	
	
	if None == kMax:
		# 10 is arbitrary, but seems to work well.
		kMax = kMin + 10
		bound = True
	else:
		bound = False
		
	nTotal = reduce(operator.add, locTotals)
	nCnot = locTotals.cnot
	nRest = locTotals.rest
	nPrep = locTotals.prepX + locTotals.prepZ
	nMeas = locTotals.measX + locTotals.prepZ
	
	kMax = min(kMax, nTotal) + 1

	cnotLoc = loccnot('A', 0, 'A', 1)
	restLoc = locrest('A', 0)
	prepLoc = locXprep('A', 0)
	measLoc = locXmeas('A', 0)
	
	prIdealCnot = noiseModel.prIdeal(cnotLoc)
	prFailCnot = noiseModel.prFail(cnotLoc)
	
	prIdealRest = noiseModel.prIdeal(restLoc)
	prFailRest = noiseModel.prFail(restLoc)
	
	prIdealPrep = noiseModel.prIdeal(prepLoc)
	prFailPrep = noiseModel.prFail(prepLoc)
	
	prIdealMeas = noiseModel.prIdeal(measLoc)
	prFailMeas = noiseModel.prFail(measLoc)
	
	nList = [nCnot, nRest, nPrep, nMeas]
	prIdealList = [prIdealCnot, prIdealRest, prIdealPrep, prIdealMeas]
	prFailList = [prFailCnot, prFailRest, prFailPrep, prFailMeas]
		
		
	results = []
	for k in range(kMin, kMax):
		iterator = PartitionIterator(k, 4, [nCnot, nRest, nPrep, nMeas])
		results += countParallel.iterParallel(iterator, constructLocPoly, [nList, prIdealList, prFailList])
	
	pr = sum(r.get() for r in results)
		
	if bound and (nTotal >= kMax):
		# Not all possible failure configurations were computed.  Bound
		# the probabiliity by ignoring the prefactor and using probabilities (instead of
		# likelyhoods)
		prIdealCnot = prIdealRest = prIdealPrep = prIdealMeas = 1
		iterator = PartitionIterator(kMax, 4, [nCnot, nRest, nPrep, nMeas]) 
		results = countParallel.iterParallel(iterator, constructLocPoly, [nList, prIdealList, prFailList])
		pr += sum(r.get() for r in results)
	
	# Simplify, if possible.
	try:
		pr = pr.simplify()
	except AttributeError:
		pass
	
	return pr


def prMinFailures(kMin, locations, noiseModel, kMax=None):
	'''
	Returns Pr[kMin <= k <= kMax], an upper bound on the probability that  between kMin and kMax
	failures (of any kind) occur at the given locations.
	If kMax=None, then returns Pr[kMin <= k] an upper bound on the probability that at least kMin
	failures occur.
	
	Arguments
	---------
	kMin		-- The minimum number of failures.
	locTotals	-- (LocationCount) The total number of each type of location.
	noiseModel	-- The noise model.
	kMax		-- (optional) The maximum number of failures.
	'''
	
	#===============================================================================
	# The calculation performed here is of the form:
	#
	# A_\vec{n} sum_\vec{k} L^k \prod_i B_i Pr[fail_i]
	#
	# A_n			-- A prefactor based on the number of each type of location.
	#				   (e.g. (1-12g)^10 (1-8g)^5 (1-4g)^8)
	# L				-- The likelyhood term from the noise model (e.g. g/(1-12g))
	# B_i			-- The binomial coefficient \choose{n_i}{k_i}
	# Pr[fail_i]	-- The probability of failure for location type i. (e.g. 8g)
	#
	# The sum is over all kMin <= k <= kMax.  The product is over all (six) location
	# types. Each product term is computed in parallel.
	#===============================================================================

	bound = Bound.UpperBound
	locTotals = locations.getTotals()
	nTotal = reduce(operator.add, locTotals)
	if None == kMax:
		# 10 is arbitrary, but seems to work well.
		kMax = kMin + 10
		bound = True
	else:
		bound = False
		
	kMax = min(kMax, nTotal) + 1
	
	logger.debug('Computing Pr[{0} <= k < {1}] for {2}, {3}'.format(kMin, kMax, locTotals, noiseModel))
	
	
	nCnot = locTotals.cnot
	nRest = locTotals.rest
	nPrepX = locTotals.prepX
	nPrepZ = locTotals.prepZ
	nMeasX = locTotals.measX
	nMeasZ = locTotals.measZ
	
	cnotLoc = loccnot('A', 0, 'A', 1)
	restLoc = locrest('A', 0)
	prepXLoc = locXprep('A', 0)
	prepZLoc = locZprep('A', 0)
	measXLoc = locXmeas('A', 0)
	measZLoc = locZmeas('A', 0)
	
	nList = [nCnot, nRest, nPrepX, nPrepZ, nMeasX, nMeasZ]
	locsList = [cnotLoc, restLoc, prepXLoc, prepZLoc, measXLoc, measZLoc]
	# Filter out location types that aren't present
	locsList = [locsList[i] for i in range(len(locsList)) if nList[i]]
	nList = [n for n in nList if n]
	
	prIdealList = [noiseModel.prIdeal(l, bound=bound) for l in locsList]
	
	weightList = []
	for l in locsList:
		nErrors = len(noiseModel.errorList(l))
		w = sum(noiseModel.getWeight(l,e, bound=bound) for e in range(nErrors))
		weightList.append(w)

	# This should look something like (1-12g)^nCnot * (1-8g)^nRest * ...
	prefactor = reduce(operator.mul, [prIdealList[i] ** nList[i] for i in range(len(prIdealList))], 1)

		
	# Check for identical weights.  These can be grouped together which
	# will significantly reduce the number of partition iterations.
	weightSet = set(weightList)
	newWeightList = []
	newNList = []
	for w in weightSet:
		n = sum(nList[i] for i in range(len(weightList)) if weightList[i] == w)
		newWeightList.append(w)
		newNList.append(n)
		
	nList = newNList
					
	pr = 0
	likelyhood = noiseModel.likelyhood(bound=bound)
	weights = boundedFailureWeights(kMin, locTotals, noiseModel, kMax, bound)
	for k in range(kMin, kMax):
		pr += weights[k-kMin] * (likelyhood ** k)
		
	logger.debug('A=%s, pr=%s', prefactor, pr)
	
	pr *= prefactor
		
	if bound and (nTotal >= kMax):
		# Not all possible failure configurations were computed.  Bound
		# the probabiliity by ignoring the prefactor and using probabilities (instead of
		# likelyhoods)
		
		iterator = PartitionIterator(kMax, len(nList), nList) 
		prFailList = [noiseModel.prFail(l) for l in locsList]
		results = countParallel.iterParallel(iterator, constructLocLikely, [nList, prFailList])
		cap = sum(r.get() for r in results)
		logger.debug('adding bounding cap: %s', cap)
		pr += cap
	
#	# Simplify, if possible.
#	try:
#		pr = pr.simplify()
#	except AttributeError:
#		pass
	
	return pr



def boundedFailureWeights(kMin, locTotals, noiseModel, kMax, bound=Bound.UpperBound):
	'''
	Calculates a vector of likelyhood weights that can be used to bound Pr[kMin <= k <= kMax].	
	Arguments
	---------
	kMin		-- The minimum number of failures.
	locTotals	-- (LocationCount) The total number of each type of location.
	noiseModel	-- The noise model.
	kMax		-- The maximum number of failures.
	'''
	
	
	nCnot = locTotals.cnot
	nRest = locTotals.rest
	nPrepX = locTotals.prepX
	nPrepZ = locTotals.prepZ
	nMeasX = locTotals.measX
	nMeasZ = locTotals.measZ
	
	cnotLoc = loccnot('A', 0, 'A', 1)
	restLoc = locrest('A', 0)
	prepXLoc = locXprep('A', 0)
	prepZLoc = locZprep('A', 0)
	measXLoc = locXmeas('A', 0)
	measZLoc = locZmeas('A', 0)
	
	nList = [nCnot, nRest, nPrepX, nPrepZ, nMeasX, nMeasZ]
	locsList = [cnotLoc, restLoc, prepXLoc, prepZLoc, measXLoc, measZLoc]
	# Filter out location types that aren't present
	locsList = [locsList[i] for i in range(len(locsList)) if nList[i]]
	nList = [n for n in nList if n]
		
	weightList = []
	for l in locsList:
		w = sum(noiseModel.getWeight(l,e, bound) for e in range(noiseModel.numErrors(l)))
		weightList.append(w)
		
	# Check for identical weights.  These can be grouped together which
	# will significantly reduce the number of partition iterations.
	weightSet = set(weightList)
	newWeightList = []
	newNList = []
	for w in weightSet:
		n = sum(nList[i] for i in range(len(weightList)) if weightList[i] == w)
		newWeightList.append(w)
		newNList.append(n)
		
	weightList = newWeightList
	nList = newNList
				
	results = [0] * (kMax - kMin)
	for k in range(kMin, kMax):
		# Calculation of individual partitions is relatively fast.  Reduce
		# the overhead by grouping partitions together.
		partitioner = PartitionIterator(k, len(nList), nList)
		slicer = SliceIterator(partitioner, 100)
		results[k-kMin] = countParallel.iterParallel(slicer, sumLocationPartitions, [nList, weightList])
	
	weights = [0] * (kMax - kMin)
	for k in range(len(weights)):
		weights[k] = sum(r.get() for r in results[k])
	
	return weights


def countsAsProbability(counts, likelyhood):
	countSums = [sum(countsK.values()) for countsK in counts]
	return sum(countSums[k] * (likelyhood ** k) for k in range(len(countSums)))





def prBadPoly(kGood, locations, noiseModel, kMax=None):
	# Count up all of the locations
	
	prBad = prMinFailures(kGood+1, locations, noiseModel, kMax)
	if int == type(prBad):
		prBad = SymPolyWrapper(sympoly1d([prBad]))
	
	return prBad
	




		

def exRecConfigs(kGoodECTotal, kGoodCnot, kGoodExRec, kGoodLECcnot):
	
	good = []
	bad = []
	maxParts = [kGoodECTotal, kGoodECTotal, kGoodCnot, kGoodECTotal,kGoodECTotal]
	
	golayCorrectable = 3
	
	for k in range(golayCorrectable + 1, kGoodExRec+1):
		for kLECa, kLECb, kC, kTECa, kTECb in PartitionIterator(k, 5, maxParts):
					
			config = [kLECa, kLECb, kC, kTECa, kTECb]
			
			if (kC > kGoodLECcnot[2] and kLECa > kGoodLECcnot[0] and kLECb > kGoodLECcnot[1]) or \
			   (kC > kGoodLECcnot[2] and kLECa > kGoodLECcnot[1] and kLECb > kGoodLECcnot[0]):
				bad.append(config)
			else:
				good.append(config)
				
	# Filter out trivial cases in which there are no failures in the Rec
	good = filter(lambda kList: (kList[2] != 0) or (kList[3] != 0) or (kList[4] != 0), good)

	print '|bad|=', len(bad), '|good|=', len(good)
	
	return good, bad


def calcPrBadVX(prep1Totals, prep2Totals, cmTotals, kGood0, kGoodCM, kGood, kMax, noise):
	'''
	Returns an upper bound on Pr[badX]
	'''
	prBadA1 = prBadPoly(kGood0, prep1Totals, noise, kGood)
	prBadA2 = prBadPoly(kGood0, prep2Totals, noise, kGood)
	return calcPrBadVerify(prep1Totals, prep2Totals, cmTotals, prBadA1, 1, prBadA2, 1, kGoodCM, kGood, kMax, noise)

def calcPrBadVZ(vx1Totals, vx2Totals, rcmTotals, prBadX1, prAcceptX1, prBadX2, prAcceptX2, kGoodRCM, kGood, kMax, noise):
	'''
	Returns an upper bound on Pr[badZ | X1, X2]
	'''
	return calcPrBadVerify(vx1Totals, vx2Totals, rcmTotals, prBadX1, prAcceptX1, prBadX2, prAcceptX2, kGoodRCM, kGood, kMax, noise)

def calcPrBadECZ(zeroTotals, rcmTotals, prBadZero, prAcceptZero, kGoodRCM, kGood, noise):
	'''
	Returns an upper bound on Pr[badECZ | accept |0>]
	
	Arguments
	---------
	veroTotals	-- The number of each type of location in |0> preparation and verification.
	rcmTotals	-- The number of each type of location in Z-error correction.
	prBadZero	-- Pr[bad | X1,X2,Z] An upper bound on the probability that |0> preparation and verification is bad
				   all verifications suceeded.
	prAcceptZero-- Pr[X1, X2, Z] A lower bound on the probability that Z-error verification and both X-error
				   verifications succeed.
	kGoodRCM	-- The maximum number of failures counted in the Z-error correction circuit.
	kGood		-- kGood for the EC
	noise		-- The noise model.
	'''
	prBadRCM = prBadPoly(kGoodRCM, rcmTotals, noise, kGood)
	#prBadECZ = prBadPoly(kMax, zeroTotals + rcmTotals, noise) / prAcceptZero
			
	prBad = prBadZero + prBadRCM# + prBadECZ
	
	return prBad

def calcPrBadECX(plusTotals, rcmTotals, prBadPlus, prAcceptPlus, kGoodRCM, kGood, noise):
	'''
	Returns upper bound on Pr[badECX | accept |+>] 
	
	Arguments
	---------
	plusTotals	-- The number of each type of location in |+> preparation and verification.
	rcmTotals	-- The number of each type of location in X-error correction.
	prBadPlus	-- Pr[bad | Z1,Z2,X] An upper bound on the probability that |+> preparation and verification is bad
				   given that all verifications succeeded.
	prAcceptPlus-- Pr[Z1, Z2, X] A lower bound on the probability that X-error verification and both Z-error
	               verifications succeed.				
	kGoodRCM	-- The maximum number of failures counted in the X-error correction circuit.
	kGood		-- kGood for the EC
	noise		-- The noise model.
	'''
	return calcPrBadECZ(plusTotals, rcmTotals, prBadPlus, prAcceptPlus, kGoodRCM, kGood, noise)

def calcPrBadEC(ecTotals, prBadECZ, prBadECX, prAcceptEC, kGoodEC, kMax, noise):
	'''
	Returns an upper bound on Pr[badEC | accept |0>, accept |+>]
	
	Arguments
	---------
	ecTotals	-- The number of each type of location in the error correction.
	prBadECZ	-- Pr[badECZ | X1, X2, Z] An upper bound on the probability that Z-error correction
	               is bad, given that all |0> ancilla verifications succeeded.
	prBadECX	-- Pr[badECX | Z1, Z2, X] An upper bound on the probability that X-error correction
	               is bad, given that all |+> ancilla verifications succeeded.
	prAcceptEC	-- Pr[accept EC] = Pr[(X1, X2, Z), (Z1, Z2, X)]  A lower bound on the probability that
				   all verifications inside the error correction (|0> verifications and |+> verifications)
				   succeed.
	kGoodEC		-- The maximum number of failures counted in the error correction.
	noise		-- The noise model.
	'''
	prBadEC = prBadPoly(kGoodEC, ecTotals, noise, kMax) / prAcceptEC
	return prBadECZ + prBadECX + prBadEC


def calcPrBadCnotExRec_LEC_CNOT_ignore(totalsEC, totalsCNOT, kGoodEC, kGoodCNOT, kGood_LEC_CNOT, prAcceptEC, noise):
	# Compute Pr of the bad configs, conditioned on the LECs accepting.  
	prIgnored = (2 ** abs(kGood_LEC_CNOT[0] - kGood_LEC_CNOT[1])) * \
				prBadPoly(kGood_LEC_CNOT[2], totalsCNOT, noise, kMax=kGoodCNOT) * \
				prBadPoly(kGood_LEC_CNOT[0], totalsEC, noise, kMax=kGoodEC) * \
				prBadPoly(kGood_LEC_CNOT[1], totalsEC, noise, kMax=kGoodEC)
	
	# Condition on acceptance.
	return prIgnored / (prAcceptEC ** 2)

def calcPrBadExRec(totalsEC, totalsGa, prBadEC, prBadGa, kGoodExRec, noise, ecPresent, prAcceptEC, prIgnored=0):	
	'''
	Returns an upper bound on Pr[badExRec | accept]
	
	Arguments
	---------
	totalsEC	-- The number of each type of location in the error correction. (Assumes all error corrections
				   are identical.)
	totalsGa	-- The number of each type of location in exRec gadget (Ga).
	prBadEC		-- Pr[badEC | acceptEC]  An upper bound on the probability that the error correction is bad
				   given that all ancilla verfications suceeded.
	prBadGa		-- Pr[badGa] An upper bound on the probability that the gadget (Ga) is bad.
	kGoodExRec	-- The maximum number of failures counted in the exRec.
	noise		-- The noise model.
	ecPresent	-- A list [LECa, LECb, TECa, TECb] that specifies which error corrections are actually
				   present in the exRec.  True -> present, False -> not present.
	prAcceptEC	-- Pr[acceptEC] A lower bound on the probability that all ancilla verifications in the error
				   correction succeed.
	badConfigs	-- (optional) A list of tuples (kLECa, kLECb, kCnot, kTECa, kTECb) of failure configurations for which
				   k <= kGoodExRec, but were not counted (because they are hard).				  
	'''
	
	numECs = sum(ecPresent)
	locTotals = numECs * totalsEC + totalsGa

	
	# When considering the whole exRec, we must condition on acceptance of all ancillas used in error correction.
	prBad = prBadPoly(kGoodExRec, locTotals, noise) / (prAcceptEC ** numECs)
		
	prBad += numECs * prBadEC + prBadGa + prIgnored 
	
	return prBad
	
	
def prBadForIgnoredConfigs(badConfigs, locTotalsEC, locTotalsCnot, noise):
	
	kMax = max([max(c) for c in badConfigs] + [-1])
	
	# For each 0 <= k <= kMax, compute an upper bound on the probability that exactly k failures occur.
	prEC = [prMinFailures(k, locTotalsEC, noise, kMax=k) for k in range(kMax+1)]
	prCnot = [prMinFailures(k, locTotalsCnot, noise, kMax=k) for k in range(kMax+1)]
	
	pr = 0
	for kLECa, kLECb, kC, kTECa, kTECb in badConfigs:
		prLECa = prEC[kLECa]
		prLECb = prEC[kLECb]
		prC = prCnot[kC]
		prTECa = prEC[kTECa]
		prTECb = prEC[kTECb]
		
		pr += prLECa * prLECb * prC * prTECa * prTECb
		
	return pr 
	


	

def calcPrBadVerify(a1Totals, a2Totals, verifyTotals, prBadA1, prAcceptA1, prBadA2, prAcceptA2, kGoodD, kGood, kMax, noise):
	'''
	Returns an upper bound on Pr[bad | accept A1, accept A2]
	
	Arguments:
	----------
	a1Totals	-- The number of each type of location in ancilla 1 preparation.
	a2Totals	-- The number of each type of location in ancilla 1 preparation.
	verifyTotals-- The number of each type of location in the verification circuit.
	prBadA1		-- Pr[bad, accept] An upper bound on the probability of the bad event and acceptance in ancilla 1 preparation. 
	prAcceptA1	-- Pr[accept] A lower bound on the probability of acceptance of ancilla 1.
	prBadA1		-- Pr[bad, accept] An upper bound on the probability of the bad event and acceptance in ancilla 2 preparation. 
	prAcceptA1	-- Pr[accept] A lower bound on the probability of acceptance of ancilla 2.
	kGoodD		-- The maximum number of failures counted in the detection circuit (rests, CNOT, measurement).
	kGood		-- The maximum number of overall failures allowed before the component is considered bad.
	kMax		-- The maximum number of failures to consider. (i.e. The maximum number of failures counted by the next largest component.)
	noise		-- The noise model.
	'''
	prBadD = prBadPoly(kGoodD, verifyTotals, noise, kGood)
	prBadTotal = prBadPoly(kGood, a1Totals + a2Totals + verifyTotals, noise, kMax)
	
	# The bad event is when:
	# 1. More than kGood failures occur in total.
	# 2. Either of the ancillas are bad, or the detection circuit is bad.
	
	# First, all of the events that are conditioned on acceptance of A1 and A2
	prBad = prBadTotal/(prAcceptA1 * prAcceptA2) + prBadA1/prAcceptA1 + prBadA2/prAcceptA2
	
	# Now add in the events that don't depend on acceptance of A1 or A2
	prBad += prBadD 
	
	return prBad


def calcThresh(countsList, 
			   rejectedX0, 
			   rejectedX1, 
			   rejectedZ, 
			   kGoodZero, 
			   kGoodV, 
			   kGoodEC, 
			   kGoodECTotal, 
			   kGoodCnot, 
			   kGoodExRec, 
			   configs):
	'''
	Computes the asymptotic threshold based on counts of logicalal 1-Rec failures.
	'''
	pseudoThresh = lambda c: calcPseudoThresh(c, 
											  rejectedX0, 
											  rejectedX1, 
											  rejectedZ, 
											  kGoodZero, 
											  kGoodV, 
											  kGoodEC, 
											  kGoodECTotal, 
											  kGoodCnot, 
											  kGoodExRec, 
											  configs,
											  pEff=pE2X)
	
	pseudoThreshList = [pseudoThresh(count) for count in countsList]
	
	return max(pseudoThreshList)
	
	
def pE2X(p):
	return 4*p/15

def pId(p):
	return p

def calcPseudoThresh(xCounts,
					 zCounts, 
					 rejectedX0, 
					 rejectedX1, 
					 rejectedZ, 
					 kGoodZero, 
					 kGoodV, 
					 kGoodEC, 
					 kGoodECTotal, 
					 kGoodCnot, 
					 kGoodExRec, 
					 configs,
					 pEff=pId):
	'''
	Computes the pseudo-threshold for a given set of counts.
	'''
	
	#TODO: rewrite for p1 polynomial.
	
	p1 = 0
	pStep = float(Fraction(1, 10000))
	p = float(Fraction(9,10000) - pStep)
	
	while p1 < pEff(p):
		p += pStep
		print 'p = ', float(p), 'pEff = ', float(pEff(p))		
		p1 = calcP1(p, xCounts, zCounts, rejectedX0, rejectedX1, rejectedZ, kGoodZero, kGoodV, kGoodEC, kGoodECTotal, kGoodCnot, kGoodExRec, configs)
		
	print p - pStep, p1
	
	return p - pStep

def calcP1(prMaligX, prMaligZ, prAccept, prBad):
	p1 = (prMaligX + prMaligZ) / prAccept + prBad		
	return p1
	


def likelyhoodPrefactorPoly(locTotals, noise):
	nCnot = locTotals.cnot
	nRest = locTotals.rest
	nPrepX = locTotals.prepX
	nPrepZ = locTotals.prepZ
	nMeasX = locTotals.measX
	nMeasZ = locTotals.measZ
	
	cnotLoc = loccnot('A', 0, 'A', 1)
	restLoc = locrest('A', 0)
	prepXLoc = locXprep('A', 0)
	measXLoc = locXmeas('A', 0)
	prepZLoc = locZprep('A', 0)
	measZLoc = locZmeas('A', 0)
	
	prefactor = noise.prIdeal(cnotLoc) ** nCnot
	prefactor *= noise.prIdeal(restLoc) ** nRest
	prefactor *= noise.prIdeal(prepXLoc) ** nPrepX
	prefactor *= noise.prIdeal(measXLoc) ** nMeasX
	prefactor *= noise.prIdeal(prepZLoc) ** nPrepZ
	prefactor *= noise.prIdeal(measZLoc) ** nMeasZ
				
	return prefactor
	

def countsToPoly(counts, locTotals, noise):
	'''
	Convert weighed error likelyhood counts into a polynomial in gamma (= p/15).
	'''
	
	prefactor = likelyhoodPrefactorPoly(locTotals, noise)
	countPr =  countsAsProbability(counts, noise.likelyhood())
	return prefactor * countPr

	

def prAcceptXPoly(rejectedCounts, kGood0, kGoodCM, locTotalsPrep1, locTotalsPrep2, locTotalsCM, noise):
	'''
	Returns a lower bound on X-verification acceptance.
	'''
		
	# Note, the Pr[bad] calculated here is not exactly the Pr[bad] calculated
	# for counts of errors that get accepted (i.e. calcPrBadZero, calcPrBadX).  
	# Here, we assume that anything that wasn't counted is rejected, and therefore bad.
	prBad1 = prMinFailures(kGood0+1, locTotalsPrep1, noise)
	prBad2 = prMinFailures(kGood0+1, locTotalsPrep2, noise)
	prBadX = prMinFailures(kGoodCM+1, locTotalsCM, noise) + prBad1 + prBad2
	
	prAccept = lowerBoundPoly(rejectedCounts, prBadX, locTotalsPrep1 + locTotalsPrep2 + locTotalsCM, noise)
	
	return prAccept
	


def prAcceptZPoly(counts, kGoodRCM, kGoodVZ, totalsRCM, totals, prAcceptX0, prAcceptX1, prBadX0, prBadX1, noise):
	'''
	Returns a lower bound on Pr[Z | X1, X1], the probability of Z-error verification acceptance.
	
	Arguments
	---------
	counts		-- A list of weighted sums for which Z-error verification fails.  Indexed by k, the number of failures.
	kGoodRCM	-- 
	kGoodVZ		--
	totalsRCM	--
	totals		--
	prAcceptX0	-- Pr[X0] A lower bound on the probability that X-error verification succeeds for ancilla 1.
	prAcceptX1	-- Pr[X1] A lower bound on the probability that X-error verification succeeds for ancilla 1.
	prBadX0		-- Pr[badX1] An upper bound on the probability of the (unconditioned) bad event for X-error verification
				   on ancilla 1.
	prBadX1		-- Pr[badX2] An upper bound on the probability of the (unconditioned) bad event for X-error verification
				   on ancilla 2.
	noise		-- The noise model.
	'''
		
	# Pr[bad] here is not exactly the same as Pr[bad] for the case of successful Z-error verification.
	prBadVZ = prMinFailures(kGoodVZ+1, totals, noise) + prMinFailures(kGoodRCM+1, totalsRCM, noise)
	prBadX = prBadX0 / prAcceptX0 + prBadX1 / prAcceptX1
	prBadZ = prBadVZ + prBadX
	
	normFactor = 1 / (prAcceptX0 * prAcceptX1)
	
	prAccept = lowerBoundPoly(counts, prBadZ, totals, noise, normFactor=normFactor)
	
	return prAccept


def lowerBoundPoly(rejectCounts, pBad, locTotals, noise, normFactor=1):
	'''
	Returns a lower bound on Pr[accept], the probability of acceptance.
	'''
	
	prReject = upperBoundPoly(rejectCounts, pBad, locTotals, noise, normFactor)		
	return 1 - prReject


def upperBoundPoly(counts, pBad, locTotals, noise, normFactor=1):
	'''
	Returns an upper bound on Pr[!accept], the probability
	that verification detects an error and aborts the ancilla preparation.
	This probability is computed by Pr[!accept] <= Pr[!accept, good] + Pr[bad]
	'''

	pr = countsToPoly(counts, locTotals, noise)
	pr = pr * normFactor + pBad
	return pr

def countResultAsPoly(result, noise):
	'''
	Returns a polynomial for the sum of the probability masses 
	corresponding to the given count result and noise model. 
	'''
	prefactor = likelyhoodPrefactorPoly(result.locTotals, noise)
	likelyhood = countsAsProbability(result.counts, noise.likelyhood())
	
	return (prefactor * likelyhood) / result.prAccept + result.prBad


def calcPrBad(ancillaPreps, settings, prX1, prX2, prZ):
	
	# Assume that all ancilla preps contain the same number of each location.
	prepTotals = ancillaPreps[0].getTotals()
	if not all([prepTotals == prep.getTotals() for prep in ancillaPreps]):
		print ancillaPreps
		raise Exception
	
	settingsEC = settings.getSubcomponent('ec')
	settingsVZ = settingsEC.getSubcomponent('vz')
	settingsVX = settingsVZ.getSubcomponent('vx')
	
	noise = settings['noise']['XZ']
	
	kGood0 = settingsVX.getSubcomponent('zero')['kGood']
	kGoodCM_VX = settingsVX.getSubcomponent('cm')['kGood']
	kGoodVX = settingsVX['kGood']
	kGoodVZ = settingsVZ['kGood']
	kGoodEC = settingsEC['kGood']
	kGoodExRec = settings['kGood']
	
	# LocationCount is initialized as:
	# LocationCount(cnot, prepX, prepZ, measX, measZ, rest)
	cmTotals = LocationCount(23, 0, 0, 0, 23, 0)
	print 'badVX1:'
	prBadX1 = calcPrBadVX(prepTotals, prepTotals, cmTotals, kGood0, kGoodCM_VX, kGoodVX, kGoodVZ, noise)
	prBadX1 /= prX1
	print 'badVX2:'
	prBadX2 = calcPrBadVX(prepTotals, prepTotals, cmTotals, kGood0, kGoodCM_VX, kGoodVX, kGoodVZ, noise)
	prBadX2 /= prX2
	
	vxTotals = 2*prepTotals + cmTotals
	rcmTotals = LocationCount(23, 0, 0, 23, 0, 23*2)
	kGoodRCM = settingsVZ.getSubcomponent('rcm')['kGood']

	print 'badVZ:'
	prBadZ = calcPrBadVZ(vxTotals, vxTotals, rcmTotals, prBadX1, prX1, prBadX2, prX2, kGoodRCM, kGoodVZ, kGoodEC, noise)
	prBadZ /= prZ
	
	prAccept = prX1 * prX2 * prZ
	zeroTotals = 2*vxTotals + rcmTotals
	rcmTotals = LocationCount(23, 0, 0, 23, 0, 23)
	kGoodRCM = settingsEC.getSubcomponent('rcm')['kGood']

	print 'badECZ:'
	prBadECZ = calcPrBadECZ(zeroTotals, rcmTotals, prBadZ, prAccept, kGoodRCM, kGoodEC, noise)
	totalsECZ = zeroTotals + rcmTotals

	plusTotals = zeroTotals.dual()
	rcmTotals = LocationCount(23, 0, 0, 0, 23, 23)
	prBadECX = calcPrBadECZ(zeroTotals, rcmTotals, prBadZ, prAccept, kGoodRCM, kGoodEC, noise)

	totalsECX = plusTotals + rcmTotals
	
	totalsEC = totalsECZ + totalsECX
	prAcceptEC = prAccept ** 2
	prBadEC = calcPrBadEC(totalsEC, prBadECZ, prBadECX, prAcceptEC, kGoodEC, kGoodExRec, noise)
	
	totalsCnot = LocationCount(23, 0, 0, 0, 0, 0)
	kGoodCnot = settings.getSubcomponent('cnot')['kGood']
	kGood_LEC_Cnot = settings['kGood-LEC-CNOT']

	ecPresent = [1]*4
	prBadCnot = prBadPoly(kGoodCnot, totalsCnot, noise, kGoodExRec)
	prIgnored = calcPrBadCnotExRec_LEC_CNOT_ignore(totalsEC, totalsCnot, kGoodEC, kGoodCnot, kGood_LEC_Cnot, prAcceptEC, noise)
	prBadExRec = calcPrBadExRec(totalsEC, totalsCnot, prBadEC, prBadCnot, kGoodExRec, noise, ecPresent, prAcceptEC, prIgnored)
	return prBadExRec

if __name__ == '__main__':
	
	logging.basicConfig(level=logging.INFO,
					format='[%(asctime)s] %(process)d %(levelname)-8s %(name)s: %(message)s')
	
	import sys
	
	nSlots = 1	
	if 1 < len(sys.argv):
		nSlots = int(sys.argv[1])
		
	print 'using {0} slots'.format(nSlots)
	
	countParallel.configureMultiProcess(nSlots)
	
	locTotals = LocationCount(1000, 40, 40, 40, 23, 23)
	prBad = prMinFailures(3, locTotals, NoiseModelXSympy(), kMax=11)
	print prBad(0.002/15)
	
	
	#settings = getMediumSettings()
#	settings = getSettings()
##	settings['noise'] = FixedPointNoise(settings['noise'], 0.001/15)
#	settings['noise'] = NoiseModelX(gmpy.mpf(0.001/15), gmpy.mpf(0.002/15))
#	settings = buildSettingsTree(settings)
#	logger.info('Computing prBad')
#	prBad = calcPrBad(settings, gmpy.mpf(0.84), gmpy.mpf(0.84), gmpy.mpf(0.68))
#	logger.info('done.')
#	print prBad(0.001/15)
#	print prBad(0.002/15)
#	logger.info('Computing derivative')
#	dPrBad = prBad.diff()
#	print dPrBad(0.001/15)
#	print dPrBad(0.002/15)
#	logger.info('done')
#	n,d = dPrBad.asNumDen()
#	logger.info('computed num/den')
#	n1 = n(0.001/15)
#	n2 = n(0.002/15)
#	d1 = d(0.002/15)
#	d2 = d(0.002/15)
#	print n1, n2
#	print d1, d2
#	print n1/d1
#	print n
#	print d
			
	