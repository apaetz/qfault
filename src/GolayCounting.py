# GolayCounting.py :
#
# Top-level executable script for counting malignant sets and computing the threshold for
# circuits based on the Golay code. 
# 
# Ben Reichardt, 5/22/2010
#
from settings import golayCountSettings
from util.polynomial import Composite, SymPolyWrapper
import logging
from counting.threshold import asymptoticThresh
from counting.levelOne import countExRecs, transformedWeights
from settings.countSettings import getTransformedSettings
from golay.ancillaPrep import getOverlapPreps, getSteaneRandomPreps
from copy import copy
from util.cache import fetchable
from counting import levelOne
import math
from component.ancilla import countVerifiedZeroZOnly, countVerifiedZeroXOnly
from component.xverify import countXVerifyXOnly
from util.plotting import plotList

logger = logging.getLogger('GolayCounting')


def computeConcatenationLevel(gamma, targetNoiseRate, level1Events, epsilon):
	
	if 15* gamma <= targetNoiseRate:
		return 0, targetNoiseRate - 15*gamma
	
	p1s = {}
	for event in level1Events.keys():
		p1s[event] = level1Events[event](gamma)
		
	logger.info('epsilon={0}'.format(epsilon))
	logger.info('p1s={0}'.format(p1s))
	
	def peff(p1, epsilon, k):
		return (epsilon ** max(4*(k-1) - 3, 0)) * p1
	
	cnotEventLabels = ('IX', 'XI', 'XX', 'IZ', 'ZI', 'ZZ')
	
	k = 1
	while True:
		pKs = {}
		for event in level1Events.keys():
			pKs[event] = peff(p1s[event], epsilon, k)
			
		cnotEvents = [pKs.pop(levelOne.maligEventLabels[label]) for label in cnotEventLabels]
		pKs['cnot'] = sum(cnotEvents)
	
		logger.info('P({0}) = {1}'.format(k, pKs))
		delta = targetNoiseRate - max(pKs.values())
		if 0 < delta:
			break
		k += 1
		
	return k, delta


def computeAncillaVerifications(k, K, A, prReject1, prReject2, deltaTarg, epsilon):
	num = math.log(deltaTarg / (K * (4*A)**(K-k)))
	if k > 2:
		den = math.log(epsilon**(4*(k-2) - 3)) + math.log(prReject2)
	if k == 2:
		den = math.log(prReject2)
	else:
		den = math.log(prReject2)
		
	return num / den

def computeGateOverhead(Mk_list, Aprep):
	Averify = 4 * Aprep + 23*3
	Ak = 1
	for m in Mk_list:
		Ak *= (9*23 + Averify * m)
		
	return Ak

def plotOverhead(zeroPrep1, zeroPrep2, zeroPrep3, zeroPrep4, settings):
	
	targetNoiseRates = [1e-12, 1e-10, 1e-9, 1e-6]
	startingNoiseRates = [1e-6, 5e-6, 1e-5, 2e-5, 4e-5, 8e-5, 1e-4, 2e-4, 4e-4, 6e-4, 8e-4, 9e-4, 1e-3, 1.1e-3, 1.2e-3]
	
	level1Events, level2Events, _, prAccept1, prAccept2 = \
		levelOneTwoBounds(zeroPrep1, zeroPrep2, zeroPrep3, zeroPrep4, settings)
	
	# TODO: Hardcoded for overlap prep
	Aprep = 118
	A = 4*Aprep + 9*23
	
	overhead_results = []
	for t in targetNoiseRates:
		overhead_results_t = []
		for p in startingNoiseRates:
			gamma = p/15
			epsilon = computeEpsilon(level1Events, level2Events, gamma)
			level, delta = computeConcatenationLevel(gamma, t, level1Events, epsilon)
			logger.info('Concatentation level for p={0} is {1}'.format(p, level))
			
			prReject1 = sum(1 - pr(gamma) for pr in prAccept1)
			prReject2 = sum(1 - pr(gamma) for pr in prAccept2)
			Mk_list = [computeAncillaVerifications(k, 
												   level, 
												   A,
												   prReject1,
												   prReject2, 
												   delta, 
												   epsilon) for k in range(1,level+1)]
			
			gateOverhead = computeGateOverhead(Mk_list, Aprep)
			logger.info('Gate overhead for t={2} p={0} is {1}'.format(p, gateOverhead, t))
			overhead_results_t.append(gateOverhead)
			
		overhead_results.append(overhead_results_t)
		
	labels = [str(t) for t in targetNoiseRates]
	plotList(startingNoiseRates, overhead_results, filename='gate-overhead-Golay', labelList=labels, xLabel='p', yLabel='gate overhead',
			 xscale='log', yscale='log')
	
def levelOneTwoBounds(zeroPrep1, zeroPrep2, zeroPrep3, zeroPrep4, settings):
		
	pair1 = (zeroPrep1, zeroPrep2)
	pair2 = (zeroPrep3, zeroPrep4)
	ancillaPairs = (pair1, pair2)
	
	settingsVZ = settings.getSubcomponent('ec').getSubcomponent('vz')
	settingsVX = settingsVZ.getSubcomponent('vx')
	noise = settings['noise']
	prAccept1 = [countXVerifyXOnly(zeroPrep1, zeroPrep2, settingsVX, noise).prAccept,
				 countXVerifyXOnly(zeroPrep3, zeroPrep4, settingsVX, noise).prAccept,
				 countVerifiedZeroZOnly(pair1, pair2, settingsVZ, noise).prAccept]
	
	GammaX, weightsX, GammaZ, weightsZ, gMax = transformedWeights(ancillaPairs, settings)

	transformedSettings = getTransformedSettings(settings, GammaX, weightsX, GammaZ, weightsZ, gMax)
	
	settingsVZ = transformedSettings.getSubcomponent('ec').getSubcomponent('vz')
	noise = transformedSettings['noise']
	prAccept2 = [Composite(countXVerifyXOnly(zeroPrep1, zeroPrep2, settingsVX, noise).prAccept, GammaX),
			 	 Composite(countXVerifyXOnly(zeroPrep3, zeroPrep4, settingsVX, noise).prAccept, GammaX),
			 	 Composite(countVerifiedZeroZOnly(pair1, pair2, settingsVZ, noise).prAccept, GammaZ)]
	
	resultsLevel2X, resultsLevel2Z = countExRecs(ancillaPairs, transformedSettings)
	
	# Level two polynomials are given in terms of Gamma.
	# Convert them to polynomials in gamma.
	level1Events = {}
	level2Events = {}
	for event in weightsX.keys():
		level1Events[event] = weightsX[event] * GammaX
		level2Events[event] = Composite(resultsLevel2X[event], GammaX)
	for event in weightsZ.keys():
		level1Events[event] = weightsZ[event] * GammaZ
		level2Events[event] = Composite(resultsLevel2Z[event], GammaZ)
		
	return level1Events, level2Events, gMax, prAccept1, prAccept2

def countAndComputeThresh(zeroPrep1, zeroPrep2, zeroPrep3, zeroPrep4, settings):

	level1Events, level2Events, gMax, _,_ = levelOneTwoBounds(zeroPrep1, zeroPrep2, zeroPrep3, zeroPrep4, settings)
		
	threshold = asymptoticThresh(level1Events, level2Events, gMax)

	logger.info('threshold gamma=%f', threshold)
	return threshold
		

def computeEpsilon(level1Events, level2Events, gamma):
	epsilon = max(level2Events[event](gamma) / level1Events[event](gamma) for event in level1Events.keys())
	return epsilon

#globalSettings = golayCountSettings.getTestSettings()
#globalSettings = golayCountSettings.getMediumSettings()
globalSettings = golayCountSettings.getSettings()	
		
if __name__ == "__main__":
	
	logging.basicConfig(level=logging.INFO,
					format='[%(asctime)s] %(process)d %(levelname)-8s %(name)s: %(message)s')
	import sys

	from counting import countParallel
	nSlots = 1	
	if 1 < len(sys.argv):
		nSlots = int(sys.argv[1])
	countParallel.configureMultiProcess(nSlots)
	
	includeRests = True
	if 3 < len(sys.argv):
		if 'false' == sys.argv[3]:
			includeRests = False
	
	zeroPrep1, zeroPrep2, zeroPrep3, zeroPrep4 = getSteaneRandomPreps(includeRests)
	if 2 < len(sys.argv):
		if 'rand' == sys.argv[2]:
			pass
		elif 'overlap' == sys.argv[2]:
			zeroPrep1, zeroPrep2, zeroPrep3, zeroPrep4 = getOverlapPreps(includeRests)
		else:
			raise Exception('invalid ancilla prep name "%s"', sys.argv[2])
	
	settings = golayCountSettings.getSettings()
	logger.info('Settings are: {0}'.format(settings))
	
#	countAndComputeThresh(zeroPrep1, zeroPrep2, zeroPrep3, zeroPrep4, copy(settings))
	plotOverhead(zeroPrep1, zeroPrep2, zeroPrep3, zeroPrep4, copy(settings))
	print 'done!'	
