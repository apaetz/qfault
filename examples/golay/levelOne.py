'''
Created on Apr 29, 2011

@author: Adam
'''

from qfault.counting.bounding import computeWeights
from qfault.counting.count_errors import CountResult
from qfault.counting.probability import countResultAsPoly, calcPrBad
from qfault.counting.threshold import pseudoThresh, findIntersection
from qfault.util.cache import fetchable
from qfault.util.counterUtils import loccnot
from qfault.util.listutils import addLists
import logging
from qfault.counting import probability
from qfault.qec.error import Pauli

logger = logging.getLogger('levelOne')

maligEventLabels = {
	'IX': r'mal$_{IX}$',
	'XI': r'mal$_{XI}$',
	'XX': r'mal$_{XX}$',
	'IZ': r'mal$_{IZ}$', 
	'ZI': r'mal$_{ZI}$',
	'ZZ': r'mal$_{ZZ}$',
	'rest0': r'mal$_{\mathrm{rest},|0\rangle}$',
	'prep0': r'mal$_{|0\rangle}$',
	'measZ': r'mal$_{\mathrm{meas},Z}$',
	'rest+': r'mal$_{\mathrm{rest},|+\rangle}$',
	'prep+': r'mal$_{|+\rangle}$',
	'measX': r'mal$_{\mathrm{meas},X}$'
}

def countExRecs(ancillaPairs, settings):	
	resultRestZero = countRestExRecZero(ancillaPairs, settings)
	resultRestPlus = countRestExRecPlus(ancillaPairs, settings)
	resultPrepZero = countPrepZeroExRec(ancillaPairs, settings)
	resultPrepPlus = countPrepPlusExRec(ancillaPairs, settings)
	resultMeasZ = countMeasZExRec(ancillaPairs, settings)
	resultMeasX = countMeasXExRec(ancillaPairs, settings)
	polysCnotX, polysCnotZ = countCnotExRec(ancillaPairs, settings)

	resultsX = [resultRestZero, resultPrepZero, resultMeasZ]
	resultsZ = [resultRestPlus, resultPrepPlus, resultMeasX]

	noise = settings['noise']
	resultsX = [countResultAsPoly(r, noise['X']) for r in resultsX] + polysCnotX
	resultsZ = [countResultAsPoly(r, noise['Z']) for r in resultsZ] + polysCnotZ
	
	labels = maligEventLabels
	cnotLabelsX = [labels['IX'], labels['XI'], labels['XX']]
	cnotLabelsZ = [labels['IZ'], labels['ZI'], labels['ZZ']]
	labelsX = [labels['rest0'], labels['prep0'], labels['measZ']] + cnotLabelsX
	labelsZ = [labels['rest+'], labels['prep+'], labels['measX']] + cnotLabelsZ
#	
#	settingsStr = str(settings).replace(' ', '').replace('.', '')
#	gMin, gMax = noise['X'].noiseRange()
##	filename = 'plot-all-' + settingsStr
##	plotPolyList(resultsX + resultsZ, gMin, gMax, filename, labelList=labelsX + labelsZ, numPoints=20)
#	
#	filename = 'plot-cnot-malig-' + settingsStr
#	plotPolyList(polysCnotX + polysCnotZ, gMin, gMax, filename, labelList=cnotLabelsX + cnotLabelsZ, numPoints=10)
	
	return dict(zip(labelsX, resultsX)), dict(zip(labelsZ, resultsZ))

@fetchable
def transformedWeights(ancillaPairs, settings):
	'''
	Calculate transformed noise model reference polynomials
	and weights.
	'''
	noise = settings['noise']

	resultsX, resultsZ = countExRecs(ancillaPairs, settings)	
	pseudothresh = cnotPseudoThresh(ancillaPairs, settings)
	
	results = {}
	results.update(resultsX)
	results.update(resultsZ)
	
		
#	settingsStr = str(settings).replace(' ', '').replace('.', '')
	_, gMax = noise['X'].noiseRange()
	gMax = min(gMax, pseudothresh)
	gMin = gMax/10
#	filename = 'plot-all-' + settingsStr
#	plotPolyList(results.values(), gMin, gMax, filename, labelList=results.keys(), numPoints=10)
#	
#	#GammaX = computeMin(resultsX.values(), gMin, gMax)
#	#ratiosX = computeRatioMaxes(resultsX.values(), GammaX, gMin, gMax)
	GammaX, ratiosX = computeWeights(resultsX.values(), gMin, gMax)
#	print 'ratiosX=', ratiosX
#	
#	#GammaZ = computeMin(resultsZ.values(), gMin, gMax)
#	#ratiosZ = computeRatioMaxes(resultsZ.values(), GammaZ, gMin, gMax)
	GammaZ, ratiosZ = computeWeights(resultsZ.values(), gMin, gMax)
#	print 'ratiosZ=', ratiosZ
#
#	settingsStr = str(settings).replace(' ', '').replace('.', '')
#	filename = 'plot-all-' + settingsStr
#	plotPolyList(results.values() + [GammaX, GammaZ], gMin, gMax, filename, labelList=results.keys() + ['GammaX', 'GammaZ'], numPoints=20)

	
	#filename = 'plot-ratios-' + settingsStr
	#ratioPolys = [Gamma * r for r in ratios]
	#plotPolyList(ratioPolys, gMin, gMax, filename, results.keys(), numPoints=20)

	ratiosX = dict(zip(resultsX.keys(), ratiosX))
	ratiosZ = dict(zip(resultsZ.keys(), ratiosZ))
	
	return GammaX, ratiosX, GammaZ, ratiosZ, gMax


def cnotPseudoThresh(ancillaPairs, settings):
	
	logger.info('Computing CNOT pseudothreshold')
	
	pair1 = ancillaPairs[0]
	pair2 = ancillaPairs[1]
	
	settingsEC = settings.getSubcomponent('ec')
	settingsVZ = settingsEC.getSubcomponent('vz')
	settingsVX = settingsVZ.getSubcomponent('vx')
	noise = settings['noise']
	
	xOnly = exRecXOnly(ancillaPairs, settings)
	zOnly = exRecZOnly(ancillaPairs, settings)	
	prX1 = countXVerifyXOnly(pair1[0], pair1[1], settingsVX, noise).prAccept
	prX2 = countXVerifyXOnly(pair2[0], pair2[1], settingsVX, noise).prAccept
	prZ = countZVerifyZOnly(pair1, pair2, settingsVZ, noise).prAccept 
	
	# For the pseudothreshold, we only care about the case in which both TECs are present
	xOnly = [x[3] for x in xOnly]
	zOnly = [z[3] for z in zOnly]
	
	countsX = addLists(*[r.counts for r in xOnly])
	countsZ = addLists(*[r.counts for r in zOnly])
	
	# Location totals and prAccept should be the same for each logical X-error.  Same for Z.
	# prBad is computed for both at the same time.
	resultX = CountResult(countsX, 0, xOnly[0].locTotals, xOnly[0].prAccept, 0)
	resultZ = CountResult(countsZ, 0, zOnly[0].locTotals, zOnly[0].prAccept, 0)
	
	noiseX = noise['X']
	noiseZ = noise['Z']
	noiseXZ = noise['XZ']
	
	gMin, gMax = noiseXZ.noiseRange()
	cnotLoc = loccnot('A', 0, 'A', 1)
	cnotWeight = sum(noiseXZ.getWeight(cnotLoc, e) for e in range(noiseXZ.numErrors(cnotLoc)))
	pMin = cnotWeight * gMin
	pMax = cnotWeight * gMax
	
	prFailX = countResultAsPoly(resultX, noiseX)
	prFailZ = countResultAsPoly(resultZ, noiseZ)
	prBad = calcPrBad(pair1 + pair2, settings, prX1, prX2, prZ)

	prFailGamma = prFailX + prFailZ + prBad
	prFail = lambda p: prFailGamma(p/cnotWeight)
		
	#thresh = pseudoThresh(prFail, pMin, pMax)
	thresh = pseudoThresh(prFail, pMin, pMax)
	

#	settingsStr = str(ancillaPairs) + str(settings)
#	settingsStr = settingsStr.replace(' ', '').replace('.', '').replace(')','').replace('(','')
#	settingsStr = settingsStr.replace('}', ']').replace('{', '[')	
#	filename = 'plot-full-CNOT-prMalig' + settingsStr
#	plotPolyList([prFail], 
#				 pMin, pMax, 
#				 filename,  
#				 numPoints=20,
#				 xLabel=r'$p$',
#				 yLabel=r'$\Pr[$CNOT malig$]$',
#				 legendLoc='upper left')

	
	print 'pseudothreshold (p) >=', thresh
	
	return thresh/cnotWeight

def pseudoThreshold(counts, locTotals, prBad, prAccept, noise, prInput=1):
	
	#logger.info('Computing CNOT pseudothreshold')
	
	gMin, gMax = noise.noiseRange()
	
	# prAccept is a lower bound on the probability of acceptance.
	# When gamma (gMax) is too high, prAccept could be negative.
	# Avoid this by explicitly checking.
	while 0 >= prAccept(gMax):
		gMax *= .95
	
	cnotLoc = loccnot('A', 0, 'A', 1)
	cnotWeight = sum(noise.getWeight(cnotLoc, e) for e in range(noise.numErrors(cnotLoc)))
	pMin = cnotWeight * gMin
	pMax = min(cnotWeight * gMax, 1)

	summedCounts = []
	for count in counts:
		summedCount = sum(count[key] for key in count.keys() if key != Pauli.I)
		summedCounts.append({None: summedCount})	
		
	print 'summed counts=', summedCounts
		
	countPoly = probability.counts_as_poly(summedCounts, locTotals, noise)
	
	prFailGamma = prInput * (countPoly / prAccept) + prBad
	prFail = lambda p: prFailGamma(p/cnotWeight)
	pthresh = counting.threshold.pseudoThresh(prFail, pMin, pMax, tolerance=1e-5)
	
	#print 'pseudothreshold (p) >=', thresh
	
	return pthresh