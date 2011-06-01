# GolayCounting.py :
#
# Top-level executable script for counting malignant sets and computing the threshold for
# circuits based on the Golay code. 
# 
# Ben Reichardt, 5/22/2010
#

import os
# By default sympy uses a global cache.  This cache seems to be buggy
# when using multiple threads.  This command turns off the cache.
os.putenv('SYMPY_USE_CACHE', 'no')
print 'Turned off the sympy cache.'

from settings import golayCountSettings
from util.polynomial import Composite
import logging
from counting.threshold import asymptoticThresh
from counting.levelOne import countExRecs, transformedWeights
from settings.countSettings import getTransformedSettings
from golay.ancillaPrep import getOverlapPreps, getSteaneRandomPreps

logger = logging.getLogger('GolayCounting')


def countAndComputeThresh(zeroPrep1, zeroPrep2, zeroPrep3, zeroPrep4, settings):
		
	pair1 = (zeroPrep1, zeroPrep2)
	pair2 = (zeroPrep3, zeroPrep4)
	ancillaPairs = (pair1, pair2)
	
	GammaX, weightsX, GammaZ, weightsZ, gMax = transformedWeights(ancillaPairs, settings)
	transformedSettings = getTransformedSettings(settings, GammaX, weightsX, GammaZ, weightsZ, gMax)
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
		
	threshold = asymptoticThresh(level1Events, level2Events, gMax)

	logger.info('threshold gamma=%f', threshold)
	return threshold
		

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
	
	logger.info('Settings are: {0}'.format(globalSettings))
	
	countAndComputeThresh(zeroPrep1, zeroPrep2, zeroPrep3, zeroPrep4, globalSettings)
	print 'done!'	
