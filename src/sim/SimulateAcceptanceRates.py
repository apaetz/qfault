
from golaySim import prepareAncillaNew, prepareAncillaOverlap
import sim.simulateUtils
from util.cache import fetchable
import logging
import math
from sim.golaySim import prepareAncilla12
from counting.countParallel import iterParallel, asyncFuncWrapper
from counting import countParallel


@fetchable
def simOverlapPrep(pMin, pMax, pStep, iters, noRests=False):
	return simAncillaPrep(prepareAncillaOverlap, pMin, pMax, pStep, iters, noRests)

@fetchable
def simSteane4Prep(pMin, pMax, pStep, iters, noRests=False):
	return simAncillaPrep(prepareAncillaNew, pMin, pMax, pStep, iters, noRests)

@fetchable
def simSteane12Prep(pMin, pMax, pStep, iters, noRests=False):
	return simAncillaPrep(prepareAncilla12, pMin, pMax, pStep, iters, noRests)



def simAncillaPrep(prepFcn, pMin, pMax, pStep, iters, noRests=False):

	errorRates = sim.simulateUtils.ErrorRates()
	p = pMin
	
	sampleData = {}
	
	while p <= pMax:
		errorRates.cnot = p
		errorRates.prepX = errorRates.prepZ = errorRates.measX = errorRates.measZ = 4/15. * errorRates.cnot
		errorRates.rest = 12/15. * errorRates.cnot
		if noRests:
			errorRates.rest = 0
		
		pool = countParallel.getPool()
		results = [pool.apply_async(prepFcn, [errorRates, 'prep0', 'Z']) for _ in xrange(iters)]
		attemptSamples = [r.get()[1] for r in results]
		
		samples = dict()
		for key in attemptSamples[0].keys():
			samples[key] = [0]*iters
			
		for i, sample in enumerate(attemptSamples):
			for key, attempts in sample.iteritems():
				samples[key][i] = attempts
				
		sampleData[p] = samples
		p += pStep
		
	return sampleData

def getStats(sample):
	N = float(len(sample))
	mean = sum(sample) / N
	sampleSigma = math.sqrt(sum((s - mean)**2 for s in sample) / (N-1))
	
	return mean, sampleSigma 

def printSampleData(data):
	for p in sorted(data.keys()):
		stats = [getStats(sample) for sample in data[p].values()]
		print p, dict(zip(data[p].keys(), stats))
	
if __name__ == '__main__':
	logging.basicConfig(level=logging.INFO)
	pMin = 0
	pMax = 2e-3
	pStep = 1e-4
	iters = 2
	
	
	#print 'Steane Random'
	#simAncillaPrep(prepareAncillaNew, pMin, pMax, pStep, iters)
	#print 'Steane Random: No rests'
	#simAncillaPrep(prepareAncillaNew, pMin, pMax, pStep, iters, noRests=True)
	
	print 'Overlap'
	data = simOverlapPrep(pMin, pMax, pStep, iters)
	printSampleData(data)
	#print 'Overlap: No rests'
	#simAncillaPrep(prepareAncillaOverlap, pMin, pMax, pStep, iters, noRests=True)
	
	#print 'Steane Reichardt'
	#simAncillaPrep(golayCode.prepareAncilla, pMin, pMax, pStep, iters)
	#print 'Steane Reichardt: No rests'
	#simAncillaPrep(golayCode.prepareAncilla, pMin, pMax, pStep, iters, noRests=True)