
import sim.simulateUtils
import math
from counting import countParallel

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
	sampleErr = sampleSigma / math.sqrt(N)
	
	# Use 1.96 for 95% confidence interval
	return mean, 1.96 * sampleErr 

def printSampleData(data):
	for p in sorted(data.keys()):
		stats = [getStats(sample) for sample in data[p].values()]
		print p, dict(zip(data[p].keys(), stats))