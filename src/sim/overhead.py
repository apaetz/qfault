
import sim.simulateUtils
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