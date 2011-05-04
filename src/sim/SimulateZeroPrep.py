'''
Created on 2010-09-10

@author: adam
'''
from golaySim import prepareAndVerify, prepareUnverifiedAncilla,\
    prepareAndVerifyNoPostselect
import util.counterUtils
from golay import golayCode
import golay.ancillaPrep
import sim.simulateUtils

prepCkts = golay.ancillaPrep.prepCircuits()
cnotsA0 = prepCkts[8]
rounds0 = [0, 6, 2, 4, 1, 3, 5]
cnotsA1 = prepCkts[70]
rounds1 = [1, 4, 5, 6, 0, 3, 2]
cnotsA2 = prepCkts[45]
rounds2 = [6, 5, 1, 3, 2, 4, 0]
cnotsA3 = prepCkts[99]
rounds3 = [0, 5, 1, 3, 6, 4, 2]

errorRates = sim.simulateUtils.ErrorRates()
errorRates.cnot = 1e-3
errorRates.prepX = errorRates.prepZ = errorRates.measX = errorRates.measZ = 4/15. * errorRates.cnot
errorRates.rest = 12/15. * errorRates.cnot
#errorRates.rest = 0

def prepA0():
    errors, _ = prepareAndVerifyNoPostselect(errorRates, 'X', 
        'a0', lambda: prepareUnverifiedAncilla(errorRates, rounds0, cnotsA0, 'a0'), 
        'a1', lambda: prepareUnverifiedAncilla(errorRates, rounds1, cnotsA1, 'a1'),
        includeMeasRest=False)
    return errors
def prepA2():
    errors, _ = prepareAndVerify(errorRates, 'X', 
        'a2', lambda: prepareUnverifiedAncilla(errorRates, rounds2, cnotsA2, 'a2'), 
        'a3', lambda: prepareUnverifiedAncilla(errorRates, rounds3, cnotsA3, 'a3'),
        includeMeasRest=False)
    return errors







if __name__ == '__main__':
    
    corrector = golayCode.Corrector()
    counts = [0] * (1<<11)
    weightCounts = [0] * 5
    
    for _ in xrange(1000):
        #errors = prepareUnverifiedAncilla(errorRates, rounds0, cnotsA0, 'a0')
        #errors = prepareUnverifiedAncilla(errorRates, rounds1, cnotsA1, 'a0')
        errors = prepA0()
        eZ = corrector.reduceError(errors['Z']['a0'], True)
        sZ = corrector.getSyndrome(eZ)
        counts[sZ] += 1
        
    sList = util.counterUtils.nonZeroIndices(counts)
    countTotal = sum(counts[s] for s in sList)
    
    prSim = [float(c)/countTotal for c in counts]


    prZErrorSim = sum(prSim[1:])
    
    print 'Pr[Z-error] (sim) = ', prZErrorSim
        
#    cutoff = 1e-5
#    for s in xrange(len(prCount)):
#        if prCount[s] >= cutoff or prSim[s] >= cutoff: 
#            print '{0},{1},{2}'.format(s, prCount[s], prSim[s])
