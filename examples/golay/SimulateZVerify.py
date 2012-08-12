'''
Created on 2010-09-10

@author: adam
'''
from golay import golayCode
from golaySim import prepareAndVerify, prepareUnverifiedAncilla, \
	prepareAndVerifyNoPostselect
from sim.simulateUtils import rest, cnot, measX
import golay.ancillaPrep
import sim.simulateUtils
import util



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
errorRates.cnot = 23e-4
errorRates.prepX = errorRates.prepZ = errorRates.measX = errorRates.measZ = 4/15. * errorRates.cnot
errorRates.rest = 12/15. * errorRates.cnot
#errorRates.rest = 0

def prepA0():
    errors, _ = prepareAndVerify(errorRates, 'X', 
        'a0', lambda: prepareUnverifiedAncilla(errorRates, rounds0, cnotsA0, 'a0'), 
        'a1', lambda: prepareUnverifiedAncilla(errorRates, rounds1, cnotsA1, 'a1'),
        includeMeasRest=True)
    return errors
def prepA2():
    errors, _ = prepareAndVerify(errorRates, 'X', 
        'a2', lambda: prepareUnverifiedAncilla(errorRates, rounds2, cnotsA2, 'a2'), 
        'a3', lambda: prepareUnverifiedAncilla(errorRates, rounds3, cnotsA3, 'a3'),
        includeMeasRest=True)
    return errors


def prepA():
    return prepareAndVerifyNoPostselect(errorRates, 'Z',
                     'a0', prepA0, 
                     'a2', prepA2,
                     includeMeasRest=False)[0]
                     
                     
def simZRCM():
    errors = {'X':{'a':0, 'b':0},'Z':{'a':0, 'b':0}}
    
    for i in range(23):
        rest(errorRates, errors, 'a', i)
        rest(errorRates, errors, 'b', i)
        cnot(errorRates, errors, 'b', i, 'a', i)
        measX(errorRates, errors, 'b', i)

    return errors

def checkZRCM():
    
    counts = [0] * (1<<11)

    for _ in xrange(1000):
        errors = simZRCM()
        eZ = corrector.hashError(errors['Z']['b'], True)
        sZ = corrector.getSyndrome(eZ)
        
        #s = (sX << 11) + sZ
        counts[sZ] += 1

    sList = util.listutils.nonZeroIndices(counts)
    countTotal = sum(counts[s] for s in sList)
    prSim = [float(c)/countTotal for c in counts]
    prZErrorSim = sum(prSim[1:])
    prAcceptSim = prSim[0]
    print 'Pr[Z-error] (sim) = ', prZErrorSim
    print 'Pr[acceptZ] (sim) = ', prAcceptSim 

if __name__ == '__main__':
    
    corrector = golayCode.Corrector()
#    checkZRCM()
#    assert False
    
    
    counts = [0] * (1<<11)
    weightCounts = [0] * 5
    
    for _ in xrange(1000):
        errors = prepA()
        eZ = corrector.hashError(errors['Z']['a2'], True)
        sZ = corrector.getSyndrome(eZ)
        counts[sZ] += 1
        
    sList = util.listutils.nonZeroIndices(counts)
    countTotal = sum(counts[s] for s in sList)
    
    prSim = [float(c)/countTotal for c in counts]
    prZErrorSim = sum(prSim[1:])
    prAcceptSim = prSim[0]
    print 'Pr[Z-error] (sim) = ', prZErrorSim
    print 'Pr[acceptZ] (sim) = ', prAcceptSim

        
#    cutoff = 1e-5
#    for s in xrange(len(prCount)):
#        if prCount[s] >= cutoff or prSim[s] >= cutoff: 
#            print '{0},{1},{2}'.format(s, prCount[s], prSim[s])
