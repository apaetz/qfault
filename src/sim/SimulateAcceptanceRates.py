
from golaySim import prepareAncillaNew, prepareAncillaOverlap
import sim.simulateUtils
from sim import golaySim


def simAncillaPrep(prepFcn, pMin, pMax, pStep, iters, noRests=False):

	errorRates = sim.simulateUtils.ErrorRates()
	p = pMin
	
	while p <= pMax:
		errorRates.cnot = p
		errorRates.prepX = errorRates.prepZ = errorRates.measX = errorRates.measZ = 4/15. * errorRates.cnot
		errorRates.rest = 12/15. * errorRates.cnot
		if noRests:
			errorRates.rest = 0
		
		golaySim.initPrepStats()
		for _ in xrange(iters):
			prepFcn(errorRates, 'prep0', 'Z')
		acceptRates = golaySim.getAcceptanceRates()
		print '{0},{1}'.format(p, acceptRates)

		p += pStep


pMin = 9e-4
pMax = 23e-4
pStep = 1e-4
iters = 10000

#print 'Steane Random'
#simAncillaPrep(prepareAncillaNew, pMin, pMax, pStep, iters)
#print 'Steane Random: No rests'
#simAncillaPrep(prepareAncillaNew, pMin, pMax, pStep, iters, noRests=True)

print 'Overlap'
simAncillaPrep(prepareAncillaOverlap, pMin, pMax, pStep, iters)
print 'Overlap: No rests'
simAncillaPrep(prepareAncillaOverlap, pMin, pMax, pStep, iters, noRests=True)

#print 'Steane Reichardt'
#simAncillaPrep(golayCode.prepareAncilla, pMin, pMax, pStep, iters)
#print 'Steane Reichardt: No rests'
#simAncillaPrep(golayCode.prepareAncilla, pMin, pMax, pStep, iters, noRests=True)