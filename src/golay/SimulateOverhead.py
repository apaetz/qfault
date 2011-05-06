from sim.golaySim import prepareAncilla12, prepareAncillaOverlap, \
	prepareAncillaNew
from sim.overhead import simAncillaPrep, getStats
from util.cache import fetchable
from util.plotting import plotList
import logging
import numpy
import math
from util import stats, listutils

@fetchable
def simOverlapPrep(pMin, pMax, pStep, iters, noRests=False):
	stats = {'prepA0':0, 'prepA2':0, 'prepA':0}
	return simAncillaPrep(prepareAncillaOverlap, pMin, pMax, pStep, iters, noRests)

@fetchable
def simSteane4Prep(pMin, pMax, pStep, iters, noRests=False):
	return simAncillaPrep(prepareAncillaNew, pMin, pMax, pStep, iters, noRests)

@fetchable
def simSteane12Prep(pMin, pMax, pStep, iters, noRests=False):
	return simAncillaPrep(prepareAncilla12, pMin, pMax, pStep, iters, noRests)

def overhead4(data, measX, measZ):
	
	X = sorted(data.keys())
	
	means = []
	errors = []
	for p in X:
		x1 = data[p]['prepA0']
		x2 = data[p]['prepA2']
		z = data[p]['prepA']
		
		overheadSample = [measX * (x1[s] + x2[s]) + measZ * z[s] for s in xrange(len(z))]
		mean, _, err = getStats(overheadSample)		
		means.append(mean)
		errors.append(err)
		
	return X, means, errors

def prAccept4(data):
	
	X = sorted(data.keys())
	
	means = []
	errors = []
	for p in X:
		x1 = data[p]['prepA0']
		x2 = data[p]['prepA2']
		z = data[p]['prepA']
		
		# For each z-error verification we have one successful X1 and X2.
		xSuccessCount = sum(z)
		
		nX1 = sum(x1)
		samplesX1 = [1]*xSuccessCount + [0]*(nX1-xSuccessCount)
		meanX1 = stats.mean(samplesX1)
		seX1 = stats.stdErr(stats.std(samplesX1), nX1)
				
		nX2 = sum(x2)
		samplesX2 = [1]*xSuccessCount + [0]*(nX2-xSuccessCount)
		meanX2 = stats.mean(samplesX2)
		seX2 = stats.stdErr(stats.std(samplesX2), nX2)

		
		nZ = sum(z)
		samplesZ = [1]*iters + [0]*(nZ-iters)
		meanZ = stats.mean(samplesZ)
		seZ = stats.stdErr(stats.std(samplesZ), nZ)
	
		se = stats.stdProduct([seX1,seX2,seZ], [meanX1,meanX2,meanZ])
		err95 = stats.err95(se) 		
				
		means.append(meanX1*meanX2*meanZ)
		errors.append(err95)
		
	return X, means, errors


def prAccept12(data):
	
	X = sorted(data.keys())
	
	means = []
	errors = []
	for p in X:
		x1Attempts = data[p]['prepA_X1']
		z1Attempts = data[p]['prepA_Z1']
		x2Attempts = data[p]['prepA_X2']
		z2Attempts = data[p]['prepA_Z2']
		x3Attempts = data[p]['prepA_X3']
		z3Attempts = data[p]['prepA_Z3']
		x4Attempts = data[p]['prepA_X4']
		z4Attempts = data[p]['prepA_Z4']

		# For each verification stage, calculate the total number of attempts
		# and the number of successes.
		x4 = (sum(x4Attempts), iters)
		z4 = (sum(z4Attempts), x4[0])
		x3 = (sum(x3Attempts), z4[0])
		z3 = (sum(z3Attempts), x3[0])
		x2 = (sum(x2Attempts), z4[0])
		z2 = (sum(z2Attempts), x4[0])
		x1 = (sum(x1Attempts), 2*z3[0]+x2[0])
		z1 = (sum(z1Attempts), x3[0]+z2[0])		
		attSucc = [x1,z1,x2,z2,x3,z3,x4,z4]
		
		samples = [[1]*succ + [0]*(att-succ) for att,succ in attSucc]
		sampleMeans = map(stats.mean, samples)
		stds = map(stats.std, samples)
		stdErrs = [stats.stdErr(stds[i], attSucc[i][0]) for i in xrange(len(stds))]
		
		err95 = stats.stdProduct(stdErrs, sampleMeans)
		
		means.append(listutils.mul(sampleMeans))
		errors.append(err95)
		
	return X, means, errors


#
#def overhead6(data, measX, measZ):
#	
#	prAcceptList = [getyvals(data, i) for i in range(1,6)]
#	X = getxvals(data)
#	
#	overheadList = []
#	for p in range(len(X)):
#		x1 = measX / prAcceptList[0][p]
#		x2 = measX / prAcceptList[1][p]
#		z1 = (x1 + x2 + measZ) / prAcceptList[2][p]
#		z2 = measX / prAcceptList[3][p]
#		x3 = (z1 + z2 + measZ) / prAcceptList[4][p]
#		
#		overheadList.append(x3)
#		
#	return X, overheadList


def overhead12(data, meas1, meas2, meas3, meas4):

	X = sorted(data.keys())
	
	def calc(x1,z1,x2,z2,x3,z3,x4,z4):
		o1 = meas1 * (x1 + z1)
		o2 = meas2 * (x2 + z2)
		o3 = meas3 * (x3 + z3)
		o4 = meas4 * (x4 + z4)
		
		return o1 + o2 + o3 + o4
	
	means = []
	errors = []
	for p in X:
		x1 = data[p]['prepA_X1']
		z1 = data[p]['prepA_Z1']
		x2 = data[p]['prepA_X2']
		z2 = data[p]['prepA_Z2']
		x3 = data[p]['prepA_X3']
		z3 = data[p]['prepA_Z3']
		x4 = data[p]['prepA_X4']
		z4 = data[p]['prepA_Z4']
		
		overheadSample = [calc(x1[s],z1[s],x2[s],z2[s],x3[s],z3[s],x4[s],z4[s]) for s in xrange(len(x1))]
		mean, _, err = getStats(overheadSample)		
		means.append(mean)
		errors.append(err)

	return X, means, errors


def qubitOverhead(samplesOverlap, samplesSteane4, samplesSteane12):
	# TODO: would be better to extract this information directly from the
	# prep circuits.
	qubitTimeX = 46*9 + 23
	qubitTimeZ = 46*2
	qubitTime1 = qubitTimeX
	qubitTime2 = 23*12
	qubitTime3 = qubitTimeZ
	qubitTime4 = qubitTime3
		
	_, meanOverlap, errorOverlap = overhead4(samplesOverlap, qubitTimeX, qubitTimeZ)
	_, meanSteane4, errorSteane4 = overhead4(samplesSteane4, qubitTimeX, qubitTimeZ)
	X, meanSteane12, errorSteane12 = overhead12(samplesSteane12, qubitTime1, qubitTime2, qubitTime3, qubitTime4)
			
	yList = [meanOverlap, meanSteane4, meanSteane12]
	yErrList = [errorOverlap, errorSteane4, errorSteane12]
	
	return X, yList, yErrList	

def cnotOverhead(samplesOverlap, samplesSteane4, samplesSteane12):

	cnotsX = 77*2 + 23
	cnotsZ = 23
	cnots1 = cnotsX
	cnots2 = 77 + 23
	cnots3 = 23
	cnots4 = cnots3
	
	_, meanSteane4, errorSteane4 = overhead4(samplesSteane4, cnotsX, cnotsZ)	
	_, meanSteane12, errorSteane12 = overhead12(samplesSteane12, cnots1, cnots2, cnots3, cnots4)
	
	cnotsX = 57*2 + 23
	X, meanOverlap, errorOverlap = overhead4(samplesOverlap, cnotsX, cnotsZ)
	
	yList = [meanOverlap, meanSteane4, meanSteane12]
	yErrList = [errorOverlap, errorSteane4, errorSteane12]
	
	return X, yList, yErrList

def prAccept(samplesOverlap, samplesSteane4, samplesSteane12):
	_, meanOverlap, errorOverlap = prAccept4(samplesOverlap)
	_, meanSteane4, errorSteane4 = prAccept4(samplesSteane4)
	X, meanSteane12, errorSteane12 = prAccept12(samplesSteane12)
		
	return X, [meanOverlap,meanSteane4,meanSteane12], [errorOverlap,errorSteane4,errorSteane12]


if __name__ == '__main__':
	logging.basicConfig(level=logging.INFO)
	pMin = 0
	pMax = 2e-3
	pStep = 1e-4
	iters = 100000
	
	samplesOverlap = simOverlapPrep(pMin, pMax, pStep, iters)
	samplesSteane4 = simSteane4Prep(pMin, pMax, pStep, iters)
	samplesSteane12 = simSteane12Prep(pMin, pMax, pStep, iters)
	
	samplesNROverlap = simOverlapPrep(pMin, pMax, pStep, iters, noRests=True)
	samplesNRSteane4 = simSteane4Prep(pMin, pMax, pStep, iters, noRests=True)
	samplesNRSteane12 = simSteane12Prep(pMin, pMax, pStep, iters, noRests=True)
	
	labels = ['Overlap-4', 'Steane-4', 'Steane-12']
	
	X, yList, yErrList = qubitOverhead(samplesOverlap, samplesSteane4, samplesSteane12)
	plotList(X, yList, yErrList, filename='plotQubitOverheadCompare', labelList=labels, xLabel='p', yLabel='Qubits', legendLoc='upper left')

	X, yList, yErrList = qubitOverhead(samplesNROverlap, samplesNRSteane4, samplesNRSteane12)
	plotList(X, yList, yErrList, filename='plotQubitOverheadNRCompare', labelList=labels, xLabel='p', yLabel='Qubits', legendLoc='upper left')
	
	
	X, yList, yErrList = cnotOverhead(samplesOverlap, samplesSteane4, samplesSteane12)
	plotList(X, yList, yErrList, filename='plotCnotOverheadCompare', labelList=labels, xLabel='p', yLabel='CNOTs', legendLoc='upper left')

	X, yList, yErrList = cnotOverhead(samplesNROverlap, samplesNRSteane4, samplesNRSteane12)
	plotList(X, yList, yErrList, filename='plotCnotOverheadNRCompare', labelList=labels, xLabel='p', yLabel='CNOTs', legendLoc='upper left')

	X, pr, err = prAccept(samplesOverlap, samplesSteane4, samplesSteane12)
	plotList(X, pr, err, filename='plotPrAcceptCompare', labelList=labels, xLabel='p', yLabel='Pr[accept]', legendLoc='upper right')