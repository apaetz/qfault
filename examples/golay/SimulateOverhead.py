from sim.golaySim import prepareAncilla12, prepareAncillaOverlap, \
	prepareAncillaNew, prepareAncilla12_alt, prepareAncilla12_alt2
from sim.overhead import simAncillaPrep
from util import stats, listutils
from util.cache import fetchable
from util.plotting import plotList
from util.stats import getStats
import logging
import sys

@fetchable
def simOverlapPrep(pMin, pMax, pStep, iters, noRests=False):
	return simAncillaPrep(prepareAncillaOverlap, pMin, pMax, pStep, iters, noRests)

@fetchable
def simSteane4Prep(pMin, pMax, pStep, iters, noRests=False):
	return simAncillaPrep(prepareAncillaNew, pMin, pMax, pStep, iters, noRests)

@fetchable
def simSteane12Prep(pMin, pMax, pStep, iters, noRests=False):
	return simAncillaPrep(prepareAncilla12, pMin, pMax, pStep, iters, noRests)

@fetchable
def simSteane12AltPrep(pMin, pMax, pStep, iters, noRests=False):
	return simAncillaPrep(prepareAncilla12_alt, pMin, pMax, pStep, iters, noRests)

@fetchable
def simSteane12Alt2Prep(pMin, pMax, pStep, iters, noRests=False):
	return simAncillaPrep(prepareAncilla12_alt2, pMin, pMax, pStep, iters, noRests)


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

def prAccept12_alt(data):
	
	X = sorted(data.keys())
	
	means = []
	errors = []
	for p in X:
		x0Attempts = data[p]['prepAX0']
		x1Attempts = data[p]['prepAX1']
		x2Attempts = data[p]['prepAX2']
		z0Attempts = data[p]['prepAZ0']
		z1Attempts = data[p]['prepAZ1']
		
		# For each verification stage, calculate the total number of attempts
		# and the number of successes.
		z1 = (sum(z1Attempts), iters)
		z0 = (sum(z0Attempts), z1[0])
		x2 = (sum(x2Attempts), 2*z0[0] + z1[0])
		x1 = (sum(x1Attempts), x2[0])
		x0 = (sum(x0Attempts), x1[0])		
		attSucc = [x0,x1,x2,z0,z1]
		
		samples = [[1]*succ + [0]*(att-succ) for att,succ in attSucc]
		sampleMeans = map(stats.mean, samples)
		stds = map(stats.std, samples)
		stdErrs = [stats.stdErr(stds[i], attSucc[i][0]) for i in xrange(len(stds))]
		
		err95 = stats.stdProduct(stdErrs, sampleMeans)
		
		means.append(listutils.mul(sampleMeans))
		errors.append(err95)
		
	return X, means, errors

def prAccept12_alt2(data):
	
	X = sorted(data.keys())
	
	means = []
	errors = []
	for p in X:
		z1Attempts = data[p]['prepAZ1']
		z2Attempts = data[p]['prepAZ2']
		x1Attempts = data[p]['prepAX1']
		x2Attempts = data[p]['prepAX2']
		x3Attempts = data[p]['prepAX3']
		
		# For each verification stage, calculate the total number of attempts
		# and the number of successes.
		x3 = (sum(x3Attempts), iters)
		x2 = (sum(x2Attempts), x3[0])
		x1 = (sum(x1Attempts), x2[0])
		z2 = (sum(z2Attempts), x3[0] + x2[0] + 2*x1[0])
		z1 = (sum(z1Attempts), z2[0])		
		attSucc = [z1,z2,x1,x2,x3]
		
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
	
	print 'meas:'
	print meas1,meas2,meas3,meas4
	
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
				
		if p == 0:
			print 'p=0:'
			print x1[0],x2[0],x3[0],x4[0],z1[0],z2[0],z3[0],z4[0]
		
		
		overheadSample = [calc(x1[s],z1[s],x2[s],z2[s],x3[s],z3[s],x4[s],z4[s]) for s in xrange(len(x1))]
		mean, _, err = getStats(overheadSample)
				
		means.append(mean)
		errors.append(err)

	return X, means, errors

def overhead12_alt(data, measVX0, measVX1, measVZ):

	X = sorted(data.keys())
	
	def calc(x0,x1,x2,z0,z1):
		o1 = measVX0 * x0
		o2 = measVX1 * (x1 + x2)
		o3 = measVZ * (z0 + z1)
		
		return o1 + o2 + o3
	
	means = []
	errors = []
	for p in X:
		x0 = data[p]['prepAX0']
		x1 = data[p]['prepAX1']
		x2 = data[p]['prepAX2']
		z0 = data[p]['prepAZ0']
		z1 = data[p]['prepAZ1']
		
		overheadSample = [calc(x0[s],x1[s],x2[s],z0[s],z1[s]) for s in xrange(len(x1))]
		mean, _, err = getStats(overheadSample)		
		means.append(mean)
		errors.append(err)

	return X, means, errors

def overhead12_alt2(data, measVZ0, measVZ1, measVX):

	X = sorted(data.keys())
	
	def calc(z1,z2,x1,x2,x3):
		o1 = measVZ0 * z1
		o2 = measVZ1 * z2
		o3 = measVX * (x1 + x2 + x3)
		
		return o1 + o2 + o3
	
	means = []
	errors = []
	for p in X:
		z1 = data[p]['prepAZ1']
		z2 = data[p]['prepAZ2']
		x1 = data[p]['prepAX1']
		x2 = data[p]['prepAX2']
		x3 = data[p]['prepAX3']
		
		overheadSample = [calc(z1[s],z2[s],x1[s],x2[s],x3[s]) for s in xrange(len(x1))]
		mean, _, err = getStats(overheadSample)		
		means.append(mean)
		errors.append(err)

	return X, means, errors

def qubitOverhead(samplesOverlap, samplesSteane4, samplesSteane12, samplesSteane12Alt,samplesSteane12Alt2):
	# TODO: would be better to extract this information directly from the
	# prep circuits.
	
	# TODO: 23*8 is inaccurate, since not all qubits need to be prepared before round 1.
	# (especially for the overlap circuit)
	prepTime = 23*8 - 1	# only 22 qubits are needed in round 1
	prepTimeOverlap = 23*8 - 9 # Some of the qubits can be prepared during round 1 or round 2.
	verifyTime = 46*2 + 23
	qubitTimeX = 2*prepTime + verifyTime
	qubitTimeXOverlap = 2*prepTimeOverlap + verifyTime
	qubitTimeZ = verifyTime
	qubitTime1 = qubitTimeX
	qubitTime2 = 23*8 + qubitTimeZ
	qubitTime3 = qubitTimeZ
	qubitTime4 = qubitTime3
		
	_, meanOverlap, errorOverlap = overhead4(samplesOverlap, qubitTimeXOverlap, qubitTimeZ)
	_, meanSteane4, errorSteane4 = overhead4(samplesSteane4, qubitTimeX, qubitTimeZ)
	_, meanSteane12, errorSteane12 = overhead12(samplesSteane12, qubitTime1, qubitTime2, qubitTime3, qubitTime4)
	_, meanSteane12Alt, errorSteane12Alt = overhead12_alt(samplesSteane12Alt, qubitTimeX, qubitTime2, qubitTimeZ)
	X, meanSteane12Alt2, errorSteane12Alt2 = overhead12_alt2(samplesSteane12Alt2, qubitTimeX, qubitTime2, qubitTimeZ)
			
	yList = [meanOverlap, meanSteane4, meanSteane12, meanSteane12Alt,meanSteane12Alt2]
	yErrList = [errorOverlap, errorSteane4, errorSteane12, errorSteane12Alt,errorSteane12Alt2]
	
	return X, yList, yErrList	

def cnotOverhead(samplesOverlap, samplesSteane4, samplesSteane12, samplesSteane12Alt,samplesSteane12Alt2):

	cnotsX = 77*2 + 23
	cnotsZ = 23
	cnots1 = cnotsX
	cnots2 = 77 + 23
	cnots3 = 23
	cnots4 = cnots3
	
	_, meanSteane4, errorSteane4 = overhead4(samplesSteane4, cnotsX, cnotsZ)	
	_, meanSteane12, errorSteane12 = overhead12(samplesSteane12, cnots1, cnots2, cnots3, cnots4)
	_, meanSteane12Alt, errorSteane12Alt = overhead12_alt(samplesSteane12Alt, cnotsX, cnots2, cnotsZ)
	_, meanSteane12Alt2, errorSteane12Alt2 = overhead12_alt2(samplesSteane12Alt2, cnotsX, cnots2, cnotsZ)
	
	cnotsX = 57*2 + 23
	X, meanOverlap, errorOverlap = overhead4(samplesOverlap, cnotsX, cnotsZ)
	
	yList = [meanOverlap, meanSteane4, meanSteane12, meanSteane12Alt, meanSteane12Alt2]
	yErrList = [errorOverlap, errorSteane4, errorSteane12, errorSteane12Alt, errorSteane12Alt2]
	
	return X, yList, yErrList

def prAccept(samplesOverlap, samplesSteane4, samplesSteane12, samplesSteane12Alt,samplesSteane12Alt2):
	_, meanOverlap, errorOverlap = prAccept4(samplesOverlap)
	_, meanSteane4, errorSteane4 = prAccept4(samplesSteane4)
	_, meanSteane12, errorSteane12 = prAccept12(samplesSteane12)
	_, meanSteane12Alt, errorSteane12Alt = prAccept12_alt(samplesSteane12Alt)
	X, meanSteane12Alt2, errorSteane12Alt2 = prAccept12_alt2(samplesSteane12Alt2)
		
	return X, [meanOverlap,meanSteane4,meanSteane12,meanSteane12Alt,meanSteane12Alt2], [errorOverlap,errorSteane4,errorSteane12,errorSteane12Alt,errorSteane12Alt2]


if __name__ == '__main__':
	logging.basicConfig(level=logging.INFO)
	pMin = 0
	pMax = 2e-3
	pStep = 1e-4
	iters = 100000
	
	if len(sys.argv) > 1:
		from counting import countParallel
		nslots = int(sys.argv[1])
		countParallel.configureMultiProcess(nslots)
	
	samplesOverlap = simOverlapPrep(pMin, pMax, pStep, iters)
	samplesSteane4 = simSteane4Prep(pMin, pMax, pStep, iters)
	samplesSteane12 = simSteane12Prep(pMin, pMax, pStep, iters)
	samplesSteane12_alt = simSteane12AltPrep(pMin, pMax, pStep, iters)
	samplesSteane12_alt2 = simSteane12Alt2Prep(pMin, pMax, pStep, iters)
	
	samplesNROverlap = simOverlapPrep(pMin, pMax, pStep, iters, noRests=True)
	samplesNRSteane4 = simSteane4Prep(pMin, pMax, pStep, iters, noRests=True)
	samplesNRSteane12 = simSteane12Prep(pMin, pMax, pStep, iters, noRests=True)
	samplesNRSteane12_alt = simSteane12AltPrep(pMin, pMax, pStep, iters, noRests=True)
	samplesNRSteane12_alt2 = simSteane12Alt2Prep(pMin, pMax, pStep, iters, noRests=True)
	
	labels = ['Overlap-4', 'Steane-4', 'Steane-12', 'Steane-12 (alt)', 'Steane-12 (alt2)']
	
	X, yList, yErrList = qubitOverhead(samplesOverlap, samplesSteane4, samplesSteane12, samplesSteane12_alt, samplesSteane12_alt2)
	plotList(X, yList[:3], yErrList[:3], filename='plotQubitOverheadCompare', labelList=labels, xLabel='p', yLabel='Qubits', legendLoc='upper left',
			xlim=[0,8000])

	X, yList, yErrList = qubitOverhead(samplesNROverlap, samplesNRSteane4, samplesNRSteane12, samplesNRSteane12_alt, samplesSteane12_alt2)
	plotList(X, yList, yErrList, filename='plotQubitOverheadNRCompare', labelList=labels, xLabel='p', yLabel='Qubits', legendLoc='upper left')
	
	
	X, yList, yErrList = cnotOverhead(samplesOverlap, samplesSteane4, samplesSteane12, samplesSteane12_alt, samplesSteane12_alt2)
	plotList(X, yList[:5], yErrList[:5], filename='plotCnotOverheadCompare', labelList=labels, xLabel='p', yLabel='CNOTs', legendLoc='upper left')

	X, yList, yErrList = cnotOverhead(samplesNROverlap, samplesNRSteane4, samplesNRSteane12, samplesNRSteane12_alt, samplesSteane12_alt2)
	plotList(X, yList, yErrList, filename='plotCnotOverheadNRCompare', labelList=labels, xLabel='p', yLabel='CNOTs', legendLoc='upper left')

	X, pr, err = prAccept(samplesOverlap, samplesSteane4, samplesSteane12, samplesSteane12_alt,samplesSteane12_alt2)
	plotList(X, pr, err, filename='plotPrAcceptCompare', labelList=labels, xLabel='p', yLabel='Pr[accept]', legendLoc='upper right')