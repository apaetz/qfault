'''
Created on 2011-04-12

@author: adam
'''
from component.exrec import exRecXOnly, countCnotExRec, exRecZOnly
from component.xverify import countXVerifyXOnly
from component.zverify import countZVerifyXOnly_uncorrected
from counting.levelOne import countExRecs, transformedWeights
from counting.probability import countResultAsPoly, countsAsProbability
from util.plotting import plotPolyList, Series, barSeries
import GolayCounting
import math
from component import xverify, zverify
from counting.countErrors import CountResult
from util import listutils

# Malignant event labels defined in levelOne.py do not necessarily match
# the labels in the paper.  This maps between the two.
maligEventLabelMap = {
	r'mal$_{IX}$': r'mal$_{IX}$',
	r'mal$_{XI}$': r'mal$_{XI}$',
	r'mal$_{XX}$': r'mal$_{XX}$',
	r'mal$_{IZ}$': r'mal$_{IZ}$', 
	r'mal$_{ZI}$': r'mal$_{ZI}$',
	r'mal$_{ZZ}$': r'mal$_{ZZ}$',
	r'mal$_{\mathrm{rest},|0\rangle}$': r'$\mathrm{mal}_{X}^{\mathrm{rest}}$',
	r'mal$_{|0\rangle}$': r'$\mathrm{mal}_{X}^{\mathrm{prep}}$',
	r'mal$_{\mathrm{meas},Z}$': r'$\mathrm{mal}_{X}^{\mathrm{meas}}$',
	r'mal$_{\mathrm{rest},|+\rangle}$': r'$\mathrm{mal}_{Z}^{\mathrm{rest}}$',
	r'mal$_{|+\rangle}$': r'$\mathrm{mal}_{Z}^{\mathrm{prep}}$',
	r'mal$_{\mathrm{meas},X}$': r'$\mathrm{mal}_{Z}^{\mathrm{meas}}$'
}
	
def mapLabels(labels):
	return [maligEventLabelMap[label] for label in labels]
	
def plotLevel1Malig(ancillaPairs, settings):
	#noise = settings['noise']

	resultsX, resultsZ = countExRecs(ancillaPairs, settings)
	
	# Sort, for nice looking legend.
	labelsX = sorted(resultsX.keys())
	valuesX = [resultsX[label] for label in labelsX]
	labelsZ = sorted(resultsZ.keys())
	valuesZ = [resultsZ[label] for label in labelsZ]
		
	GammaX, ratiosX, GammaZ, ratiosZ, gMax = transformedWeights(ancillaPairs, settings)

	for key in ratiosX:
		ratiosX[key] = int(math.ceil(ratiosX[key]))
	for key in ratiosZ:
		ratiosZ[key] = int(math.ceil(ratiosZ[key]))
		
	print 'ratiosX=', ratiosX
	print 'ratiosZ=', ratiosZ

	settingsStr = str(ancillaPairs) + str(settings)
	settingsStr = settingsStr.replace(' ', '').replace('.', '').replace(')','').replace('(','')
	settingsStr = settingsStr.replace('}', ']').replace('{', '[')		
		

	gMin = gMax/500
	
	ratiosXStr = str(ratiosX.values()).replace(' ', '').replace('.', '').replace(')','').replace('(','')
	ratiosZStr = str(ratiosZ.values()).replace(' ', '').replace('.', '').replace(')','').replace('(','')	
	
	# For comparison, we want to plot X and Z events using the same y-axis limits.	
	yMax = max(p(gMax) for p in valuesX+valuesZ) * 1.05
	# pylot chokes on mpf (bug?), so cast here.
	yMax = float(yMax)
	ylim = (0,yMax)

	numPoints = 20
	
	filename = 'plot-level1-maligX-' + settingsStr
	plotPolyList(valuesX,# + [GammaX], 
				 gMin, gMax, 
				 filename, 
				 labelList=mapLabels(labelsX) + [r'$\Gamma_X$'], 
				 numPoints=numPoints,
				 xLabel=r'$\gamma$',
				 legendLoc='upper left',
				 ylim = ylim)
	

	filename = 'plot-level1-maligZ-' + settingsStr
	plotPolyList(valuesZ,# + [GammaZ], 
				 gMin, gMax, 
				 filename, 
				 labelList=mapLabels(labelsZ) + [r'$\Gamma_Z$'], 
				 numPoints=numPoints,
				 xLabel=r'$\gamma$',
				 legendLoc='upper left',
				 ylim = ylim)
	
	
	filename = 'plot-level1-tranformedX-' + ratiosXStr + settingsStr
	valuesXTransformed = [GammaX * ratiosX[label] for label in labelsX]
	plotPolyList(valuesXTransformed, 
				 gMin, gMax, 
				 filename, 
				 labelList=mapLabels(labelsX), 
				 numPoints=numPoints,
				 xLabel=r'$\gamma$',
				 legendLoc='upper left')

	filename = 'plot-level1-tranformedZ-' + ratiosZStr + settingsStr
	valuesZTransformed = [GammaZ * ratiosZ[label] for label in labelsZ]
	plotPolyList(valuesZTransformed, 
				 gMin, gMax, 
				 filename, 
				 labelList=mapLabels(labelsZ), 
				 numPoints=numPoints,
				 xLabel=r'$\gamma$',
				 legendLoc='upper left')

	
	#filename = 'plot-ratios-' + settingsStr
	#ratioPolys = [Gamma * r for r in ratios]
	#plotPolyList(ratioPolys, gMin, gMax, filename, results.keys(), numPoints=20)	
	
def plotXZCorrection(ancillaPreps, settings):
	pair1, pair2 = ancillaPreps
	noise = settings['noise']
	settingsVX = settings.getSubcomponent('ec').getSubcomponent('vz').getSubcomponent('vx')
	
	zResult = xverify.countXVerifyZOnly_uncorrected(pair1[0], pair1[1], settingsVX, noise)
	corrections = xverify.getXZCorrections(pair1[0], pair1[1], settingsVX, noise, zResult.locTotals)
	
	# We want to plot corrections by fault-order.  Syndrome corrections are not important.
	uncorrectedZero = [c[0] for c in zResult.counts]
	uncorrectedNZ = [sum(c)-c[0] for c in zResult.counts]
	correctedZero = [c[0] for c in corrections]
	correctedNZ = [sum(c)-c[0] for c in corrections]
		
	corrResults = [CountResult([0]*(k) + [correctedNZ[k]], None, zResult.locTotals, 1, 0) for k in range(len(correctedNZ))]
	corrResultsZero = [CountResult([0]*(k) + [correctedZero[k]], None, zResult.locTotals, 1, 0) for k in range(len(correctedZero))]
	uncorrResults = [CountResult([0]*(k) + [uncorrectedNZ[k]], None, zResult.locTotals, 1, 0) for k in range(len(uncorrectedNZ))]
	uncorrResultsZero = [CountResult([0]*(k) + [uncorrectedZero[k]], None, zResult.locTotals, 1, 0) for k in range(len(uncorrectedZero))]
	
	corrPolys = [countResultAsPoly(r, noise['Z']) for r in corrResults]
	corrPolysZero = [countResultAsPoly(r, noise['Z']) for r in corrResultsZero]
	uncorrPolys = [countResultAsPoly(r, noise['Z']) for r in uncorrResults]
	uncorrPolysZero = [countResultAsPoly(r, noise['Z']) for r in uncorrResultsZero]

	g = .001/15
	corrBar = [-pr(g) for pr in corrPolys]
	corrBarZero = [-pr(g) for pr in corrPolysZero]
	uncorrBar = [pr(g) for pr in uncorrPolys]
	uncorrBarZero = [pr(g) for pr in uncorrPolysZero]
	
	corrected = listutils.addLists(corrPolys)#, corrPolysZero)
	uncorrected = listutils.addLists(uncorrPolys)#, uncorrPolysZero)
	
	print 'g=', g
	gUncorr = [pr(g) for pr in uncorrected]
	gCorr = [pr(g) for pr in corrected]
	print 'uncorrected=', gUncorr 
	print 'corrected=', gCorr
	print 'ratio=', sum(gCorr) / sum(gUncorr)
	
	corrSeries = Series(range(len(corrBar)), corrBar, 'correction')
	corrSeriesZero = Series(range(len(corrBarZero)), corrBarZero, 'correction (0)')

	uncorrSeries = Series(range(len(uncorrBar)), uncorrBar, 'uncorrected')
	uncorrSeriesZero = Series(range(len(uncorrBarZero)), uncorrBarZero, 'uncorrected (0)')
	
	barSeries([uncorrSeries, corrSeries, uncorrSeriesZero, corrSeriesZero], filename='bar-XZ-correction', colors=['g','r', 'b', 'k'])

	
	gMin = 0
	_, gMax = noise['Z'].noiseRange()
	labels = ['k=1', 'k=2', 'k=3'] + ['0','1', '2', '3', '4','5','6','7','8']
	ylabel = r'Pr[K=k,best,$\neg$accept]'
	plotPolyList(corrPolys[1:] + uncorrPolys, gMin, gMax, 'plot-XZ-corrections-', labelList=labels, xLabel=r'$\gamma$', yLabel=ylabel, numPoints=20, legendLoc='upper left')


def plotXZCorrection_X(ancillaPreps, settings):
	pair1, pair2 = ancillaPreps
	noise = settings['noise']
	settingsVZ = settings.getSubcomponent('ec').getSubcomponent('vz')
	
	xResult = zverify.countZVerifyXOnly_uncorrected(pair1, pair2, settingsVZ, noise)
	corrections = zverify.getXZCorrections(pair1, pair2, settingsVZ, noise, xResult.locTotals)
	
	# We want to plot corrections by fault-order.  Syndrome corrections are not important.
	uncorrectedZero = [c[0] for c in xResult.counts]
	uncorrectedNZ = [sum(c)-c[0] for c in xResult.counts]
	correctedZero = [c[0] for c in corrections]
	correctedNZ = [sum(c)-c[0] for c in corrections]
		
	corrResults = [CountResult([0]*(k) + [correctedNZ[k]], None, xResult.locTotals, 1, 0) for k in range(len(correctedNZ))]
	corrResultsZero = [CountResult([0]*(k) + [correctedZero[k]], None, xResult.locTotals, 1, 0) for k in range(len(correctedZero))]
	uncorrResults = [CountResult([0]*(k) + [uncorrectedNZ[k]], None, xResult.locTotals, 1, 0) for k in range(len(uncorrectedNZ))]
	uncorrResultsZero = [CountResult([0]*(k) + [uncorrectedZero[k]], None, xResult.locTotals, 1, 0) for k in range(len(uncorrectedZero))]
	
	corrPolys = [countResultAsPoly(r, noise['Z']) for r in corrResults]
	corrPolysZero = [countResultAsPoly(r, noise['Z']) for r in corrResultsZero]
	uncorrPolys = [countResultAsPoly(r, noise['Z']) for r in uncorrResults]
	uncorrPolysZero = [countResultAsPoly(r, noise['Z']) for r in uncorrResultsZero]

	g = .001/15
	
	corrected = listutils.addLists(corrPolys)#, corrPolysZero)
	uncorrected = listutils.addLists(uncorrPolys)#, uncorrPolysZero)
	
	print 'X-error corrections for Z-verification:'
	print 'g=', g
	gUncorr = [pr(g) for pr in uncorrected]
	gCorr = [pr(g) for pr in corrected]
	print 'uncorrected=', gUncorr 
	print 'corrected=', gCorr
	print 'ratio=', sum(gCorr) / sum(gUncorr)	
	
	
def plotPrAccept(ancillaPreps, settings):
	pair1, pair2 = ancillaPreps

	noise = settings['noise']
	gMin = 0
	_, gMax = noise['X'].noiseRange()
	settingsEC = settings.getSubcomponent('ec')
	settingsVZ = settingsEC.getSubcomponent('vz')
	settingsVX = settingsVZ.getSubcomponent('vx')
#	
	resultX1 = countXVerifyXOnly(pair1[0], pair1[1], settingsVX, noise)	
	resultX2 = countXVerifyXOnly(pair2[0], pair2[1], settingsVX, noise)
	resultZ = countZVerifyXOnly_uncorrected(pair1, pair2, settingsVZ, noise)
	
	plots = [resultX1.prAccept, resultX2.prAccept, resultZ.prAccept]
	labels = [r'Pr[accept$^{(1)}$]', 
			  r'Pr[accept$^{(2)}$]', 
			  r'Pr[accept | accept$^{(1)}$, accept$^{(2)}$]']
	plotPolyList(plots, gMin, gMax, 'plot-prAccept', labelList=labels, xLabel=r'$\gamma$', numPoints=20)
	

def plotOverhead(ancillaPreps, settings):
	pair1, pair2 = ancillaPreps

	noise = settings['noise']
	gMin, gMax = noise['X'].noiseRange()
	settingsEC = settings.getSubcomponent('ec')
	settingsVZ = settingsEC.getSubcomponent('vz')
	settingsVX = settingsVZ.getSubcomponent('vx')
#	
	resultX1 = countXVerifyXOnly(pair1[0], pair1[1], settingsVX, noise)	
	resultX2 = countXVerifyXOnly(pair2[0], pair2[1], settingsVX, noise)
	resultZ = countZVerifyXOnly_uncorrected(pair1, pair2, settingsVZ, noise)
		
	# Overhead = qubits * time / Pr[accept]
	qubitsTimeX = 9 * 46 + 23
	qubitsTimeZ = 46*2
	overhead = qubitsTimeX / resultX1.prAccept + \
			   qubitsTimeX / resultX2.prAccept + \
			   qubitsTimeZ / resultZ.prAccept

	plotPolyList([overhead], gMin, gMax, 'overheadPlot', xLabel=r'$\gamma$', yLabel='qubit overhead', numPoints=10)



def plotCnotExRecDetails(ancillaPairs, settings, gMaxAlt=1.):
	noise = settings['noise']
	gMin = 0
	_, gMax = noise['X'].noiseRange()
	gMax = min(gMaxAlt, gMax)
	
	zOnly = exRecZOnly(ancillaPairs, settings)	
	xOnly = exRecXOnly(ancillaPairs, settings)
	xMaxes, zMaxes = countCnotExRec(ancillaPairs, settings)
	
	badChars = ['.', ',', ' ', ')', '(', '}', '{']
	pairStr = str(ancillaPairs)
	for char in badChars:
		pairStr = pairStr.replace(char, '')
	
	errorStrX = ['IX', 'XI', 'XX']
	errorStrZ = ['IZ', 'ZI', 'ZZ']
	ecLabels = ['--', '-B', 'A-', 'AB', 'max']
	
	def plotDetail(polys, filename, error, ylabel='', yscale='linear', legend='upper left'):
		plotPolyList(polys, gMin, gMax, filename, labelList=ecLabels, numPoints=20, legendLoc=legend, xLabel=r'$\gamma$', yLabel=ylabel, yscale=yscale)
		
	for error in reversed(range(3)):
		esX = errorStrX[error]
		prBadXList = [result.prBad for result in xOnly[error]]
		plotDetail(prBadXList, 'plot-prBad-cnot-' + esX + pairStr, error, ylabel='Pr[bad]', yscale='log')
		polysX = [countResultAsPoly(r, noise['X']) for r in xOnly[error]] + [xMaxes[error]]
		plotDetail(polysX, 'plot-cnot-exrec-malig-' + esX + pairStr, error, ylabel=r'Pr[mal$_{' + esX + r'}$]')
		plotDetail(polysX, 'plot-log-cnot-exrec-malig-' + esX + pairStr, error, ylabel=r'Pr[mal$_{' + esX + r'}$]', yscale='log')

		#print 'Pr[bad](0.00014)=', [p(0.00014) for p in prBadXList]
		#prAcceptXList = [result.prAccept for result in  xOnly[error]]
		#plotDetail(prAcceptXList, 'plot-prAccept-cnot-' + errorStrX[error] + pairStr, error, legend='upper right')
		
		esZ = errorStrZ[error]
		polysZ = [countResultAsPoly(r, noise['Z']) for r in zOnly[error]] + [zMaxes[error]]
		plotDetail(polysZ, 'plot-cnot-exrec-malig-' + esZ + pairStr, error, ylabel=r'Pr[mal$_{' + esZ + r'}$]')
		plotDetail(polysZ, 'plot-log-cnot-exrec-malig-' + esZ + pairStr, error, ylabel=r'Pr[mal$_{' + esZ + r'}$]', yscale='log')
		prBadZList = [result.prBad for result in zOnly[error]]
		plotDetail(prBadZList, 'plot-prBad-cnot-' + esZ + pairStr, error, ylabel='Pr[bad]', yscale='log')
		#prAcceptZList = [result.prAccept for result in  xOnly[error]]
		#plotDetail(prAcceptZList, 'plot-prAccept-cnot-' + errorStrZ[error] + pairStr, error, legend='upper right')


def xverify_errorCount(preps, settings):
	from counting.countErrors import countXerrorsZero, propagateAndReduceZeroX
	pair1, pair2 = preps
	noise = settings['noise']
	settings0 = settings.getSubcomponent('ec').getSubcomponent('vz').getSubcomponent('vx').getSubcomponent('zero')

	prepA = propagateAndReduceZeroX(pair1[0])
	countsA1 = [countXerrorsZero(k, prepA, 'A', noise) for k in range(settings0['kGood']+1)]
	#countsA2 = [countXerrorsZero(k, pair1[1], 'B', noise) for k in range(settingsVX['kGood']+1)]

	numErrors = [listutils.numNonzeroEntries(c) for c in countsA1]
	print 'numErrors=', numErrors

if __name__ == '__main__':
	settings = GolayCounting.globalSettings
	#preps = GolayCounting.getSteaneRandomPreps()
	preps = GolayCounting.getOverlapPreps()
	pairs = ((preps[0], preps[1]), (preps[2], preps[3]))
	
	#GammaX, weightsX, GammaZ, weightsZ, gMax = transformedWeights(pairs, settings)
	#settings = getTransformedSettings(settings, GammaX, weightsX, GammaZ, weightsZ, gMax)
	
	#configureMultiProcess(12)
	
	#xverify_errorCount(pairs, settings)
	#plotXZCorrection(pairs, settings)
	#plotXZCorrection_X(pairs, settings)
	plotPrAccept(pairs, settings)
	plotCnotExRecDetails(pairs, settings)
	plotLevel1Malig(pairs, settings)
	
	#preps = [prep.filterAgainst('rest') for prep in preps]
	#preps = GolayCounting.getOverlapPreps()
	#pairs = ((preps[0], preps[1]), (preps[2], preps[3]))
		
	#plotCnotExRecDetails(pairs, settings, gMaxAlt=0.000002)
	#plotLevel1Malig(pairs, settings)
	pass