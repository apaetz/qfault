'''
Created on 2011-04-12

@author: adam
'''
from component.xverify import countXVerifyXOnly
from component.zverify import countZVerifyXOnly_uncorrected
from counting.levelOne import countExRecs, transformedWeights
from util.plotting import plotPolyList
import GolayCounting
import math
	
	
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
	
	filename = 'plot-level1-tranformedX-' + ratiosXStr + settingsStr
	valuesXTransformed = [GammaX * ratiosX[label] for label in labelsX]
	plotPolyList(valuesXTransformed, 
				 gMin, gMax, 
				 filename, 
				 labelList=labelsX, 
				 numPoints=20,
				 xLabel=r'$\gamma$',
				 legendLoc='upper left')

		
	filename = 'plot-level1-maligX-' + settingsStr
	plotPolyList(valuesX + [GammaX], 
				 gMin, gMax, 
				 filename, 
				 labelList=labelsX + [r'$\Gamma_X$'], 
				 numPoints=20,
				 xLabel=r'$\gamma$',
				 legendLoc='upper left')
	

	filename = 'plot-level1-maligZ-' + settingsStr
	plotPolyList(valuesZ + [GammaZ], 
				 gMin, gMax, 
				 filename, 
				 labelList=labelsZ + [r'$\Gamma_Z$'], 
				 numPoints=20,
				 xLabel=r'$\gamma$',
				 legendLoc='upper left')

	filename = 'plot-level1-tranformedZ-' + ratiosZStr + settingsStr
	valuesZTransformed = [GammaZ * ratiosZ[label] for label in labelsZ]
	plotPolyList(valuesZTransformed, 
				 gMin, gMax, 
				 filename, 
				 labelList=labelsZ, 
				 numPoints=20,
				 xLabel=r'$\gamma$',
				 legendLoc='upper left')

	
	#filename = 'plot-ratios-' + settingsStr
	#ratioPolys = [Gamma * r for r in ratios]
	#plotPolyList(ratioPolys, gMin, gMax, filename, results.keys(), numPoints=20)	
	
def plotPrAccept(ancillaPreps, settings):
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
	
	plots = [resultX1.prAccept, resultX2.prAccept, resultZ.prAccept]
	labels = [r'$\Pr[X_1]$', 
			  r'$\Pr[X_2]$', 
			  r'$\Pr[Z | X_1, X_2]$']
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



if __name__ == '__main__':
	settings = GolayCounting.globalSettings
	preps = GolayCounting.getSteaneRandomPreps()
	#preps = GolayCounting.getOverlapPreps()
	pairs = ((preps[0], preps[1]), (preps[2], preps[3]))
	
	#GammaX, weightsX, GammaZ, weightsZ, gMax = transformedWeights(pairs, settings)
	#settings = getTransformedSettings(settings, GammaX, weightsX, GammaZ, weightsZ, gMax)
	
	#configureMultiProcess(12)
	
	#plotCnotExRecDetails(pairs, settings, gMaxAlt=0.000002)
	plotLevel1Malig(pairs, settings)
	
	preps = GolayCounting.getSteaneRandomPreps()
	preps = [prep.filterAgainst('rest') for prep in preps]
	#preps = GolayCounting.getOverlapPreps()
	pairs = ((preps[0], preps[1]), (preps[2], preps[3]))
		
	#plotCnotExRecDetails(pairs, settings, gMaxAlt=0.000002)
	plotLevel1Malig(pairs, settings)