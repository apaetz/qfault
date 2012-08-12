'''
Created on 2010-08-20

@author' :  adam
'''
from qfault.adt.tree import TreeNode
from copy import copy
from qfault.noise import NoiseModelZeroSympy, TransformedNoiseModelXSympy, \
	TransformedNoiseModelZSympy
import logging
import math
from qfault.counting import levelOne

logger = logging.getLogger('settings')

class PrettyFloat(float):
	'''
	For shorter filenames.
	'''
	
	def __str__(self):
		return '{0:g}'.format(self)
	

	
class ComponentSettings(TreeNode):
	
	def getSubcomponent(self, name):
		return self.child(name)
	
	def __repr__(self):
		return str(self)

def buildSettingsTree(settings):
	# Build the tree from the bottom up.
	noise = settings['noise']
	zero = ComponentSettings('zero', {'kGood': settings['kGood0']})
	cnotMeas = ComponentSettings('cm', {'kGood': settings['kGoodVX']}) 
	verifyX = ComponentSettings('vx', 
							    {'kGood': settings['kGoodVXTotal'],
								 'kBest': settings['kBestVX']},
								children=[zero,cnotMeas])
	
	restCnotMeasVerify = ComponentSettings('rcm', {'kGood': settings['kGoodVZ']})
	verifyZ = ComponentSettings('vz',
							    {'kBest': settings['kBestVZ'],
								 'kGood': settings['kGoodECTotal']},
								children=[verifyX,restCnotMeasVerify])

	restCnotMeasEC = ComponentSettings('rcm', {'kGood': settings['kGoodECX']})	
	ec = ComponentSettings('ec',
						   {'kGood': settings['kGoodECTotal']},
						   children=[verifyZ, restCnotMeasEC])
	
	cnot = ComponentSettings('cnot', {'kGood': settings['kGoodExRecCnot']})
	
	exRec = ComponentSettings('exRec',
							  {'kGood': settings['kGoodExRec'],
							   'kGood-LEC-CNOT': settings['kGoodLECcnot'],
							   'pMin': settings['pMin'],
							   'pMax': settings['pMax'],
							   'noise': noise},
							  children=[ec, cnot])
	
	return exRec





def getTransformedSettings(settings, GammaX, weightsX, GammaZ, weightsZ, gMax):
	
		for key in weightsX:
			weightsX[key] = int(math.ceil(weightsX[key]))
		for key in weightsZ:
			weightsZ[key] = int(math.ceil(weightsZ[key]))
			
		labels = levelOne.maligEventLabels
		
		prepZ = weightsX[labels['prep0']]
		measZ = weightsX[labels['measZ']]
		rest0 = weightsX[labels['rest0']]
		cnotIX = weightsX[labels['IX']]
		cnotXI = weightsX[labels['XI']]
		cnotXX = weightsX[labels['XX']]
		
		prepX = weightsZ[labels['prep+']]
		measX = weightsZ[labels['measX']]
		restPlus = weightsZ[labels['rest+']]
		cnotIZ = weightsZ[labels['IZ']]
		cnotZI = weightsZ[labels['ZI']]
		cnotZZ = weightsZ[labels['ZZ']]
		
		
		# Adjust the noise settings.
		settings = copy(settings)
		gMin, _ = settings['noise']['X'].noiseRange()
		GMinX = GammaX(gMin)
		GMaxX = GammaX(gMax)
		GMinZ = GammaZ(gMin)
		GMaxZ = GammaZ(gMax)
		logger.info('Gmin=%f,%f : Gmax=%f,%f', GMinX, GMinZ, GMaxX, GMaxZ)
		
		noiseDict = {'XZ': NoiseModelZeroSympy(),
					 'X':  TransformedNoiseModelXSympy(prepZ, measZ, rest0, cnotIX, cnotXI, cnotXX, GMinX, GMaxX),
					 'Z':  TransformedNoiseModelZSympy(prepX, measX, restPlus, cnotIZ, cnotZI, cnotZZ, GMinZ, GMaxZ)}
		settings['noise'] = noiseDict
		settingsEC = settings.getSubcomponent('ec')
		settingsVZ = settingsEC.getSubcomponent('vz')
		settingsVX = settingsVZ.getSubcomponent('vx')
		
		# There is no XZ correlation information.
		settingsVX['kBest'] = 0
		settingsVZ['kBest'] = 0
		
		return settings	