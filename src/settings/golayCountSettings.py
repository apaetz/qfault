'''
Created on 2010-08-20

@author' :  adam
'''
from noise import NoiseModelXSympy, NoiseModelXZSympy
from settings.countSettings import PrettyFloat, buildSettingsTree
import gmpy
import logging

logger = logging.getLogger('settings')
	
# Note, pMin is obsolete.
pMin = PrettyFloat(0.0009)
pMax = PrettyFloat(0.0023)


# These need to be multi-precision because powers of
# the binomial failure/ideal probabilities can have extremely large
# coefficients.  So large, that they cannot be converted into double
# types.
gMin = gmpy.mpf(pMin/15)
gMax = gmpy.mpf(pMax/15)

noiseDict = {'XZ': NoiseModelXZSympy(gMin, gMax),
			 'X':  NoiseModelXSympy(gMin, gMax),
			 'Z':  NoiseModelXSympy(gMin, gMax)
			}

countSettings =  {
	'kGood0' : 4,
	'kGoodVX' : 4,
	'kGoodVXTotal': 6,
	'kBestVX': 3,
	'kGoodVZ' : 4,
	'kBestVZ': 1,
	'kGoodECX' : 3,
	'kGoodECZ' : 3,
	'kGoodECTotal' : 11,
	'kGoodExRecCnot' : 2,
	'kGoodLECcnot' : (3,3,1), # LEC-A, LEC-B, CNOT
	'kGoodExRec' : 25,
	'noise': noiseDict,
	'pMin': pMin,
	'pMax': pMax
}

testSettings = {
	'kGood0' : 1,
	'kGoodVX' : 1,
	'kGoodVXTotal': 1,
	'kBestVX': 1,
	'kGoodVZ' : 1,
	'kBestVZ': 1,
	'kGoodECX' : 1,
	'kGoodECZ' : 1,
	'kGoodECTotal' : 2,
	'kGoodExRecCnot' : 1,
	'kGoodExRec' : 4,
	'kGoodLECcnot' : (2,2,1), # LEC-A, LEC-B, CNOT
	'noise': noiseDict,
	'pMin': pMin,
	'pMax': pMax
}

testSettingsMedium = {
	'kGood0' : 2,
	'kGoodVX' : 2,
	'kGoodVXTotal': 4,
	'kBestVX': 1,
	'kGoodVZ' : 2,
	'kBestVZ': 1,
	'kGoodECX' : 2,
	'kGoodECZ' : 2,
	'kGoodECTotal' : 6,
	'kGoodExRecCnot' : 1,
	'kGoodExRec' : 9,
	'noise': noiseDict,
	'pMin': pMin,
	'pMax': pMax
}

def getSettings():
	return buildSettingsTree(countSettings)

def getTestSettings():
	return buildSettingsTree(testSettings)

def getMediumSettings():
	return buildSettingsTree(testSettingsMedium)

if __name__ == '__main__':
	settings = buildSettingsTree(countSettings)
	print settings
	