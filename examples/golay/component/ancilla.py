'''
 ancilla.py

 This file contains functions for counting X and Z errors on an encoded |0> or |+> ancilla
 that has been verified against X and Z errors.
 
 These functions act primarily as wrappers around the Z-error verification component
 functions in zverify.py.  The difference is that these functions set the resulting
 probabilities assuming that all verification stages have accepted (i.e.,  
 prBad = Pr[bad | accept] and prAccept = Pr[X_1,X_2,Z]).
'''

from util.counterUtils import isFaultTolerant
from component.xverify import countXVerifyXOnly, countXVerifyZOnly_uncorrected
from component.zverify import countZVerifyXOnly, countZVerifyZOnly
from counting.countErrors import CountResult
from golay import golayCode

corrector = golayCode.Corrector()

def countVerifiedZeroXOnly(prepPairA, prepPairB, settings, noise):
	'''
	Computes the weighted X-error counts for a verified encoded |0>.
	
	Only X errors are considered by this function.  To count
	Z errors see countVerifiedZOnly.
	
	@param prepPairA: The encoding circuit pair for X-error verification on block A.
	@type prepPairA:  tuple  - (Locations, Locations)
	@param prepPairB: The encoding circuit pair for X-error verification on block B.
	@type prepPairB:  tuple  - (Locations, Locations)
	@param settings:  Z-error verification settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''

	settingsVX = settings.getSubcomponent('vx')
	resultX1 = countXVerifyXOnly(prepPairA[0], prepPairA[1], settingsVX, noise)
	resultX2 = countXVerifyXOnly(prepPairB[0], prepPairB[1], settingsVX, noise)	
	resultZ = countZVerifyXOnly(prepPairA, prepPairB, settings, noise)

	# Pr[bad | X1, X2, Z] <= Pr[bad | X1,X2] / Pr[Z]
	prBad = resultZ.prBad / resultZ.prAccept

	# Pr[X1, X2, Z] >= Pr[X1]*Pr[X2]*Pr[Z]
	prAccept = resultX1.prAccept * resultX2.prAccept * resultZ.prAccept

	
	# Perform some sanity checks.
	if not isFaultTolerant(resultZ.counts, 3, corrector):
		raise Exception
	
	return CountResult(resultZ.counts, None, resultZ.locTotals, prAccept, prBad)

def countVerifiedZeroZOnly(prepPairA, prepPairB, settings, noise):
	'''
	Computes the weighted Z-error counts for a verified encoded |0>.
	
	Only Z errors are considered by this function.  To count
	X errors see countVerifiedZOnly.
	
	@param prepPairA: The encoding circuit pair for X-error verification on block A.
	@type prepPairA:  tuple  - (Locations, Locations)
	@param prepPairB: The encoding circuit pair for X-error verification on block B.
	@type prepPairB:  tuple  - (Locations, Locations)
	@param settings:  Z-error verification settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''

	settingsVX = settings.getSubcomponent('vx')
	resultX1 = countXVerifyZOnly_uncorrected(prepPairA[0], prepPairA[1], settingsVX, noise)
	resultX2 = countXVerifyZOnly_uncorrected(prepPairB[0], prepPairB[1], settingsVX, noise)	
	resultZ = countZVerifyZOnly(prepPairA, prepPairB, settings, noise)

	# Pr[bad | X1, X2, Z] <= Pr[bad | X1,X2] / Pr[Z]
	prBad = resultZ.prBad / resultZ.prAccept

	# Pr[X1, X2, Z] >= Pr[X1]*Pr[X2]*Pr[Z | X1, X2]
	prAccept = resultX1.prAccept * resultX2.prAccept * resultZ.prAccept

	
	# Perform some sanity checks.
	if not isFaultTolerant(resultZ.counts, 3, corrector):
		raise Exception
	
	return CountResult(resultZ.counts, None, resultZ.locTotals, prAccept, prBad)

def countVerifiedPlusXOnly(prepPairA, prepPairB, settings, noise):
	'''
	Computes the weighted X-error counts for a verified encoded |+>.
	
	Only X errors are considered by this function.  To count
	Z errors see countVerifiedZOnly.
	
	@param prepPairA: The encoding circuit pair for X-error verification on block A.
	@type prepPairA:  tuple  - (Locations, Locations)
	@param prepPairB: The encoding circuit pair for X-error verification on block B.
	@type prepPairB:  tuple  - (Locations, Locations)
	@param settings:  Z-error verification settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''
		
	# |+> is simply the dual of |0>
	resultZero = countVerifiedZeroZOnly(prepPairA, prepPairB, settings, noise)
	
	return CountResult(resultZero.counts, 
					   None, 
					   resultZero.locTotals.dual(), 
					   resultZero.prAccept, 
					   resultZero.prBad)

def countVerifiedPlusZOnly(prepPairA, prepPairB, settings, noise):
	'''
	Computes the weighted Z-error counts for a verified encoded |+>.
	
	Only Z errors are considered by this function.  To count
	X errors see countVerifiedZOnly.
	
	@param prepPairA: The encoding circuit pair for X-error verification on block A.
	@type prepPairA:  tuple  - (Locations, Locations)
	@param prepPairB: The encoding circuit pair for X-error verification on block B.
	@type prepPairB:  tuple  - (Locations, Locations)
	@param settings:  Z-error verification settings.
	@type settings:   ComponentSettings
	@param noise:     The (X-only, Z-only, XZ) noise models
	@type noise:      dict
	@rtype:           CountResult             
	'''
	
	# |+> is simply the dual of |0>
	resultZero = countVerifiedZeroXOnly(prepPairA, prepPairB, settings, noise)
	
	return CountResult(resultZero.counts, 
					   None, 
					   resultZero.locTotals.dual(), 
					   resultZero.prAccept, 
					   resultZero.prBad)

	