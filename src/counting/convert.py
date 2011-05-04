'''
Created on 2010-12-22

This file contains functions for converting weighted error counts of one type,
to error counts of another type.

@author: adam
'''
import counting.probability
import math
import logging

logger = logging.getLogger('counting.convert')

def zCountsFromXZCounts(xzCounts):
	'''
	Re-indexes 23-qubit (one block) encoded |0> XZ error counts into counts indexed only by Z-error syndrome.
	
	@param xzCounts: XZ error counts. Counts must be indexed by [k][X-syndrome (12-bits) | Z-syndrome (11-bits)].
	@type xzCounts: list
	@return: A list of counts indexed by [k][Z-syndrome (11-bits)]  
	'''
	mask11 = (1<<11) - 1
	zCounts = [0] * len(xzCounts)
	
	for k in range(len(zCounts)):
		zCountsK = [0] * (1<<11)
		
		# Sum over the X-syndromes
		for s, count in xzCounts[k].iteritems():
			zCountsK[s & mask11] += count
			
		zCounts[k] = zCountsK
		
	return zCounts

def xCountsFromXZCounts(xzCounts):
	'''
	Re-indexes 23-qubit (one block) encoded |0> XZ error counts into counts indexed only by X-error syndrome.
	
	@param xzCounts: XZ error counts. Counts must be indexed by [k][X-syndrome (12-bits) | Z-syndrome (11-bits)].
	@type xzCounts: list
	@return: A list of counts indexed by [k][X-syndrome (12-bits)]  
	'''

	xCounts = [0] * len(xzCounts)
	
	for k in range(len(xCounts)):
		xCountsK = [0] * (1<<12)
		
		# Sum over the Z-syndromes
		for s, count in xzCounts[k].iteritems():
			xCountsK[s >> 11] += count
			
		xCounts[k] = xCountsK
		
	return xCounts


def rescaleXZCounts(xzCounts, xzTotals, xTotals, type, noise):
	'''
	Converts XZ error counts into Z-only error counts.
	
	Error counts are coefficients of likelyhood terms of the form (gamma / (1-W*gamma)).
	The value of W is different for XZ counts than it is for X-only or Z-only counts.
	In order to combine (e.g. add, subtract, convolve) XZ counts with X-/Z-only counts
	the XZ counts must be rescaled according to the likelyhood term of the other
	count type.
	
	This function rescales the XZ counts so that they can be combined with the other
	type of counts, while maintaining an upper bound on the overall probability
	calculation. 
		
	@param xzCounts: The XZ error counts indexed as [k][X/Z syndrome]
	@type xzCounts: list
	@param xzTotals: The location count totals for the XZ error counts
	@type xzTotals: LocationCount
	@param xTotals: The location count totals for the X-/Z-only error counts.
	@type xTotals: LocationCount
	@param type: Either 'X' for X-only or 'Z' for Z-only
	@type type: string
	@param noise: The noise models
	@type noise: dict
	@return: The scaled error counts indexed as [k][X/Z syndrome], and an estimate of the error
			introduced by the scaling indexed by [k]
	@rtype: tuple  
	'''

	def scaleCount(c, scalingFactor):
		scaled = scalingFactor * c
		return long(math.floor(scaled))
	
	def scalingError(c, scalingFactor):
		scaled = scalingFactor * c
		intScaled = long(math.floor(scaled))
		return scaled - intScaled
	
	noiseX = noise[type]
	noiseXZ = noise['XZ']
	gMin, gMax = noiseX.noiseRange()
	
	AXZ = counting.probability.likelyhoodPrefactorPoly(xzTotals, noiseXZ)
	AX = counting.probability.likelyhoodPrefactorPoly(xTotals, noiseX)
	print 'gmin=', gMin, 'gmax=', gMax
	print 'AXZ=', AXZ(gMax), AXZ(gMin)
	print 'AX=', AX(gMax), AX(gMin)
	scalingCoeff = max(AXZ(gMax) / AX(gMax), AXZ(gMin) / AX(gMin))
	logger.debug('scalingCoeff=%f', scalingCoeff)
	
	scalingConstant = noiseX.likelyhood() / noiseXZ.likelyhood()
	scalingConstant = max(scalingConstant(gMax), scalingConstant(gMin))
	print 'scalingConstant=', scalingConstant
	logger.debug('scalingConstant=%f', scalingConstant)
	
	scaledCounts = [0] * len(xzCounts)
	scalingErrors = [0] * len(xzCounts)
	for k in range(len(xzCounts)):
		scalingFactor = (scalingConstant ** k) * scalingCoeff
		scaledCounts[k] = [scaleCount(c, scalingFactor) for c in xzCounts[k]]
		scalingErrors[k] = float(sum(scalingError(c, scalingFactor) for c in xzCounts[k]))
		
	return scaledCounts, scalingErrors


def normalizedCounts(counts, prAccept, nC, nPrep, nMeas, nR, gMin, gMax):
	'''
	Returns error counts normalized so that the probability of error based on any one count is <= 1.
	'''
	
	nPM = nPrep + nMeas
	
	# Lower bound on likelyhood prefactor.
	A = pow(1-12*gMax, nC) * pow(1-8*gMax, nR) * pow(1-4*gMax, nPM)
	
	preFactor = prAccept / A
	
	# Upper bound
	#gFactor = (1 - 12*gMax) / (4 * gMin)
	gFactor = (1 - 12*gMax) / (gMin)

	normCounts = [0] * len(counts)
	for k in range(len(counts)):
		countMax = preFactor * (gFactor ** k)
		
		# We want to maintain counts as integers
		countMax = long(math.ceil(countMax))
		
		normCounts[k] = [min(c, countMax) for c in counts[k]]

	return normCounts

