'''
This file contains functions for bounding polynomials that represent malignant event probabilities.

Created on 2010-10-24

@author: adam
'''
from counting.countParallel import iterParallel
from util import listutils
from util.plotting import evalExprList
import logging


logger = logging.getLogger('bounding')

def samplePolys(polyList, xMin, xMax, numPoints):
	'''
	Returns a sample of y-values of the given polynomials over 
	the range [xMin, xMax].
	This function is designed to sample functions that behave
	exponentially over the interval.  It samples many points
	at the beginning and end of the interval.
	'''
#	delta = float(xMax-xMin)/numPoints
#	X = [xMin + i*delta for i in range(numPoints)]
	# Typical curves are exponential, so we want to plot more points
	# at the larger end of the range.
	#X = [xMin] * numPoints
	#for i in range(1,numPoints):
	#	X[i] = X[i-1] + (xMax - X[i-1]) * (i/float(numPoints-1))**10
	exp = [2 - 1.5*i/float(numPoints-1) for i in range(numPoints)]
	X = [xMin + (xMax-xMin) * (i/float(numPoints-1))**exp[i] for i in range(numPoints)]
	
	# Sample each polynomial in parallel
	sampleResults = iterParallel(polyList, evalExprList, [X])
	
	return [r.get() for r in sampleResults]


def computeMax(polyList, xMin, xMax, numPoints=1000):
	'''
	Returns a polynomial that upper bounds all of the given polynomials over the
	specified range.  All polynomials in the list must be monotone nondecreasing.
	'''
	
	Ylist = samplePolys(polyList, xMin, xMax, numPoints) 
			
	# Find the polynomial with the largest value at xMax.
	# Use this polynomial as a reference for which to
	# construct a scaled polynomial that upper bounds
	# all of the others.
	refIndex = 0
	for i, Y in enumerate(Ylist):
		if Y[-1] > Ylist[refIndex][-1]:
			refIndex = i
			
	refY = Ylist[refIndex]
	
	otherIndices = set(range(len(Ylist)))
	otherIndices.remove(refIndex)
	offset = max(Ylist[i][0] for i in otherIndices)
	
	scalingFactor = max([computeUBScalingFactor(Ylist[i], refY) for i in otherIndices])
	
	logger.info('refIndex=%s, scalingFactor=%s, offset=%s', refIndex, scalingFactor, offset)
	return polyList[refIndex] * scalingFactor + offset



def computeUBScalingFactor(Y, refY):
	'''
	Let P be the polynomial represented by Y, and refP the polynomial
	represented by the reference refY. This function returns a scalar 
	alpha such that alpha*refP >= P.
	
	Both Y and refY must be monotone non-decreasing.
	
	Y 		-- An ordered list of y-values for P.
	refY 	-- The corresponding y-values of refP.
	'''
				
	ratios = [Y[x+1]/refY[x] for x in range(len(Y)-1)]
	maxIndex = listutils.imax(ratios)
	scalingFactor = max(ratios)
	logger.info('scalingFactor=%s, maxIndex=%s', scalingFactor, maxIndex)
	
	return scalingFactor
	
	

def computeWeights(polys, xMin, xMax, numPoints=1000):
	'''
	Returns a reference polynomial P, and a set of weights {alpha_i}
	such that alpha_i*P >= polys[i] over the interval [xMin,xMax].
	All polynomials in polys must be monotone non-decreasing over the
	interval [0,xMax].
	'''

	polyY = samplePolys(polys, xMin, xMax, numPoints)		
	derivs = [p.diff() for p in polys]
	derivMaxes = [deriv(xMax) for deriv in derivs]
			
	# Find the polynomial with the maximum derivative (at xMax).
	# Then, create a reference by successively dividing it
	# so that it touches the smallest of the polynomials.  
	refPolyIndex = listutils.imax(derivMaxes)
	logger.info('refPoly=%s', refPolyIndex)
	
	# Next, find the maximum value of any of the polynomials
	# (other than the reference) at xMin.
	otherIndices = set(range(len(polyY)))
	otherIndices.remove(refPolyIndex)
	offsets = [polyY[i][0] for i in otherIndices]
	#offset = max(polyY[i][0] for i in otherIndices)
	#logger.info('offset=%s', offset)			
				
	maxY = polyY[refPolyIndex]
	ratios = [maxY[-1] / Y[-1] for Y in polyY]
	ratio = max(ratios) * 4. # Multiply (shrink the reference) to reduce error introduced by integer ratios.
	logger.info('ratio=%s', ratio)

	# Scale the reference.
	refY = [y/ratio for y in maxY]
	
	# Ratios at the beginning of the range can be large.  To avoid
	# overweighting, we can add a small offset to reduce ratios
	# for small values of x.  This offset also ensures maximality
	# over the interval [0,xMin].
	minWeights = [polyY[i][0]/refY[0] for i in otherIndices]  # All of the weights will be at least this large.
	
	# This offset needs to obey weight[i] * refOffset >= poly[i](xMin) for each
	# index i (other than the polynomial from which the reference was constructed).
	refOffset = max(offsets[i]/minWeights[i] for i in range(len(offsets)))
	logger.info('refOffset=%s', refOffset)
	
	refY = [y + refOffset for y in refY]
	refPoly = polys[refPolyIndex] / ratio + refOffset
	
	weights = [computeUBScalingFactor(Y,refY) for Y in polyY]
	weights[refPolyIndex] = ratio	
		
	return refPoly, weights