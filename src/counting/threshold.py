'''
Created on Apr 29, 2011

@author: Adam
'''
import logging
logger = logging.getLogger('threshold')

def pseudoThresh(prFail, pMin, pMax, tolerance=1e-6):
	prFail0 = lambda p: p
	prFail1 = lambda p: prFail(prFail0(p))
	return findIntersection(prFail0, prFail1, pMin, pMax, tolerance)


def asymptoticThresh(level1Events, level2Events, gMax):
		
	thresholds = []
	for locName in level1Events.keys():
		logger.info('Computing threshold for %s', locName)
		pr1 = level1Events[locName]
		pr2 = level2Events[locName]
		
#		plots = [prFail1, prFail2]
#		labels = [locName + '1', locName + '2']
#		plotPolyList(plots, gMin, gMax, locName + 'level1and2', labelList=labels, xLabel=r'$\gamma$', numPoints=10)
#		plotPolyList(plots, GammaZ(gMin), GammaZ(gMax), locName + 'level1and2', labelList=labels, xLabel=r'$\gamma$', numPoints=10)
		
		thresh = findIntersection(pr1, pr2, 0., gMax)
		logger.info('%s threshold = %s', locName, thresh)
		thresholds.append(thresh)

	threshold = min(thresholds)
	
	return threshold



def findIntersection(prFail0, prFail1, pMin, pMax, tolerance=1e-8):
	'''
	Computes a lower bound on the threshold over the inverval [pMin, pMax] to within tolerance.
	Computation is done via recursive binary search.
	'''
	
	p = (pMax + pMin) / 2
	p0 = prFail0(p) 
	p1 = prFail1(p)
	logger.info('pMin=%f, pMax=%f', pMin, pMax)
	logger.info('p={0}, p0={1}, p1={2}'.format(p, p0, p1))
	diff = p0 - p1
	if diff < 0:
		# threshold is less than p
		
		if (p - pMin < tolerance):
			# We're really close to pMin and haven't reached the pseudothreshold.
			if prFail1(pMin) > prFail0(pMin):
				# threshold is less than pMin. Give up.
				return None
			return pMin

		return findIntersection(prFail0, prFail1, pMin, p, tolerance)
	
	if (diff <= tolerance) or (pMax - p < tolerance):
		# threshold may be greater than p, but is very close.
		return p
	
	
	# Pseudothreshold is greater than p
	return findIntersection(prFail0, prFail1, p, pMax, tolerance)