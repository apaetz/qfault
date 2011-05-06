'''
Created on 2011-05-06

@author: adam
'''
import math
from util import listutils

def getStats(sample):
	N = float(len(sample))
	mean = sum(sample) / N
	sampleSigma = math.sqrt(sum((s - mean)**2 for s in sample) / (N-1))
	sampleErr = sampleSigma / math.sqrt(N)
	
	# Use 1.96 for 95% confidence interval
	return mean, sampleSigma, 1.96 * sampleErr 

def mean(sample):
	return sum(sample) / float(len(sample))

def std(sample):
	m = mean(sample)
	N = float(len(sample))
	return math.sqrt(sum((s - m)**2 for s in sample) / (N-1))

def stdErr(std, N):
	return std / math.sqrt(N)

def err95(se):
	return 1.96 * se
	
def stdProduct(stds, vals):
	
	prod = listutils.mul(vals)
	varTerms = [(prod*s/vals[i])**2 for i,s in enumerate(stds)]
	return math.sqrt(sum(varTerms))
		
	