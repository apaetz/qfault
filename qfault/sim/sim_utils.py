#!/usr/bin/python
# by Ben Reichardt, created 12/28/2009
# 
# Utilities for simulating exRecs.
# 


import random

def rest(errorRates, errors, block, bit):
	# introduce new errors
	rate = errorRates.rest
	r = random.random()		# random float from (0,1)
	if r < rate: 
		# might be nice to be able to simulate X and Z errors separately
		r = 1 + int(3 * r / rate)	# r is now an integer in {1,2,3}, can use its two bits to assign errors
		if r&1: 
			errors['X'][block] ^= 1<<bit
		if r>>1&1:
			errors['Z'][block] ^= 1<<bit
		#print "Error! rest", util.counterUtils.errorsToStr(errors, 23, False), block, bit
	return errors	

def measX(errorRates, errors, block, bit):
	if random.random() < errorRates.measX:
		errors['Z'][block] ^= 1<<bit
		#print "Error! measX", util.counterUtils.errorsToStr(errors, 23, False), block, bit
	return errors

def measZ(errorRates, errors, block, bit):
	if random.random() < errorRates.measZ:
		errors['X'][block] ^= 1<<bit
		#print "Error! measZ", util.counterUtils.errorsToStr(errors, 23, False), block, bit
	return errors

def prepX(errorRates, errors, block, bit):
	if random.random() < errorRates.prepX:
		errors['Z'][block] ^= 1<<bit
		#print "Error! prepX", util.counterUtils.errorsToStr(errors, 23, False), block, bit
	return errors

def prepZ(errorRates, errors, block, bit):
	if random.random() < errorRates.prepZ:
		errors['X'][block] ^= 1<<bit
		#print "Error! prepZ", util.counterUtils.errorsToStr(errors, 23, False), block, bit
	return errors
	
def cnot(errorRates, errors, block1, bit1, block2, bit2):
	# propagate existing errors
	b1, b2 = errors['X'][block1], errors['X'][block2]
	b2 = b2 ^ ((b1>>bit1&1)<<bit2)
	errors['X'][block2] = b2
	b1, b2 = errors['Z'][block1], errors['Z'][block2]
	b1 = b1 ^ ((b2>>bit2&1)<<bit1)
	errors['Z'][block1] = b1
	
	# introduce new errors
	rate = errorRates.cnot
	r = random.random()		# random float from (0,1)
	#if block1 == 'a1' and block2 == 'a1' and bit1 == 1 and bit2 == 19:
	#	r = rate * 9/float(15)
	#print "cnot", util.counterUtils.errorsToStr(errors, 23, False), block1, bit1, block2, bit2
	if r < rate: 
		# might be nice to be able to simulate X and Z errors separately
		r = 1 + int(15 * r / rate)	# r is now an integer in {1,2,...,15}, can use its four bits to assign errors
		if r&1: 
			errors['X'][block1] ^= 1<<bit1
		if r>>1&1:
			errors['Z'][block1] ^= 1<<bit1
		if r>>2&1:
			errors['X'][block2] ^= 1<<bit2
		if r>>3&1:
			errors['Z'][block2] ^= 1<<bit2
		#print "Error! cnot", util.counterUtils.errorsToStr(errors, 23, False), block1, bit1, block2, bit2
	return errors

#def checkMalignantUsingDecoder(combinedErrors):
#	"""Checks whether the decoders can be pulled past, as in the AGP paper.
#	
#	It applies the decoders to the first and last checkpoints, and sees if they agree. 
#	"""
#	decode = lambda block: [decoders[etype][combinedErrors[etype][block]] for etype in ['X', 'Z']]
#	Da1, Db1, Da3, Db3 = map(decode, ['Da1', 'Db1', 'Da3', 'Db3'])
#	# having decoded the blocks at each checkpoint, we see whether the decoders commute
#	# for example, if Da1 has a logical X error, and Db1 doesn't, then Da3 and Db3 both should
#	isXMalignant = Da3[0] != Da1[0] or Db3[0] != (Da1[0] ^ Db1[0])
#	isZMalignant = Da3[1] != (Da1[1] ^ Db1[1]) or Db3[1] != Db1[1]
#	isMalignantUsingDecoder = isXMalignant or isZMalignant
#		# it is also interesting to see what location sets are malignant just for X errors, or 
#		# are simultaneously malignant for X and Z errors; see the following commented-out lines
#	#isMalignantUsingDecoder = isXMalignant
#	#isMalignantUsingDecoder = isXMalignant and isZMalignant
#	#if isMalignantUsingDecoder: 
#	#	print Da1, Db1, Da3, Db3
#	#	print util.counterUtils.errorsToStr(combinedErrors, 9)
#	return isMalignantUsingDecoder

class ErrorRates():
	def __str__(self):
		return "cnot= %.6f, rest= %.6f, prepX= %.6f, prepZ= %.6f, measX= %.6f, measZ = %.6f" \
			% (self.cnot, self.rest, self.prepX, self.prepZ, self.measX, self.measZ)
	
# to run the doctests, run python or python -v directly on this script
if __name__ == "__main__":
	import doctest
	doctest.testmod()