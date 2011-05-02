# golayCode.py :
# This is Python code for setting up the 23-qubit Golay code.  
# One interim application is to simulate CNOT exRecs with depolarizing noise 
# to estimate how high the threshold might be.  It should be possible to get 
# good upper bounds on the probability of each syndrome to make the first level 
# of simulation entirely rigorous.  Then we could plug in to a more conservative 
# analysis.  How good this can be will depend on how the CNOT exRec simulation 
# ends up working out.  
#
# Portions of the Golay code code are adapted from code by Robert Morelos-Zaragoza 
# (robert@spectra.eng.hawaii.edu).
#
# Ben Reichardt, 12/28/2009
#

from util.counterUtils import weight, etostr, readBinaryString, parity
import operator
import random
import util.counterUtils

def permuteBinaryString(string, permutation):
	"""Permutes the bits of the string according the input permutation.
	
	Note that bits beyond the length of the permutation will be cleared.
	>>> etostr(permuteBinaryString(readBinaryString('1...11'), [0, 3, 2, 1, 5, 4]), 6)
	'011001'
	"""
	return reduce(operator.ixor, [((string>>i)&1)<<sigmai for i,sigmai in enumerate(permutation)])
def permuteListByInverse(list, permutation):
	"""Permutes the input list by the *inverse* of the input permutation.
	
	>>> permuteListByInverse(['a','b','c','d'], [0, 2, 3, 1])
	['a', 'c', 'd', 'b']
	"""
	return [list[i] for i in permutation]
def permuteList(list, permutation):
	"""Permutes the input list, like Mathematica's Permute by with indices starting at 0.  
	
	Allocates new memory.  
	>>> permuteList(['a','b','c','d'], [0, 2, 3, 1])
	['a', 'd', 'b', 'c']
	"""
	result = [0] * len(permutation)
	for i,sigmai in enumerate(permutation):
		result[sigmai] = list[i]
	return result
def permutationPower(permutation, power):
	"""Raises the input permutation to the input *nonnegative* power.  Slow."""
	result = range(len(permutation))
	p = power
	while p > 0: 
		result = permuteList(result, permutation)
		p -= 1
	return result

def permuteSyndrome(syndrome, permutation):
	"""Applies a qubit permutation to an error equivalence class represented by its 11-bit syndrome.
	
	An 11-bit syndrome corresponds to some 23-qubit error-equivalence class.  This permutes the 23-qubit 
	error and then reduces back to an 11-bit syndrome.  The main property it uses is that an error equivalent 
	to a given syndrome can be found by simply interpreting that syndrome as a 23-bit string (because the 
	matrix of stabilizers we are using begins with the identity)."""
	return getSyndrome(permuteBinaryString(syndrome, permutation))
def permuteSyndromesByInverse(syndromeIndexedList, permutation):
	"""Permutes the input 2^{11}-long list, indexed by syndrome, according to the inverse of the input permutation.
	
	(Using the inverse slightly simplifies the code.)
	"""
	return [syndromeIndexedList[permuteSyndrome(s, permutation)] for s in range(2**11)]

class BitPermutations:
	"""Class for generating pseudo-random permutations that preserve the code's structure.
	
	See '2010-07-22 Golay code permutation symmetry group M23.nb' for details.
	"""
	def __init__(self):
		a1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0]
		a2 = [0, 1, 16, 12, 3, 5, 8, 17, 2, 6, 11, 22, 13, 18, 19, 14, 9, 10, 4, 21, 15, 20, 7]
		self.a1powers = [permutationPower(a1, p) for p in range(23)]
		self.a2powers = [permutationPower(a2, p) for p in range(5)]
		self.numSteps = 30
	def pseudorandomBitPermutation(self):
		permutation = range(23)
		for _ in range(self.numSteps): 
			pow1 = random.randrange(23)
			pow2 = random.randrange(5)
			permutation = permuteListByInverse(permutation, self.a1powers[pow1])
			permutation = permuteListByInverse(permutation, self.a2powers[pow2])
		return permutation

#  There are two different problems that we need to be able to solve:
#  1. find the minimum correction needed for a given error pattern (given set of failed parity checks)
#  2. find the minimum weight error equivalent to the given one
# The first problem is important for implementation, the second for "debugging."
# To be equivalent, an error must have the same pattern of failed parity checks, and must also have the same
# overall parity. The minimum correction fails the same parity checks as the given error, so is a candidate
# to be equivalent.  However, it might differ in overall parity.  Therefore, solving the first problem does
# not give a solution to the second.
# A solution to the second problem has no real meaning to the first problem since it uses information not
# available to the first algorithm, namely implicitly the correct overall parity.  The first algorithm only
# knows parity checks, not the true error.

# This solves the second problem.  If logicallystabilized is set to true, then it solves the first problem.
# It does so by brute force enumeration over all the stabilizer elements.
# This incorrectly reduces ...1...1.....1..11.1..1 -- an element of the codespace with odd parity (i.e., an
# element of the superposition of the logical |1>).
# I think this is because it should be enumerating over the generators, not the parity checks (stabilizers).
# For this problem, unlike for the Hamming code, the two are not equivalent.

#def reduceError(e, logical=False): 
#	"""Reduces the input X or Z error modulo the stabilizers and possibly the logical operator.  
#	
#	Returns a minimum-Hamming-weight representative.
#	"""
#	stabilizers = [2422785, 4845570, 3610628, 7221256, 7903248, 5621792, 1009728, 2019456, 4038912, 8077824, 5405696]
#	# 01001001111100000000001 = 2422785
#	# 10010011111000000000010 = 4845570
#	# 01101110001100000000100 = 3610628
#	# 11011100011000000001000 = 7221256
#	# 11110001001100000010000 = 7903248
#	# 10101011100100000100000 = 5621792
#	# 00011110110100001000000 = 1009728
#	# 00111101101000010000000 = 2019456
#	# 01111011010000100000000 = 4038912
#	# 11110110100001000000000 = 8077824
#	# 10100100111110000000000 = 5405696
#	if (logical):
#		stabilizers.append((1<<23)-1)
#	return util.counterUtils.reduceError(e, stabilizers)


stabGenerators = map(readBinaryString, [
		"1.1..1..11111..........",
		"1111.11.1....1.........",
		".1111.11.1....1........",
		"..1111.11.1....1.......",
		"...1111.11.1....1......",
		"1.1.1.111..1.....1.....",
		"1111...1..11......1....",
		"11.111...11........1...",
		".11.111...11........1..",
		"1..1..11111..........1.",
		".1..1..11111..........1"
	])

#def reduceError(e, logical=False):
#	c1 = corrections[getSyndrome(e)]
#	if not logical: 
#		return c1
#	c2 = corrections[getSyndrome(e ^ ((1<<23)-1))]
#	w1 = weight(c1)
#	w2 = weight(c2)
#	return c1 if w1 < w2 else c2

def getSyndrome(pattern): 
	"""Compute the syndrome corresponding to the given pattern.
	
	Returns which parity checks are violated, an 11 bit number from 0 to (1<<10)-1.
	That is, the remainder after dividing the pattern (when considering it as 
	the vector representation of a polynomial) by the generator polynomial, GENPOL.
	In the program this pattern has several meanings: (1) pattern = information
	bits, when constructing the encoding table; (2) pattern = error pattern,
	when constructing the decoding table; and (3) pattern = received vector, to
	obtain its syndrome in decoding.
	"""
	X22    = 1<<22				# vector representation of X^22
	X11    = 1<<11				# vector representation of X^11
	MASK12 = (1<<23)-(1<<11)	# auxiliary vector for testing
	GENPOL = 0xc75				# generator polynomial, g(x) = x^11+x^10+x^6+x^5+x^4+x^2+1
	
	aux = X22
	if pattern >= X11:
		while pattern & MASK12: 
			while not (aux & pattern): 
				aux = aux >> 1
			pattern ^= (aux/X11) * GENPOL
	return pattern

# ---------------------------------------------------------------------
#                  Generate DECODING TABLE
# 
# An entry to the decoding table is a syndrome and the resulting value
# is the most likely error pattern. First an error pattern is generated.
# Then its syndrome is calculated and used as a pointer to the table
# where the error pattern value is stored.
# Effectively, this takes an error pattern modulo the stabilizer, and 
# modulo the logical operators (decoding_table[get_syndrome(err)] should 
# be the same as either reduceerrorGolay(err) or reduceerrorGolaylookup(err).
# ---------------------------------------------------------------------
def generateCorrectionsTable(): 
	"""
	It would take a while to precompute the corrections for all 2^23 different X errors.  
	Instead, we precompute the corrections for all 2^11 syndromes only.  This can be done 
	by simply computing the syndromes for all 0-, 1-, 2- or 3-bit errors, and using that 
	the code is perfect.
	"""
	corrections = [-1] * (1<<11)
	
	corrections[0] = 0
	for weight in range(1, 4): 
		for t in util.counterUtils.SubsetIterator(range(23), weight): 
			tint = reduce(lambda x, y: x+y, map(lambda x: 1<<x, t))
			corrections[getSyndrome(tint)] = tint
	return corrections

class Corrector:
	
	syndromeMask = (1<<11) - 1
	errorMask = (1<<23) - 1
	
	"""Class that initializes correction and decoding lookup tables, exporting corresponding functions."""
	def __init__(self):
		self.corrections = generateCorrectionsTable()
	def correctError(self, e):
		return self.corrections[getSyndrome(e)]
	def correctXError(self, e):
		return self.correctError(e)
	def correctZError(self, e):
		return self.correctError(e)
	#def decodeError(self, e):
	#	return bool(weight(e ^ self.correctError(e), 23) % 2)
	def decodeError(self, e, type):
		#return self.decodeError(e)
		return bool(weight(e ^ self.correctError(e), 23) % 2)
	def decodeSyndrome(self, s):
		e = self.getError(s)
		return bool(weight(e ^ self.corrections[s & self.syndromeMask], 23) % 2)
	def reduceError(self, e, logical=False): 
		"""This is the same as reduceErrorGolay, except preprocessed according to the syndrome."""
		c1 = self.correctError(e)
		if not logical: 
			return c1
		c2 = self.correctError(e ^ ((1<<23)-1))
		w1 = weight(c1)
		w2 = weight(c2)
		return c1 if w1 < w2 else c2
	def getSyndrome(self, e):
		"""Returns the 11-bit error syndrome."""
		return getSyndrome(e)
	def getLogicalSyndrome(self, e):
		"""Returns a 12-bit error syndrome, the parity followed by the 11-bit error syndrome."""
		s = self.getSyndrome(e)
		w = weight(e, 23)
		s ^= (w&1)<<11
		return s
	
	def getError(self, s):
		'''Returns an error that generates the given (possibly logical) syndrome'''
		sMask = self.syndromeMask
		parityBit = bool(s & (sMask + 1))
		e = s & sMask
		if parity(e, 11) != parityBit:
			e ^= self.errorMask
		
		return e
	
	def correctLogicalSyndrome(self, s):
		'''Returns a syndrome corresponding to the appropriate correction for the given logical syndrome.'''
		# Strip off the parity bit
		s &= self.syndromeMask
		correction = self.corrections[s]
		return self.getLogicalSyndrome(correction)
	
	def correctSyndrome(self, s):
		s &= self.syndromeMask
		correction = self.corrections[s]
		return self.getSyndrome(correction)
		
	
cnots = [
	[12,14,18,15,11,21,13], # target of CNOT gates in time steps 1-7 from source qubit 0
	[19,13,15,16,12,22,14], # .. from qubit 1 (21 in reversed-index matrices, sorry)
	[20,16,12,18,17,11,21], # 2	
	[21,12,19,13,22,18,17], # 3
	[15,19,11,22,21,20,12], # 4
	[16,18,22,11,15,14,20], # 5
	[18,11,16,17,14,13,19], # 6
	[17,15,14,12,20,19,18], # 7
	[13,21,20,19,18,16,15], # 8
	[14,22,21,20,19,17,16], # 9
	[22,20,17,14,13,12,11] # 10
]

newCnots = [0] * len(cnots[0])
for rnd in range(len(cnots[0])):
	thisRound = [0] * len(cnots)
	for sbit in range(len(cnots)):
		thisRound[sbit] = (sbit, cnots[sbit][rnd])
		
	newCnots[rnd] = thisRound

cnots = newCnots



def errorWeights():
	"""Computes the minimum-weight errors in every error equivalence class.
	
	Returns a 2^11 x 2 array containing for each syndrome the minimum weight of any error 
	that generates that syndrome, with weight either even (first entry) or odd (second entry).
	This is quite slow, so I have now hard-coded the answer into the first line.
	
	Note that there are (23 choose w) different weight-w Z errors, for w in {0,1,2,3}.  
	There are (23 choose w) different weight-w X errors, for w in {0,1,2,3}, and (23 choose 7-w)
	weight-w X errors, for w in {4,5,6,7}.
	
	>>> errorWeights()[0]		# (1<<23)-1 has trivial syndrome and odd weight
	[0, 7]
	"""
	return [[0, 7], [6, 1], [6, 1], [2, 5], [6, 1], [2, 5], [2, 5], [4, 3], [6, 1], [2, 5], [2, 5], [4, 3], [2, 5], 
		[4, 3], [4, 3], [4, 3], [6, 1], [2, 5], [2, 5], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [6, 1], [2, 5], [2, 5], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [2, 5], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [6, 1],
		[2, 5], [2, 5], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [6, 1], [2, 5], 
		[2, 5], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [2, 5], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [6, 1], [2, 5], [2, 5], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], 
		[4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3],
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3],
		[2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [2, 5], [6, 1], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3],
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [6, 1], 
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3],
		[4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [6, 1], [2, 5], [2, 5], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [4, 3],
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3],
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3],
		[4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], 
		[4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [2, 5], [6, 1], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [2, 5], [4, 3], [6, 1], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], 
		[4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [2, 5], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [6, 1], [2, 5], 
		[4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [6, 1], [2, 5], [2, 5], 
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3],
		[4, 3], [4, 3], [2, 5], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3],
		[4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3],
		[4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [6, 1], [4, 3], [2, 5], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [2, 5], [2, 5], 
		[6, 1], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5],
		[4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [2, 5], [2, 5], [6, 1], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [2, 5], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3],
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], 
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3],
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [6, 1], [2, 5], [4, 3], [4, 3], [2, 5],
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [2, 5], [2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [6, 1], [2, 5], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [2, 5], [6, 1], [4, 3], [4, 3], [4, 3], [2, 5], 
		[4, 3], [4, 3], [4, 3], [2, 5], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [2, 5], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [2, 5], [4, 3], [4, 3], [4, 3], [6, 1], [2, 5], [2, 5], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], 
		[4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [2, 5], [4, 3], [4, 3], [4, 3], [4, 3], 
		[2, 5], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3], [4, 3]]
	weights = [[23,23] for s in range(2**11)]
	for error in range(2**23):
		if not error % (1<<17):
			print etostr(error, 23)
		w = weight(error)
		s = getSyndrome(error)
		weights[s][w%2] = min(weights[s][w%2], w)
	return weights

# ---------------------------------------------------------------------
#  The following code is all for verifying that we indeed have a [[23,1,7]] code, with no bugs.  
# ---------------------------------------------------------------------

def generateGolayEncodingTable(): 
	encodingTable = [-1] * 4096
	for pattern in range(4096): 
		temp = pattern << 11									# multiply information by X^{11}
		encodingTable[pattern] = temp + getSyndrome(temp)		# add redundancy
	print "Golay [[23,1,7]] generator matrix is:"
	for j in range(12): 
		print etostr(encodingTable[1 << j], 23)
	return encodingTable

def printParities(err): 
	print etostr(getSyndrome(err))

#def checkminimumdistance(): 
#	minweight = 23;
#	for pattern in range(1, 4096): 
#		if weight(encodingTable[pattern], 23) < minweight: 
#			minweight = weight(encodingTable[pattern], 23)
#	print "minimum weight of codeword is %d" % minweight
#
#def checkweightdistribution(): 
#	weights = [0] * 24
#	for pattern in range(4096): 
#		weights[weight(encodingTable[pattern], 23)] += 1
#	print "Weight distribution is:"
#	for i in range(24):
#		if weights[i]: 
#			print "%d: %d" % (i, weights[i])
#
#def checklinearity(): 
#	result = True
#	for pattern in range(4096): 
#		code = 0;
#		for j in range(12): 
#			if pattern & (1<<j):
#				code ^= encodingTable[1<<j]
#		if encodingTable[pattern]-code:
#			result = False
#	print "Code is %slinear" % "" if result else "not "
#
#def checkGolayCodeSetup(): 
#	print "\nChecking Golay code setup..."
#	checkminimumdistance()
#	checkweightdistribution()
#	checklinearity()
#	
#	print "Logical X/Z is:"
#	print etostr(encodingTable[4096-1])
#
#	test = readBinaryString("...1...1.....1..11.1..1")
#	test = readBinaryString(".....1.1.......1.......")
#	test = readBinaryString("......1......111x......")
#	print "syndrome= %d" % getSyndrome(test)
#	a1, a2, a3 = reduceError(test, False), reduceError(test, False), corrections[getSyndrome(test)]
#	print "\ntest= %d, lookup= %d, brute= %d, table= %d" % (test, a1, a2, a3)
#	for t in [test, a1, a2, a3]: 
#		print etostr(t)
#		printParities(t)
#	for k in range(2048): 
#		if corrections[k] == test: 
#			print "error is lowest weight possible for the syndrome"
#			print "k= %d" % k
#			print etostr(k)
#	print "violated parity checks are "
#	print etostr(getSyndrome(test))
#
#	print "lowest weight rep'n of the error, logical op + error are "
#	print etostr(reduceError(test, False))
#	print etostr(reduceError(test ^ ((1<<23)-1), False))
#	
#	print "\n"
#	q = [0] * 11
#	q[10] = ((1<<12)+(1<<10)+(1<<7)+(1<<4)+(1<<3)+(1<<2)+(1<<1)+1) << 10
#	for k in range(9, -1, -1): 
#		q[k] = q[k+1] >> 1
#		if q[k] & (1<<10):
#			q[k] ^= q[10]
#	for k in range(11):
#		print etostr(q[k])
#		printParities(q[k])
#		print q[k]
#	print ""
#	for k in range(12): 
#		printParities(encodingTable[1<<k])
#	print "...Finished checking Golay code setup\n"
#	
#	print "-- Just one last check; let me be sure I wrote down the 11 X and 12 Z stabilizers of encoded |0> correctly..."
#	checks = map(readBinaryString, [
#		"1.1..1..11111..........",
#		"1111.11.1....1.........",
#		".1111.11.1....1........",
#		"..1111.11.1....1.......",
#		"...1111.11.1....1......",
#		"1.1.1.111..1.....1.....",
#		"1111...1..11......1....",
#		"11.111...11........1...",
#		".11.111...11........1..",
#		"1..1..11111..........1.",
#		".1..1..11111..........1"
#	])
#	for check in checks: 
#		print reduceError(check, False)
#	print "\n"
#
#	checks = map(readBinaryString, [
#		"11...111.1.1...........",
#		".11...111.1.1..........",
#		"1111.11.1....1.........",
#		".1111.11.1....1........",
#		"..1111.11.1....1.......",
#		"11.11..11.......1......",
#		".11.11..11.......1.....",
#		"..11.11..11.......1....",
#		"11.111...11........1...",
#		"1.1.1..1.11.........1..",
#		"1..1..11111..........1.",
#		"1...111.1.1...........1"
#	])
#	for check in checks: 
#		print reduceError(check, True)
#	print "\n\n"
#	
#	print "\n\n\nrunning more tests...\n"
#	correction 	= readBinaryString("......1.......1.x...x..")
#	test 		= readBinaryString("......1.......1.1...1..")
#	test ^= correction
#	print etostr(reduceError(test, False))
#
#def makeparitycheckmatrix(): 
#	p = [0] * 11
#	h = (1<<12)+(1<<10)+(1<<7)+(1<<4)+(1<<3)+(1<<2)+(1<<1)+1
#	p[0] = ((1<<12)+(1<<10)+(1<<7)+(1<<4)+(1<<3)+(1<<2)+(1<<1)+1) << 10
#	for i in range(1, 11): 
#		p[i] = p[i-1] >> 1
#		if p[i] & (1<<10):
#			p[i] ^= p[0]
#	print "Golay parity check matrix is:"
#	for i in range(11): 
#		print etostr(p[i])
#		#print etostr(p[i] ^ (1 << (10-i)))
#		#print "%d %d" % (p[i], p[i] ^ (1 << (10-i)))
#	print "Inner product of generator with check matrix:"
#	pattern = 1
#	while pattern < 4096: 
#		for j in range(11): 
#			print weight(p[j] & encodingTable[pattern], 23)
#		pattern *= 2
#	print ""
#
#	logical = (1<<23)-1
#	for i in range(11): 
#		logical ^= p[i]
#	print "Logical X/Z is:"
#	print etostr(logical)
#	print logical
#
#	for i in range(11): 
#		print weight(p[i] & logical, 23)
#
#	print "Check matrix including logical X/Z is:"
#	print etostr(logical)
#	for i in range(11): 
#		c = p[i]
#		if (c >> 11) & 1:
#			c ^= logical
#		print etostr(c)
#	print ""

#encodingTable = generateGolayEncodingTable()
#corrections = generateCorrectionsTable()
#makeparitycheckmatrix()
#checkGolayCodeSetup()


# to run the doctests, run python or python -v directly on this script
if __name__ == "__main__":
	import doctest
	doctest.testmod()