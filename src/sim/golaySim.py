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

from sim.simulateUtils import cnot, measX, measZ, prepX, prepZ, rest
from util.counterUtils import weight, etostr, readBinaryString, parity
import golay.ancillaPrep
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


def prepareUnverifiedAncilla(errorRates, roundOrder, cnots, blockName='a'):
	"""Simulates preparation of a logical |0> ancilla with the input round order.
	
	The round order should be a permutation of range(7).
	"""
	b = blockName
	errors = {'X':{b:0},'Z':{b:0}}	# initialize the errors to all zero
	touched = set()
	
	# If a qubit is ever the target of a CNOT, then it is prepared as |0>,
	# otherwise it is prepared as |+>.
	prepZQubits = set([tgt for rnd in cnots for _, tgt in rnd])
	prepXQubits = set(range(23)) - prepZQubits
	
	for i in prepXQubits:
		prepX(errorRates, errors, b, i)
	for i in prepZQubits:
		prepZ(errorRates, errors, b, i)

	for rnd in roundOrder:			# simulate each of the seven rounds sequentially
		restingQubits = set(touched)
		for sbit, tbit in cnots[rnd]:
			cnotBits = set([sbit, tbit])
			touched.update(cnotBits)
			restingQubits -= cnotBits
			cnot(errorRates, errors, b, sbit, b, tbit)
	
		# in each round, some qubits may not touched (there are at most 11 CNOTs touching 22 qubits)
		# insert a rest location for those qubits, provided they have been initialized
		for rbit in restingQubits:
			rest(errorRates, errors, b, rbit)
				
	return errors

encodingSequences = [
	[ 0, 1, 3, 2, 4, 5, 6 ], 
	[ 0, 1, 3, 2, 4, 5, 6 ], 
	[ 6, 5, 4, 2, 3, 1, 0 ],
	[ 6, 5, 4, 2, 3, 1, 0 ],
	[ 5, 1, 6, 3, 4, 2, 0 ], 
	[ 0, 2, 4, 3, 6, 1, 5 ]
]

prepStats = {
			 'prepA0': {'success': 0, 'attempts': 0},
			 'prepA0v': {'success': 0, 'attempts': 0},
			 'prepA2': {'success': 0, 'attempts': 0},
			 'prepA4': {'success': 0, 'attempts': 0},
			 'prepA': {'success': 0, 'attempts': 0},
			 'prepA_X1': {'success': 0, 'attempts': 0},
			 'prepA_Z1': {'success': 0, 'attempts': 0},
			 'prepA_X2': {'success': 0, 'attempts': 0},
			 'prepA_Z2': {'success': 0, 'attempts': 0},
			 'prepA_X3': {'success': 0, 'attempts': 0},
			 'prepA_Z3': {'success': 0, 'attempts': 0},
			 'prepA_X4': {'success': 0, 'attempts': 0},
			 'prepA_Z4': {'success': 0, 'attempts': 0},
			 }

def initPrepStats():
	global prepStats
	prepStats = {
			 'prepA0': {'success': 0, 'attempts': 0},
			 'prepA0v': {'success': 0, 'attempts': 0},
			 'prepA2': {'success': 0, 'attempts': 0},
			 'prepA4': {'success': 0, 'attempts': 0},
			 'prepA': {'success': 0, 'attempts': 0},
			 'prepA_X1': {'success': 0, 'attempts': 0},
			 'prepA_Z1': {'success': 0, 'attempts': 0},
			 'prepA_X2': {'success': 0, 'attempts': 0},
			 'prepA_Z2': {'success': 0, 'attempts': 0},
			 'prepA_X3': {'success': 0, 'attempts': 0},
			 'prepA_Z3': {'success': 0, 'attempts': 0},
			 'prepA_X4': {'success': 0, 'attempts': 0},
			 'prepA_Z4': {'success': 0, 'attempts': 0}, 
			 }
	

def printAcceptanceRates(prepStats):
	totalRate = 1.0
	for name, stats in prepStats.items():
		if 0 == stats['attempts']:
			continue
		rate = float(stats['success']) / stats['attempts']
		totalRate *= rate
		print name, rate
	
	print 'Pr[accept] =', totalRate
	
def getAcceptanceRate():
	global prepStats
	totalRate = 1.0
	for _, stats in prepStats.items():
		if 0 == stats['attempts']:
			continue
		rate = float(stats['success']) / stats['attempts']
		totalRate *= rate
		#print name, rate
	
	return totalRate

def getAcceptanceRates():
	global prepStats
	
	rates = {}
	totalRate = 1.0
	for name, stats in prepStats.items():
		if 0 == stats['attempts']:
			continue
		rate = float(stats['success']) / stats['attempts']
		totalRate *= rate
		rates[name] = rate
		#print name, rate
	
	rates['total'] = totalRate
	return rates

def prepareAncilla(errorRates, name='a0', eigenstate='Z'):
	"""This will simulate preparation of a verified ancilla.
	
	Note that it assumes that the error model is symmetrical under X and Z errors.  Thus to prepare an encoded 
	|+>, it prepares an encoded |0> and swaps the X and Z errors.  
	"""

	def prepA0():
		errors, attempts = prepareAndVerify(errorRates, 'X', 
			'a0', lambda: prepareUnverifiedAncilla(errorRates, encodingSequences[0], cnots, 'a0'), 
			'a1', lambda: prepareUnverifiedAncilla(errorRates, encodingSequences[1], cnots, 'a1'))
		prepStats['prepA0']['attempts'] += attempts
		prepStats['prepA0']['success'] += 1
		return errors
	def prepA2():
		errors, attempts = prepareAndVerify(errorRates, 'X', 
			'a2', lambda: prepareUnverifiedAncilla(errorRates, encodingSequences[2], cnots, 'a2'), 
			'a3', lambda: prepareUnverifiedAncilla(errorRates, encodingSequences[3], cnots, 'a3'))
		prepStats['prepA2']['attempts'] += attempts
		prepStats['prepA2']['success'] += 1
		return errors
	def prepA4():
		errors, attempts = prepareAndVerify(errorRates, 'Z', 
			'a4', lambda: prepareUnverifiedAncilla(errorRates, encodingSequences[4], cnots, 'a4'), 
			'a5', lambda: prepareUnverifiedAncilla(errorRates, encodingSequences[5], cnots, 'a5'))
		prepStats['prepA4']['attempts'] += attempts
		prepStats['prepA4']['success'] += 1
		return errors
	def prepA0v():
		errors, attempts = prepareAndVerify(errorRates, 'Z', 'a0', prepA0, 'a2', prepA2)
		prepStats['prepA0v']['attempts'] += attempts
		prepStats['prepA0v']['success'] += 1
		return errors
	
	errors, attempts = prepareAndVerify(errorRates, 'X', 'a0', prepA0v, 'a4', prepA4)
	prepStats['prepA']['attempts'] += attempts
	prepStats['prepA']['success'] += 1
	
	
	#TODO: this is a hack
#	if 0 == (prepStats['prepA']['success'] % 500):
#		printAcceptanceRates(prepStats)
	
	if eigenstate == 'Z':	# might as well clean out the other ancilla blocks
		errors = {'X':{name:errors['X']['a0']}, 'Z':{name:errors['Z']['a0']}}
	else: #eigenstate == 'X' --- assuming that the error model is symmetrical, we can just swap the X and Z errors
		errors = {'X':{name:errors['Z']['a0']}, 'Z':{name:errors['X']['a0']}}
	#print "prepared ancilla:", util.counterUtils.errorsToStr(errors, 23, True)	
	return errors



def prepareAncilla12(errorRates, name='a0', eigenstate='Z'):
	"""This will simulate preparation of a verified ancilla.
	
	Note that it assumes that the error model is symmetrical under X and Z errors.  Thus to prepare an encoded 
	|+>, it prepares an encoded |0> and swaps the X and Z errors.  
	"""
	
	roundOrder = range(7)
	
	def prepUnverified(blockName):
		return lambda: prepareUnverifiedAncilla(errorRates, roundOrder, cnots, blockName)

	def prepA_1(XorZ, name1, name2):
		errors, attempts = prepareAndVerify(errorRates, XorZ, 
										    name1, prepUnverified(name1), 
											name2, prepUnverified(name2))
		key = 'prepA_' + XorZ + '1'
		prepStats[key]['attempts'] += attempts
		prepStats[key]['success'] += 1
		return errors
	def prepA_2(XorZ):
		errors, attempts = prepareAndVerify(errorRates, XorZ, 
										    'a1', lambda: prepA_1(XorZ, 'a1', 'a2'), 
											'a3', prepUnverified('a3'))
		key = 'prepA_' + XorZ + '2'
		prepStats[key]['attempts'] += attempts
		prepStats[key]['success'] += 1
		return errors
	def prepA_Z3():
		errors, attempts = prepareAndVerify(errorRates, 'Z', 
										    'a0', lambda: prepA_1('X', 'a0', 'a1'), 
											'a2', lambda: prepA_1('X', 'a2', 'a3'))
		key = 'prepA_Z3'
		prepStats[key]['attempts'] += attempts
		prepStats[key]['success'] += 1
		return errors
	
	def prepA_X3():
		errors, attempts = prepareAndVerify(errorRates, 'X', 
										    'a0', prepA_Z3, 
											'a1', lambda: prepA_1('Z', 'a1', 'a2'))
		key = 'prepA_X3'
		prepStats[key]['attempts'] += attempts
		prepStats[key]['success'] += 1
		return errors
	def prepA_Z4():
		errors, attempts = prepareAndVerify(errorRates, 'Z', 
										    'a0', prepA_X3, 
											'a1', lambda: prepA_2('X'))
		key = 'prepA_Z4'
		prepStats[key]['attempts'] += attempts
		prepStats[key]['success'] += 1
		return errors
	def prepA_X4():
		errors, attempts = prepareAndVerify(errorRates, 'X', 
										    'a0', prepA_Z4, 
											'a1', lambda: prepA_2('Z'))
		key = 'prepA_X4'
		prepStats[key]['attempts'] += attempts
		prepStats[key]['success'] += 1
		return errors			
	
	errors = prepA_X4()
	
	if eigenstate == 'Z':	# might as well clean out the other ancilla blocks
		errors = {'X':{name:errors['X']['a0']}, 'Z':{name:errors['Z']['a0']}}
	else: #eigenstate == 'X' --- assuming that the error model is symmetrical, we can just swap the X and Z errors
		errors = {'X':{name:errors['Z']['a0']}, 'Z':{name:errors['X']['a0']}}
	#print "prepared ancilla:", util.counterUtils.errorsToStr(errors, 23, True)	
	return errors
	

def prepareAncillaNew(errorRates, name='a0', eigenstate='Z'):
	"""This will simulate preparation of a verified ancilla.
	
	Note that it assumes that the error model is symmetrical under X and Z errors.  Thus to prepare an encoded 
	|+>, it prepares an encoded |0> and swaps the X and Z errors.  
	"""

	prepCkts = golay.ancillaPrep.prepCircuits()
	cnotsA0 = permuteList(prepCkts[8], [0, 6, 2, 4, 1, 3, 5])
	cnotsA1 = permuteList(prepCkts[70], [1, 4, 5, 6, 0, 3, 2]) 
	cnotsA2 = permuteList(prepCkts[45], [6, 5, 1, 3, 2, 4, 0])
	cnotsA3 = permuteList(prepCkts[99], [0, 5, 1, 3, 6, 4, 2])
	
	return prepareAncilla4(cnotsA0, cnotsA1, cnotsA2, cnotsA3, errorRates, name, eigenstate)


goodSteaneRandomPreps = [
		[[8, 70, 45, 99], [[0, 6, 2, 4, 1, 3, 5], 
						   [1, 4, 5, 6, 0, 3, 2], 
						   [6, 5, 1, 3, 2, 4, 0], 
						   [0, 5, 1, 3, 6, 4, 2]]],
		[[29, 44, 76, 99], [[0, 5, 4, 6, 2, 3, 1], 
						    [3, 1, 0, 2, 5, 6, 4]], 
						    [6, 2, 5, 0, 1, 4, 3], 
						    [2, 1, 6, 0, 3, 4, 5]],
		[[30, 59, 63, 80], [[4, 2, 3, 5, 6, 0, 1], 
						    [3, 2, 6, 1, 4, 0, 5]], 
						    [3, 6, 0, 4, 1, 2, 5], 
						    [3, 0, 6, 5, 4, 1, 2]],
		[[30, 59, 76, 99], [[4, 2, 3, 5, 6, 0, 1], 
						    [3, 2, 6, 1, 4, 0, 5]], 
						    [6, 2, 5, 0, 1, 4, 3], 
						    [2, 1, 6, 0, 3, 4, 5]],
		[[45, 99, 90, 91], [[6, 5, 1, 3, 2, 4, 0], 
						    [0, 5, 1, 3, 6, 4, 2]], 
						    [3, 5, 2, 4, 0, 6, 1], 
						    [5, 6, 4, 1, 3, 0, 2]],
		[[76, 99, 90, 91], [[6, 2, 5, 0, 1, 4, 3], 
						    [2, 1, 6, 0, 3, 4, 5]], 
						    [3, 5, 2, 4, 0, 6, 1], 
						    [5, 6, 4, 1, 3, 0, 2]]
	]

def prepareAncillaSteaneRandom(index, errorRates, name='a0', eigenstate='Z'):
	""" 
	"""
	indices, rounds = goodSteaneRandomPreps[index]
	ancillas = [permuteList(indices[i], rounds[i]) for i in range(4)]
	return prepareAncilla4(ancillas[0], ancillas[1], ancillas[2], ancillas[3], errorRates, name, eigenstate)
	
def prepareAndVerify(errorRates, XorZ, a0name, a0prep, a1name, a1prep, includeMeasRest=False):
	"""Given functions for preparing two ancillas, repeatedly call them to verify for XorZ-type errors."""
	attempts = 0
	while True:
		attempts += 1
		a0 = a0prep()
		a1 = a1prep()
		errors = {	'X':{a0name:a0['X'][a0name], a1name:a1['X'][a1name]}, 
					'Z':{a0name:a0['Z'][a0name], a1name:a1['Z'][a1name]}}	
		if XorZ == 'X':
			sourceblock, targetblock = a0name, a1name
			meas = measZ
		else:	# XorZ == 'Z'
			sourceblock, targetblock = a1name, a0name
			meas = measX
		for i in range(23):
			cnot(errorRates, errors, sourceblock, i, targetblock, i)
			meas(errorRates, errors, a1name, i)
			if includeMeasRest:
				rest(errorRates, errors, a0name, i)
		if not getSyndrome(errors[XorZ][a1name]):
			return errors, attempts
		# NOTE: The following commmented-out code is to study the case where we do not try to verify against Z
		# errors (this breaks the symmetry between preparing encoded |0> and |+>).  I used this code, instead of the 
		# above two lines, to study whether our analysis could study X errors completely separately from Z errors. 
		# Intuitively, the major problem with that approach is that Z-error verification eliminates many X errors 
		# as well (at least half of them, and perhaps as many as 3/4 for two-qubit failure locations).
		# For 10,000 trials at 2e-3 CNOT error rate, this code gave a logical X failure rate of 1.2e-3.  This compares 
		# to about 1e-4 to 4e-4 using X and Z verification.  While significant, if accurate, these numbers aren't
		# completely deadly.  A factor of 5 lost in the logical error rate corresponds to a factor of 5^{1/4} in the
		# physical error rate, or ~33%.  (Interestingly, logical X and logical Z errors almost never coincide [I
		# haven't ever seen it, anyway], so splitting up errors at the next level should be perfectly fine.)
		#if eigenstate is 'Z':	#### DEBUG TEMP
		#	if XorZ is 'Z':
		#		return errors	# I am not verifying against Z errors
		#	if XorZ is 'X':
		#		if not getSyndrome(errors[XorZ][a1name]):
		#			return errors
		#if eigenstate is 'X':
		#	if XorZ is 'X':
		#		return errors	# I am not verifying against Z errors
		#	if XorZ is 'Z':
		#		if not getSyndrome(errors[XorZ][a1name]):
		#			return errors	
	
def prepareAndVerifyNoPostselect(errorRates, XorZ, a0name, a0prep, a1name, a1prep, includeMeasRest=False):
	"""Given functions for preparing two ancillas, repeatedly call them to verify for XorZ-type errors."""
	attempts = 0
	while True:
		attempts += 1
		a0 = a0prep()
		a1 = a1prep()
		errors = {	'X':{a0name:a0['X'][a0name], a1name:a1['X'][a1name]}, 
					'Z':{a0name:a0['Z'][a0name], a1name:a1['Z'][a1name]}}	
		if XorZ == 'X':
			sourceblock, targetblock = a0name, a1name
			meas = measZ
		else:	# XorZ == 'Z'
			sourceblock, targetblock = a1name, a0name
			meas = measX
		for i in range(23):
			cnot(errorRates, errors, sourceblock, i, targetblock, i)
			meas(errorRates, errors, a1name, i)
			if includeMeasRest:
				rest(errorRates, errors, a0name, i)
		return errors, attempts


def prepareAncillaOverlap(errorRates, name='a0', eigenstate='Z'):
	"""This will simulate preparation of a verified ancilla.
	
	Note that it assumes that the error model is symmetrical under X and Z errors.  Thus to prepare an encoded 
	|+>, it prepares an encoded |0> and swaps the X and Z errors.  
	"""
	from golay import GolayOverlap
	
	perms = GolayOverlap.bestXZset
	cnotsA0 = GolayOverlap.getOverlapPrep(perms[0][0])
	cnotsA1 = GolayOverlap.getOverlapPrep(perms[0][1])
	cnotsA2 = GolayOverlap.getOverlapPrep(perms[1][0])
	cnotsA3 = GolayOverlap.getOverlapPrep(perms[1][1])

	return prepareAncilla4(cnotsA0, cnotsA1, cnotsA2, cnotsA3, errorRates, name, eigenstate)

def prepareAncilla4(cnotsA0, cnotsA1, cnotsA2, cnotsA3, errorRates, name='a0', eigenstate='Z'):
	"""This will simulate preparation of a verified ancilla.
	
	Note that it assumes that the error model is symmetrical under X and Z errors.  Thus to prepare an encoded 
	|+>, it prepares an encoded |0> and swaps the X and Z errors.  
	"""
	rounds = range(7)
	
	def prepA0():
		errors, attempts = prepareAndVerify(errorRates, 'X', 
			'a0', lambda: prepareUnverifiedAncilla(errorRates, rounds, cnotsA0, 'a0'), 
			'a1', lambda: prepareUnverifiedAncilla(errorRates, rounds, cnotsA1, 'a1'),
			includeMeasRest=True)
		prepStats['prepA0']['attempts'] += attempts
		prepStats['prepA0']['success'] += 1
		return errors
	def prepA2():
		errors, attempts = prepareAndVerify(errorRates, 'X', 
			'a2', lambda: prepareUnverifiedAncilla(errorRates, rounds, cnotsA2, 'a2'), 
			'a3', lambda: prepareUnverifiedAncilla(errorRates, rounds, cnotsA3, 'a3'),
			includeMeasRest=True)
		prepStats['prepA2']['attempts'] += attempts
		prepStats['prepA2']['success'] += 1
		return errors
		
	errors, attempts = prepareAndVerify(errorRates, 'Z', 'a0', prepA0, 'a2', prepA2, includeMeasRest=True)
	prepStats['prepA']['attempts'] += attempts
	prepStats['prepA']['success'] += 1
	
	if eigenstate == 'Z':	# might as well clean out the other ancilla blocks
		errors = {'X':{name:errors['X']['a0']}, 'Z':{name:errors['Z']['a0']}}
	else: #eigenstate == 'X' --- assuming that the error model is symmetrical, we can just swap the X and Z errors
		errors = {'X':{name:errors['Z']['a0']}, 'Z':{name:errors['X']['a0']}}
	#print "prepared ancilla:", util.counterUtils.errorsToStr(errors, 23, True)	
	return errors



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