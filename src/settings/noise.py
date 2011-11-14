'''
Class representations of noise models.

Created on Nov 7, 2010

@author: adam
'''
from abc import abstractmethod, ABCMeta
from qec.error import Pauli
from util.polynomial import sympoly1d, SymPolyWrapper
from util.iteration import SubsetIterator
import warnings

class Bound(object):
	UpperBound = 0
	LowerBound = 1


errorListX = {
		'prepZ': (Pauli.X,),
		'measZ': (Pauli.X,),	
		'rest':  (Pauli.X,),
		'cnot':  (Pauli.I+Pauli.X, Pauli.X+Pauli.I, Pauli.X+Pauli.X)
}

errorListZ = {
		'prepX': (Pauli.Z,),
		'measX': (Pauli.Z,),	
		'rest':  (Pauli.Z,),
		'cnot':  (Pauli.I+Pauli.Z, Pauli.Z+Pauli.I, Pauli.Z+Pauli.Z)
}

errorListXZ = {
		'prepZ': (Pauli.X,),
		'measZ': (Pauli.X,),
		'prepX': (Pauli.Z,),
		'measX': (Pauli.Z,),	
		'rest':  (Pauli.X, Pauli.Z, Pauli.Y),
		'cnot':  (Pauli.I+Pauli.X,
				  Pauli.I+Pauli.Z,
				  Pauli.I+Pauli.Y,
				  Pauli.X+Pauli.I,
				  Pauli.X+Pauli.X,
				  Pauli.X+Pauli.Z,
				  Pauli.X+Pauli.Y,
				  Pauli.Z+Pauli.I,
				  Pauli.Z+Pauli.X,
				  Pauli.Z+Pauli.Z,
				  Pauli.Z+Pauli.Y,
				  Pauli.Y+Pauli.I,
				  Pauli.Y+Pauli.X,
				  Pauli.Y+Pauli.Z,
				  Pauli.Y+Pauli.Y)
}			


class NoiseModel(object):
	'''
	Abstract base class for concrete noise models.  
	Initilizer takes a noise strength interval [gMin,gMax] as arguments.
	'''
	__metaclass__ = ABCMeta
	
	def __init__(self, gMin=0, gMax=1):
		self._gMin = gMin
		self._gMax = gMax
				
	
	def numErrors(self, loc):
		'''
		Returns the number of different errors that can occur
		at location loc.
		
		DEPRECATED
		Use errorList(), instead.
		'''
		warnings.warn('Use errorList(), instead.', category=DeprecationWarning)
		return len(self.errorList(loc))
	
	@abstractmethod
	def errorList(self, loc):
		'''
		Returns a list of errors that can occur (with non-zero probability)
		at the given location.
		'''

	def noiseRange(self):
		'''
		Returns the interval [gMin,gMax] over which the noise model is defined.
		'''
		return (self._gMin, self._gMax)
	
	def __str__(self):
		return '[{0}, {1}]'.format(self._gMin, self._gMax)
					
	def __repr__(self):
		return str(self)
	
class DepolarizingNoiseModel(NoiseModel):
	'''
	Abstract base class for depolarizing noise models.
	'''
	
	def __init__(self, gMin=0, gMax=1):
		super(DepolarizingNoiseModel,self).__init__(gMin,gMax)
	
	@abstractmethod
	def getWeight(self, loc, error, bound):
		'''
		Returns the "weight" associated with the given error for
		location loc.  The probability of the error is calculated
		as Pr[error] = weight * gamma, where gamma is the noise
		strength.
		'''
			
	@abstractmethod
	def prFail(self, loc, bound):
		'''
		Returns the probability that location loc will fail, expressed
		in terms of the noise strength.
		'''
	
	def prIdeal(self, loc, bound):
		'''
		Returns the probability that location loc acts ideally,
		i.e, that it does not fail.
		'''
		if Bound.UpperBound == bound:
			bound = Bound.LowerBound
		else:
			bound = Bound.UpperBound
			
		return 1 - self.prFail(loc, bound)
		
class DepolarizingNoiseModelSympy(DepolarizingNoiseModel):
	'''
	Abstract sub-class of depolarizing noise in which probabilities are
	represented by Sympy symbolic expressions.
	'''
							
	def prFail(self, loc, bound):
		weights = [self.getWeight(loc, error, bound) for error in self.errorList(loc)]
		coeffs = [sum(weights), 0]
		sympoly = sympoly1d(coeffs)		
		return SymPolyWrapper(sympoly)
	
class CountingNoiseModel(DepolarizingNoiseModelSympy):
	
	def getWeight(self, loc, error, bound=Bound.UpperBound):
		return 1
	
	def likelyhood(self, bound=Bound.UpperBound):
		# 1 / (1-g)
		return SymPolyWrapper(1 / sympoly1d([-1, 1]))
	
class CountingNoiseModelX(CountingNoiseModel):
	
	def errorList(self, loc):
		try:
			return errorListX[loc['type']]
		except KeyError:
			return []
		
class CountingNoiseModelZ(CountingNoiseModel):
	
	def errorList(self, loc):
		try:
			return errorListZ[loc['type']]
		except KeyError:
			return []

		
#class UpperBoundNoiseModelSympy(DepolarizingNoiseModelSympy):
#	'''
#	Special case of depolarizing noise in which upper bounds on the
#	error probabilities are known, but an upper bound on the
#	probability that a location does *not* fail is unknown.
#	'''
#	
#	def prIdeal(self, loc):
#		return SymPolyWrapper(sympoly1d([1]))
	
class NoiseModelMarginalSympy(DepolarizingNoiseModelSympy):
	'''
	Marginal noise model for considering X errors (or alternatively Z errors)
	independently of Z errors (resp. X errors).
	'''
	
	def __init__(self, errorList, gMin=0, gMax=1):
		super(NoiseModelMarginalSympy, self).__init__(gMin, gMax)
		self._errorList = errorList
	
	def getWeight(self, loc, error, bound=Bound.UpperBound):
		# Weights are the same regardless of the bound type.
		if loc['type'] == 'rest':
			return 8
		return 4
	
	def likelyhood(self, bound):
		# g/(1-12g): Upper bound by dividing by 1-12g in all cases.
		if Bound.UpperBound == bound:
			return  SymPolyWrapper(sympoly1d([1,0]) / sympoly1d([-12, 1]))
		
		# g/(1-4g): Upper bound by dividing by 1-4g in all cases.
		return  SymPolyWrapper(sympoly1d([1,0]) / sympoly1d([-4, 1]))
	
	def errorList(self, loc):
		try:
			return self._errorList[loc['type']]
		except KeyError:
			return []
	
	def __str__(self):
		return 'w=4.r=8'
	
class NoiseModelXSympy(NoiseModelMarginalSympy):
	
	def __init__(self, gMin=0, gMax=1):
		super(NoiseModelXSympy, self).__init__(errorListX, gMin, gMax)
	
class NoiseModelZSympy(NoiseModelMarginalSympy):
	
	def __init__(self, gMin=0, gMax=1):
		super(NoiseModelZSympy, self).__init__(errorListZ, gMin, gMax)
	
	
class NoiseModelXZSympy(DepolarizingNoiseModelSympy):
	'''
	Concrete realization of the depolarizing noise model in which
	noise weights, likelyhoods, etc. are given in terms of *lower* bounds.
	'''
	
	def getWeight(self, loc, error, bound=Bound.UpperBound):
		'''
		For XZ counts we want lower bounds and so the likelyhoods are calculated as g/(1-4g) instead
		of g/(1-15g).  As a result, the weights are a bit different.
		'''
		
		# Weights are the same regardless of bound type
		if loc['type'] == 'cnot':
			return 1
		return 4
	
	def likelyhood(self, bound):
		if Bound.LowerBound == bound:
			# g/(1-4g): Lower bound by dividing by 1-4g in all cases.
			return SymPolyWrapper(sympoly1d([1,0]) / sympoly1d([-4, 1]))
		
		# g/(1-15g): Upper bound by dividing by 1-15g in all cases.
		return SymPolyWrapper(sympoly1d([1,0]) / sympoly1d([-15, 1]))
	
	def errorList(self, loc):
		try:
			return errorListXZ[loc['type']]
		except KeyError:
			return []
	
	def __str__(self):
		return 'lower-w=4.c=1'

class TransformedNoiseModelSympy(DepolarizingNoiseModelSympy):
	'''
	Abstract class for a transformed noise model, i.e., one
	in which error weights of similar errors may be drastically
	different.
	'''
			
	def getWeight(self, loc, error):
		return self._weights[loc['type']][error]
	
	def errorList(self, loc):
		return self._weights[loc['type']].keys()
	
	def likelyhood(self):
		return SymPolyWrapper(sympoly1d([1,0]))
	
	def __str__(self):
		s = ''
		for w in self._weights.values():
			s += str(w) + '.'
			
		return s[:-1]

	
class TransformedNoiseModelXSympy(TransformedNoiseModelSympy):
	'''
	Marginal transformed error model for X errors.
	'''
	
	def __init__(self, prepZ, measZ, rest, cnotIX, cnotXI, cnotXX, gMin, gMax):
		super(TransformedNoiseModelXSympy, self).__init__(gMin, gMax)
		self._weights = {'prepZ': {Pauli.X: prepZ},
						 'measZ': {Pauli.X: measZ},
						 'rest':  {Pauli.X: rest},
						 'cnot':  {Pauli.I+Pauli.X: cnotIX, 
								   Pauli.X+Pauli.I: cnotXI, 
								   Pauli.X+Pauli.X: cnotXX}
						}

class TransformedNoiseModelZSympy(TransformedNoiseModelSympy):
	'''
	Marginal transformed error model for Z errors.
	'''
	
	def __init__(self, prepX, measX, rest, cnotIZ, cnotZI, cnotZZ, gMin, gMax):
		super(TransformedNoiseModelZSympy, self).__init__(gMin, gMax)
		self._weights = {'prepX': {Pauli.Z: prepX},
						 'measX': {Pauli.Z: measX},
						 'rest':  {Pauli.Z: rest},
						 'cnot':  {Pauli.I+Pauli.Z: cnotIZ, 
								   Pauli.Z+Pauli.I: cnotZI, 
								   Pauli.Z+Pauli.Z: cnotZZ}}
		
			
class NoiseModelZeroSympy(DepolarizingNoiseModelSympy):
	'''
	Special case in which the noise strength is zero.
	'''
	
	def getWeight(self, loc, error):
		return 0
	
	def errorList(self, loc):
		return []
	
	def likelyhood(self):
		return SymPolyWrapper(sympoly1d([0]))
	
	def prIdeal(self, loc):
		return SymPolyWrapper(sympoly1d([1]))
	
	def __str__(self):
		return '0'
		

class FixedPointNoise(NoiseModel):
	'''
	Special case noise model in which the noise strength is fixed
	a particular value.
	Constructed as FixedPointNoise(noiseModel, point), where noiseModel
	is the underlying noise model to use, and point is the value
	at which the noise strength is fixed.
	'''
	
	def __init__(self, noiseModel, point):
		super(FixedPointNoise, self).__init__(point,point)
		self._noise = noiseModel
		self._point = point
		
	def getWeight(self, loc, error):
		return self._noise.getWeight(loc, error)
	
	def likelyhood(self):
		return self._noise.likelyhood()(self._point)
	
	def numErrors(self, loc):
		return self._noise.numErrors(loc)
				
	def prFail(self, loc):
		return self._noise.prFail(loc)(self._point)
	
	def prIdeal(self, loc):
		return self._noise.prIdeal(loc)(self._point)