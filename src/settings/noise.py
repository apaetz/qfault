'''
Class representations of noise models.

Created on Nov 7, 2010

@author: adam
'''
from util.counterUtils import errorRangeX, errorRangeXZ
from abc import abstractmethod, ABCMeta
from util.polynomial import sympoly1d, SymPolyWrapper


class NoiseModel(object):
	'''
	Abstract base class for concrete noise models.  
	Initilizer takes a noise strength interval [gMin,gMax] as arguments.
	'''
	__metaclass__ = ABCMeta
	
	def __init__(self, gMin=0, gMax=1):
		self._gMin = gMin
		self._gMax = gMax
				
	@abstractmethod
	def numErrors(self, loc):
		'''
		Returns the number of different errors that can occur
		at location loc.
		'''

	def noiseRange(self):
		'''
		Returns the interval [gMin,gMax] over which the noise model is defined.
		'''
		return (self._gMin, self._gMax)
					
	def __repr__(self):
		return str(self)
	
class DepolarizingNoiseModel(NoiseModel):
	'''
	Abstract base class for depolarizing noise models.
	'''
	
	def __init__(self, gMin=0, gMax=1):
		super(DepolarizingNoiseModel,self).__init__(gMin,gMax)
	
	@abstractmethod
	def getWeight(self, loc, error):
		'''
		Returns the "weight" associated with the given error for
		location loc.  The probability of the error is calculated
		as Pr[error] = weight * gamma, where gamma is the noise
		strength.
		'''
	
	@abstractmethod
	def prFail(self, loc):
		'''
		Returns the probability that location loc will fail, expressed
		in terms of the noise strength.
		'''
	
	def prIdeal(self, loc):
		'''
		Returns the probability that location loc acts ideally,
		i.e, that it does not fail.
		'''
		return 1 - self.prFail(loc)
		
class DepolarizingNoiseModelSympy(DepolarizingNoiseModel):
	'''
	Abstract sub-class of depolarizing noise in which probabilities are
	represented by Sympy symbolic expressions.
	'''
							
	def prFail(self, loc):
		weights = [self.getWeight(loc, error) for error in range(self.numErrors(loc))]
		coeffs = [sum(weights), 0]
		sympoly = sympoly1d(coeffs)		
		return SymPolyWrapper(sympoly)
		
class UpperBoundNoiseModelSympy(DepolarizingNoiseModelSympy):
	'''
	Special case of depolarizing noise in which upper bounds on the
	error probabilities are known, but an upper bound on the
	probability that a location does *not* fail is unknown.
	'''
	
	def prIdeal(self, loc):
		return SymPolyWrapper(sympoly1d([1]))
	
class NoiseModelXSympy(DepolarizingNoiseModelSympy):
	'''
	Marginal noise model for considering X errors (or alternatively Z errors)
	independently of Z errors (resp. X errors).
	'''
	
	def getWeight(self, loc, error):
		if loc['type'] == 'rest':
			return 8
		return 4
	
	def likelyhood(self):
		# g/(1-12g): Upper bound by dividing by 1-12g in all cases.
		return  SymPolyWrapper(sympoly1d([1,0]) / sympoly1d([-12, 1]))
	
	def numErrors(self, loc):
		return errorRangeX(loc)
	
	def __str__(self):
		return 'w=4.r=8'
	
	
class NoiseModelXZSympy(DepolarizingNoiseModelSympy):
	'''
	Concrete realization of the depolarizing noise model in which
	noise weights, likelyhoods, etc. are given in terms of *lower* bounds.
	'''
	
	def getWeight(self, loc, error):
		'''
		For XZ counts we want lower bounds and so the likelyhoods are calculated as g/(1-4g) instead
		of g/(1-12g).  As a result, the weights are a bit different.
		'''
		if loc['type'] == 'cnot':
			return 1
		return 4
	
	def likelyhood(self):
		# g/(1-4g): Lower bound by dividing by 1-4g in all cases.
		return SymPolyWrapper(sympoly1d([1,0]) / sympoly1d([-4, 1]))
	
	def numErrors(self, loc):
		return errorRangeXZ(loc)
	
	def __str__(self):
		return 'w=4.c=1'

class TransformedNoiseModelSympy(UpperBoundNoiseModelSympy):
	'''
	Abstract class for a transformed noise model, i.e., one
	in which error weights of similar errors may be drastically
	different.
	'''
			
	def getWeight(self, loc, error):
		return self._weights[loc['type']][error]
	
	def numErrors(self, loc):
		return errorRangeX(loc)
	
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
		self._weights = {'prepZ': [prepZ],
						 'measZ': [measZ],
						 'rest':  [rest],
						 'cnot':  [cnotIX, cnotXI, cnotXX]}

class TransformedNoiseModelZSympy(TransformedNoiseModelSympy):
	'''
	Marginal transformed error model for Z errors.
	'''
	
	def __init__(self, prepX, measX, rest, cnotIZ, cnotZI, cnotZZ, gMin, gMax):
		super(TransformedNoiseModelZSympy, self).__init__(gMin, gMax)
		self._weights = {'prepX': [prepX],
						 'measX': [measX],
						 'rest':  [rest],
						 'cnot':  [cnotIZ, cnotZI, cnotZZ]}
		
			
class NoiseModelZeroSympy(DepolarizingNoiseModelSympy):
	'''
	Special case in which the noise strength is zero.
	'''
	
	def getWeight(self, loc, error):
		return 0
	
	def numErrors(self, loc):
		return 0
	
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