from util.counterUtils import sumLocs
import operator
from util import counterUtils


class Locations(object):
	'''
	Container for a set of circuit locations.  The primary motivation for this
	class is to allow a name to be attached to the list of locations.
	'''
	
	def __init__(self, list, name):
		self.list = list
		self.name = name
		
	def __str__(self):
		return self.name
	
	def __repr__(self):
		return str(self)

	def __getitem__(self, index):
		return self.list[index]
	
	def __len__(self):
		return len(self.list)
	
	def __add__(self, other):
		return Locations(self.list + other.list, str(self) + ':' + str(other))
			
	def filterAgainst(self, type):
		return Locations(filter(lambda x: x['type'] != type, self.list), self.name + '-' + type)
		
	def getTotals(self):
		return LocationCount._fromDict(sumLocs(self))
	
	def blocknames(self):
		return counterUtils.allBlocks(self)
	
	
class LocationCount(tuple):
		'LocationCount(cnot, prepX, prepZ, measX, measZ, rest)' 

		__slots__ = () 

		_fields = ('cnot', 'prepX', 'prepZ', 'measX', 'measZ', 'rest') 

		def __new__(cls, cnot, prepX, prepZ, measX, measZ, rest):
			return tuple.__new__(cls, (cnot, prepX, prepZ, measX, measZ, rest)) 

		@classmethod
		def _make(cls, iterable, new=tuple.__new__, len=len):
			'Make a new LocationCount object from a sequence or iterable'
			result = new(cls, iterable)
			if len(result) != 6:
				raise TypeError('Expected 6 arguments, got %d' % len(result))
			return result 
		
		@classmethod
		def _fromDict(cls, dict):
			cnot = dict.get('cnot', 0)
			prepX = dict.get('prepX', 0)
			prepZ = dict.get('prepZ', 0)
			measX = dict.get('measX', 0)
			measZ = dict.get('measZ', 0)
			rest = dict.get('rest', 0)
			
			return LocationCount(cnot, prepX, prepZ, measX, measZ, rest)
			

		def __repr__(self):
			return 'LocationCount(cnot=%r, prepX=%r, prepZ=%r, measX=%r, measZ=%r, rest=%r)' % self 

		def _asdict(self):
			'Return a new dict which maps field names to their values'
			return {'cnot': self[0], 'prepX': self[1], 'prepZ': self[2], 'measX': self[3], 'measZ': self[4], 'rest': self[5]} 

		def _replace(self, **kwds):
			'Return a new LocationCount object replacing specified fields with new values'
			result = self._make(map(kwds.pop, ('cnot', 'prepX', 'prepZ', 'measX', 'measZ', 'rest'), self))
			if kwds:
				raise ValueError('Got unexpected field names: %r' % kwds.keys())
			return result 

		def __getnewargs__(self):
			return tuple(self) 

		cnot = property(operator.itemgetter(0))
		prepX = property(operator.itemgetter(1))
		prepZ = property(operator.itemgetter(2))
		measX = property(operator.itemgetter(3))
		measZ = property(operator.itemgetter(4))
		rest = property(operator.itemgetter(5))


		def __mul__(self, scalar):
			new = [scalar * c for c in self]
			return LocationCount(*new)

		def __rmul__(self, scalar):
			return self * scalar

		def __add__(self, other):
			newList = [self[i] + other[i] for i in range(len(self))]
			return LocationCount(*newList)
			
		def __sub__(self, other):
			newList = [self[i] - other[i] for i in range(len(self))]
			return LocationCount(*newList)			
		
		def filterOut(self, str):
			new = LocationCount(*self)
			if 'X' == str:
				new.prepX = 0
				new.measX = 0
			if 'Z' == str:
				new.prepZ = 0
				new.measZ = 0
			return new
		
		def dual(self):
			'''
			Returns location counts for the equivalent dual circuit.
			(i.e. swaps prepX <-> prepZ, measX <-> measZ)
			'''
			return LocationCount(self.cnot, self.prepZ, self.prepX, self.measZ, self.measX, self.rest)


