'''
Created on 2010-03-12

@author: adam
'''
import gzip
import logging
import cPickle
import json as json
import os.path
import functools
import copy
import sys
import shelve

logger = logging.getLogger('util.cache')

fetchEnabled = True
memoEnabled = True

def enableFetch(enable=True):
	global fetchEnabled
	fetchEnabled = enable
	logger.info('Fetching enabled=' + str(fetchEnabled))
	
def enableMemo(enable=True):
	global memoEnabled
	memoEnabled = enable
	logger.info('Memos enabled=' + str(memoEnabled))


class memoize(object):
	'''
	Code copied from http://code.activestate.com/recipes/52201
	'''


	def __init__(self, function):
		'''
		Constructor
		'''
		self.func = function
		self.memo = {}
		
	def __call__(self, *args, **kwargs):
		key = self.get_key(args, kwargs)
		return self._fetch(key, args, kwargs)
	
	def _fetch(self, key, args, kwargs):
		if not memoEnabled:
			return self.func(*args, **kwargs)	

		try:
			return self.getMemo(key)
		except KeyError:
			self.setMemo(key, self.func(*args, **kwargs))
		return self.getMemo(key)
	
	def hasMemo(self, key):
		return self.memo.has_key(key)
	
	def setMemo(self, key, result):
		logger.debug('setting memo for {0}'.format(key))
		self.memo[key] = result
		
	def getMemo(self, key):
		logger.debug('getting memo for {0}'.format(key))
		return copy.copy(self.memo[key])
	
	def get_key(self, args, kwargs):
		key = [0] * (len(args) + len(kwargs))
		for i, arg in enumerate(args + tuple(kwargs.values())):
			try:
				key[i] = arg.__hash__()
			except TypeError:
				key[i] = str(arg)
		
		return tuple(key) + tuple(kwargs)
	
	def _methodCall(self, obj, *args, **kwargs):
		funcName = ''.join([repr(obj), '.', self.func.func_name])
		key = self.get_key(tuple([funcName]) + args, kwargs)
		return self._fetch(key, tuple([obj]) + args, kwargs)
	
	def __get__(self, obj, objtype):
		'''
		Defining __get__ causes memoize to become a descriptor, which means
		that, when @memoize is applied to an instance method __get__ will
		be called.  The return value of __get__ is then called to obtain the
		function value.
		In this case (when the function is an instance method) we also
		want to include the instance data in the fetch key.  This is accomplished
		by calling _methodCall with the appropriate object.
		'''
		return functools.partial(self._methodCall, obj)
	
class memoizeMutable(memoize):
	
	def __init__(self, function):
		super(memoizeMutable, self).__init__(function)
	
	def get_key(self, args):
		return cPickle.dumps(args)
	
class memoizeFetchable(memoize):
	def __init__(self, function):
		super(memoizeFetchable, self).__init__(function)
		self.dm = DataManager()
		
	def hasMemo(self, key):
		if super(memoizeFetchable, self).hasMemo(key):
			return True
		
		return self.dm.exists(self.getDMKey(key))
	
	def getDMKey(self, key):
		argStr = reduce(lambda s, arg: s + '.' + str(arg), key, '')
		return self.func.func_name + argStr
	
	def setMemo(self, key, result):
		self.dm.save(result, self.getDMKey(key))
		return super(memoizeFetchable, self).setMemo(key, result)
	
	def getMemo(self, key):
		if not super(memoizeFetchable, self).hasMemo(key):
			result = self.dm.load(self.getDMKey(key))
			super(memoizeFetchable, self).setMemo(key, result)
			
		return super(memoizeFetchable, self).getMemo(key)
		
		
	

		
class fetchable(object):
	'''
	Decorator that automatically loads the return value of the function from a file, when available.
	If a file containing the return value does not exist, the function is exectuted normally and
	the return value is saved to a file.  The filename is based on the function name and the
	arguments passed to the function.
	'''
	
	def __init__(self, func):
		self.func = func
		self.fetched = False
		#self.dm = DataManager()
	
	def __call__(self, *args, **kwargs):
		key = self.get_key(self.func.func_name, args, kwargs)
		return self._fetch(key, args, kwargs)
	
	def _methodCall(self, obj, *args, **kwargs):
		funcName = repr(obj) + '.' + self.func.func_name
		key = self.get_key(funcName, args, kwargs)
		return self._fetch(key, tuple([obj]) + args, kwargs)
	
	def _fetch(self, key, args, kwargs):
		if not fetchEnabled:
			return self.func(*args, **kwargs)
			
		dm = DataManager()	
		try:
			data = dm.load(key)
			self.fetched = True
		except IOError:
			logger.debug('Fetch of {0} failed. Computing from scratch.'.format(key))
			data = self.func(*args, **kwargs)
			dm.save(data, key)
				
		return data

	
	@staticmethod
	def get_key(funcName, args, kwargs):
		args = list(args) + kwargs.values()
		key =  reduce(lambda s, arg: s + '.' + str(arg), args, funcName)
		
		# TODO: create a proper string filter class.
		badChars = [' ', '(', ')', '[', ']', '{', '}', '\'', '|', '>']
		for c in badChars:
			key = key.replace(c,'')
			
		# Windows (NTFS) does not allow colons.	
		key = key.replace(':', ';')
		
		return key
		
#	def fetch(self, key):
#		data = self.dm.load(key)
#		self.fetched = True
#		return data
	
	def __get__(self, obj, objtype):
		'''
		Defining __get__ causes fetchable to become a descriptor, which means
		that, when @fetchable is applied to an instance method __get__ will
		be called.  The return value of __get__ is then called to obtain the
		function value.
		In this case (when the function is an instance method) we also
		want to include the instance data in the fetch key.  This is accomplished
		by calling _methodCall with the appropriate object.
		'''
		return functools.partial(self._methodCall, obj)
		
	


class DataManager(object):
	'''
	Abstraction for saving and reading files.
	'''	
	
	# TODO: parameterize the default data dir.
	def __init__(self, dataDir = (os.path.pardir + os.path.sep + 'data' + os.path.sep)):
		'''
		Constructor
		'''		
		self.dataDir = dataDir
		self.fileExt = '.txt.gz'
		
		self.save = self.savepickle
		self.load = self.loadpickle
		
		if not os.path.exists(self.dataDir):
			os.mkdir(self.dataDir)
			
		self.lookup = self.initializeLookup()
			
	def initializeLookup(self):
		filename = self.dataDir + 'datafile-lookup.txt'
		lookup = shelve.open(filename)
		return lookup
			
#	def savejson(self, obj, key=None):
#		if None == key:
#			key = self.constructKey(obj)
#			
#		filename = self.constructFilename(key)
##		logger.debug('Saving {0} to {1}'.format(key, filename))
#		outfile = gzip.open(filename, 'wb', 1)
#		outfile.write(json.dumps(obj))
#		outfile.close()
#		
#		logger.debug('Saved {0} to {1}'.format(key, filename))
#	
#		return key
	
	def savepickle(self, obj, key=None):
		if None == key:
			key = self.constructKey(obj)
			
		filename = self.constructFilename(key)
#		logger.debug('Saving {0} to {1}'.format(key, filename))
		outfile = gzip.open(filename, 'wb', 1)
		outfile.write(cPickle.dumps(obj, 2))
		outfile.close()
		
		logger.debug('Saved {0} to {1}'.format(key, filename))
	
		return key
	
	def constructKey(self, obj):
		return str(obj)
	
	def constructFilename(self, key):
		try:
			lookupVal = self.lookup[key]
		except KeyError:
			lookupVal = len(self.lookup)
			self.lookup[key] = lookupVal
			# Write the new entry to disk
			self.lookup.sync() 
		
		return self.dataDir + 'datafile.' + str(lookupVal) + self.fileExt
	
#	def loadjson(self, key):
#		filename = self.constructFilename(key)
##		logger.debug('Loading {0} from {1}'.format(key, filename))
#		infile = gzip.open(filename, 'rb')
#		obj = json.loads(infile.read())
#		infile.close()
#		
#		logger.debug('Loaded {0} from {1}'.format(key, filename))
#		
#		return obj
	
	def loadpickle(self, key):
		filename = self.constructFilename(key)
#		logger.debug('Loading {0} from {1}'.format(key, filename))
		infile = gzip.open(filename, 'rb')
		obj = cPickle.loads(infile.read())
		infile.close()
		
		logger.debug('Loaded {0} from {1}'.format(key, filename))
		
		return obj
	
	def exists(self, key):
		filename = self.constructFilename(key)
		return os.path.exists(filename)
	
if __name__ == '__main__':
	
	@fetchable
	def fetchFunc(a1, a2, a3='arg3'):
		print a1, a2, a3
	
	fetchFunc(1,2)
		  
		
		
		