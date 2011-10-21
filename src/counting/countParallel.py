'''
Functions and classes for parallel processing.
Created on 2010-08-30

@author: adam
'''
from counting.countErrors import convolveDict

print 'loading countParallel.'

from multiprocessing.pool import Pool
from util.counterUtils import convolveCounts, PartitionIterator
from util.listutils import addLists, addDicts
import logging
import math
import time


logger = logging.getLogger('countParallel')
defaultPool = None

class DummyPool(object):
	'''
	A Pool adapter that synchronizes all asynchronous calls.  Useful as a substitute for a real Pool
	when only 1 CPU is available (or Pool() would otherwise be initialized with only one worker process).
	'''


	def __init__(self, initializer=None, initargs=()):
		'''
		Constructor
		'''
		
		#self.pool = Pool(1, initializer, initargs)
		
	def apply(self, func, args=(), kwds={}):
		'''
		Equivalent of `apply()` builtin
		'''
		return self.pool.apply(func, args, kwds)

	def map(self, func, iterable, chunksize=None):
		'''
		Equivalent of `map()` builtin
		'''
		return self.pool.map(func, iterable, chunksize)
		
	def imap(self, func, iterable, chunksize=1):
		'''
		Equivalent of `itertools.imap()` -- can be MUCH slower than `Pool.map()`
		'''
		return self.pool.imap(func, iterable, chunksize)

	def imap_unordered(self, func, iterable, chunksize=1):
		'''
		Like `imap()` method but ordering of results is arbitrary
		'''
		return self.pool.imap_unordered(func, iterable, chunksize)

	def apply_async(self, func, args=(), kwds={}, callback=None):
		'''
		Asynchronous equivalent of `apply()` builtin
		'''
		result = apply(func, args, kwds)
		if None != callback:
			callback(result)
		return DummyResult(result)

	def map_async(self, func, iterable, chunksize=None, callback=None):
		'''
		Asynchronous equivalent of `map()` builtin
		'''
		result = self.pool.map(func, iterable, chunksize)
		if None != callback:
			callback(result)
		return DummyResult(result)


	def close(self):
		return self.pool.close()

	def terminate(self):
		return self.pool.terminate()
		
	def join(self):
		return self.pool.join()
	
class DummyResult(object):
	
	def __init__(self, result):
		self._result = result
		
	def get(self):
		return self._result
	
	def ready(self):
		return True


def enableFetchableMultiprocess():
	'''
	Code taken from: http://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods
	'''
	def _pickle_method(method):
		func_name = method.im_func.__name__
		obj = method.im_self
		cls = method.im_class
		return _unpickle_method, (func_name, obj, cls)
	
	def _unpickle_method(func_name, obj, cls):
		for cls in cls.mro():
			try:
				func = cls.__dict__[func_name]
			except KeyError:
				pass
			else:
				break
		return func.__get__(obj, cls)
	
	import copy_reg
	import types
	copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
	
def configureMultiProcess(nSlots):
	enableFetchableMultiprocess()
	
	logger.info('Configuring pool of {0} workers'.format(nSlots-1))
	if 1 >= nSlots:
		pool = DummyPool()
	else:
		pool = Pool(processes=(nSlots-1))
		
	setPool(pool)
	setSlots(nSlots)
		
	return pool

def setSlots(slots):
	global nSlots
	nSlots = slots
	
def numSlots():
	return nSlots

def getPool():
	if None == defaultPool:
		return Pool(1)
	return defaultPool

def setPool(pool):
	global defaultPool
	defaultPool = pool


convolve = lambda *args, **kwargs: convolveParallel(defaultPool, *args, **kwargs)

def asyncFuncWrapper(argList):
	asyncFcn = argList[0]
	return asyncFcn(*argList[1:])

# TODO: generalize to allow an arbitrary number of counts (instead of just 2).	
def convolveParallel(pool, counts0, counts1, partMax0=None, partMax1=None, kMax=None, convolveFcn=convolveDict, extraArgs=[], splitListsInto=[2,2]):
	
	if None == partMax0:
		partMax0 = len(counts0)-1
	if None == partMax1:
		partMax1 = len(counts1)-1
	if None == kMax:
		kMax = partMax0 + partMax1
		
	kMax = min(kMax, len(counts0) + len(counts1) - 2)
	
	resultsList = distributeConvolutions(pool, counts0, counts1, partMax0, partMax1, kMax, convolveFcn, extraArgs, splitListsInto)
			
	# Collate the results
	convolved = [0] * (kMax+1)

	sleeptime = 0
	while (len(resultsList)):
		time.sleep(sleeptime)
		sleeptime = 1
		
		done = [False] * len(resultsList)
		readyResults = [[] for _ in range(kMax+1)]
		for i, (k, result) in enumerate(resultsList):
			if result.ready():
				readyResults[k].append(result)
				done[i] = True
				
		# Remove the completed convolutions
		resultsList = [resultsList[i] for i in range(len(resultsList)) if not done[i]]
		
		resultsToProcess = sum(len(r) for r in readyResults)
		if resultsToProcess:
			logger.info('Now processing {0} results. {1} results remain to be convolved.'.format(resultsToProcess, len(resultsList)))
		
		for k in range(len(readyResults)):
			ready = readyResults.pop(0)
			if len(ready) == 0:
				continue
			
			ready = [r.get() for r in ready]
			if type(ready[0]) is dict:
				addFunc = addDicts
			else:
				addFunc = addLists
			
			logger.info('Summing contents of {0} lists of size {1} for k={2}'.format(len(ready), len(ready[0]), k))
			if 0 != convolved[k]:
				ready += [convolved[k]]
			
			convolved[k] = addFunc(*ready)
				

			
		
#	for k in range(len(resultsList)):
#		results = resultsList.pop(0)
#		results = results.get()
#		logger.info('Summing contents of {0} lists of size {1} for k={2}'.format(len(results), len(results[0]), k))
#		if type(results[0]) is dict:
#			convolved.append(addDicts(*results))
#		else:
#			convolved.append(addLists(*results))
#		logger.info('Completed parallel convolution for k={0}'.format(k))
	
	return convolved

def distributeConvolutions(pool, counts0, counts1, partMax0, partMax1, kMax, convolveFcn, extraArgs, splitListsInto):
	count0Parts = [splitContents(count, splitListsInto[0]) for count in counts0]
	count1Parts = [splitContents(count, splitListsInto[1]) for count in counts1]

	# Distribute the work
	results = []
	for k in reversed(range(kMax+1)):
		logger.info('Parallel Convolving {0} in {1}x{2} parts'.format(k, splitListsInto[0], splitListsInto[1]))
		
		for k0, k1 in PartitionIterator(k, 2, [partMax0, partMax1]):
			for count0Part in count0Parts[k0]:
				for count1Part in count1Parts[k1]:
					# Compute each convolution in parallel
					args = [count0Part, count1Part] + extraArgs
					result = pool.apply_async(convolveFcn, args)
					results.append((k,result))
		
	return results

def splitContents(listOrDict, numParts=2):
	if type(listOrDict) is dict:
		return splitDictContents(listOrDict, numParts)
	return splitListContents(listOrDict, numParts)

def splitDictContents(d, numParts=2):
	'''
	>>> splitDictContents({1:'a', 2:'b', 3:'c'}, 2)
	[{1: 'a', 2: 'b'}, {3: 'c'}]
	'''
	
	# Special case to avoid costly d.items() call below, if possible.
	if 1 == numParts:
		return [d]
	
	keys = d.keys()
	
	# If the number of elements in the dict is less than the
	# requested number of parts, then some of the parts will be empty.
	nonEmptyParts = min(len(keys), numParts)
	if 0 != nonEmptyParts:
		chunksize = int(math.ceil(len(keys) / float(nonEmptyParts)))
	else:
		chunksize = 1
		
	parts = []
	for i in xrange(0, len(keys), chunksize):
		end = min(i + chunksize, len(keys))
		parts.append({keys[j]: d[keys[j]] for j in xrange(i, end)})
		
	# Append any empty parts
	parts += [{}] * (numParts - nonEmptyParts)
			
	return parts

def splitListContents(l, numParts=2):
	'''
	Split the values of a list into the specified number of parts (slices).  Each slice
	is placed into a list of size equal to the original list, but for which the other
	entries have been set to zero.
	>>> l = [0,0,0,1,2,3,4,5,6,7,8,9,10]
	>>> splitListContents(l, 4)
	[[0,0,0,1,2,3,0,0,0,0,0,0,0], [0,0,0,0,0,0,4,5,6,0,0,0,0], [0,0,0,0,0,0,0,0,0,7,8,9,10]]
	'''
	if 1 == numParts:
		return [l]
	
	listLen = len(l)
	partLen = listLen / numParts
	
	parts = []
	for partNum in range(numParts):
		offset = partNum * partLen
		if partNum == (numParts-1):
			# This is the last part.  Be sure
			# to get the entire list
			end = listLen
		else:
			end = offset + partLen
		logger.info('Splitting into part {0} [{1}:{2}]'.format(partNum, offset, end))
		contents = l[offset:end]
		
		# If the contents of this part are all zeros, then it is useless.
		if any(contents):
			part = [0]*(offset) + contents + [0]*(listLen - end)
			parts.append(part)
		
	return parts


def packTasks(nSlots, tasks, costs):
	'''
	Assign tasks (somewhat) evenly among the given number of slots based on the given cost for
	each task. (Think 1-D bin packing, but without hard constraints.)
	'''
	assert len(tasks) == len(costs)
	
	#TODO: this technique is simple, but doesn't pack the slices very evenly
	totalCost = sum(costs)
	targetCost = totalCost / nSlots
	
	logger.info('totalCost={0}, targetCost={1}'.format(totalCost, targetCost))
	
	taskSlices = []
	for _ in range(nSlots):
		slice = []
		cost = 0
		while (cost < targetCost) and (0 < len(tasks)):
			c = tasks.pop()
			thisCost = costs.pop()
			if targetCost < cost + thisCost:
				tasks.append(c)
				costs.append(thisCost)
				break
			
			cost += thisCost
			slice.append(c)
		
		taskSlices.append(slice)
		
	# There may be a few left over
	#print 'leftovers=', len(tasks)

	i = 0
	for t in tasks:
		taskSlices[i].append(t)
		i = (i+1) % len(taskSlices)
		
	# Eliminate any empty slices
	taskSlices = filter(len, taskSlices)
		
	return taskSlices



def iterParallel(iterator, iterFunc, extraArgs=[], callback=None, pool=None):	
	if None == pool:
		pool = defaultPool
				
	results = [pool.apply_async(iterFunc, [i] + extraArgs, callback=callback) for i in iterator]
	return results


if __name__ == '__main__':
	import doctest
	doctest.testmod()
		





