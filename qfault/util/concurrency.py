'''
Functions and classes for parallel processing.
Created on 2010-08-30

@author: adam
'''
import multiprocessing
import logging
from qfault.util import listutils
import itertools

__all__ = ['enable_concurrency', 'map_concurrent']


logger = logging.getLogger('count_parallel')

def initialize_concurrency(number_of_slots=multiprocessing.cpu_count()):
    '''
    Initialize the concurrency module with the given number of 
    slots (processes).  By default, it sets the number of slots
    to equal the number of CPUs.
    It is intended that this function be called once, at the
    beginning of execution.  However, it may be called any
    number of times to reset the number of slots.
    '''
    global total_slots
    
    _enable_concurrent_pickle()
    
    if 0 == number_of_slots:
        pool =DummyPool()
    else:
        logger.info('Configuring pool of {0} workers'.format(number_of_slots))
        pool = multiprocessing.Pool(processes=(number_of_slots))
        
    _set_pool(pool)
    _set_slot_count(max(1,number_of_slots))

class Noop(object):
    
    def __call__(self, arg):
        return arg
    
def map_concurrent(function, collection):
    '''
    Concurrently maps all of the elements in the 
    collection according to 'function'. Behavior
    is equivalent to map(collection) but processing
    is done concurrently on slices of the input.
    
    >>> initialize_concurrency(1)
    >>> import operator
    >>> import functools
    >>> negate = functools.partial(operator.mul, -1)
    >>> list(map_concurrent(negate, range(4)))
    [0, -1, -2, -3]
    '''
    mapped = mapreduce_concurrent(function, Noop(), collection)
    return itertools.chain(*mapped)

class SliceMapReduce(object):
    
    def __init__(self, map_func, reduce_func):
        self._map = map_func
        self._reduce = reduce_func
        
    def __call__(self, _slice):
        mapped = map(self._map, _slice)
        return self._reduce(mapped)

def mapreduce_concurrent(map_func, reduce_func, collection):
    '''
    Concurrently maps all of the elements in the 
    collection according to 'map_func', then reduces
    the result according to reduce_func. Behavior
    is equivalent to reduce(map(collection)), but
    the processing is done concurrently on slices of the
    input.
    
    >>> initialize_concurrency(1)
    >>> import operator
    >>> import functools
    >>> negate = functools.partial(operator.mul, -1)
    >>> mapreduce_concurrent(negate, sum, range(4))
    -6
    '''   
    slice_map = SliceMapReduce(map_func, reduce_func)
    
    pool = _get_pool()
    map_result = pool.map(slice_map,
                          listutils.equal_chop(collection, 
                                               _slot_count()))
    return reduce_func(map_result)
    

def _enable_concurrent_pickle():
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
    
def _set_slot_count(count):
    global _num_slots
    _num_slots = count
    
def _slot_count():
    global _num_slots
    return _num_slots

def _get_pool():
    global _the_pool
    return _global_pool

def _set_pool(pool):
    global _global_pool
    _global_pool = pool
    

class DummyPool(object):
    '''
    '''

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
        result = map(func, iterable, chunksize)
        if None != callback:
            callback(result)
        return DummyResult(result)
    
    def map(self, func, iterable, chunksize=None):
        '''
        Equivalent of `map()` builtin
        '''
        return map(func, iterable)

class DummyResult(object):
    
    def __init__(self, result):
        self._result = result
        
    def get(self):
        return self._result
    
    def ready(self):
        return True
    
    
_global_pool = DummyPool()
_num_slots = 1

if __name__ == '__main__':
    import doctest
    doctest.testmod()
        


