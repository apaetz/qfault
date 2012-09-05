'''
Created on 2011-08-29

@author: adam
'''
import math
import itertools


### NOTE: There are built-in itertools functions that do this kind of thing (http://docs.python.org/py3k/library/itertools.html#itertools.combinations_with_replacement)
class SubsetIterator:
    """Class for iterating over all subsets of a given size of (a shallow copy of) the given list.
    
    The subsets are given in lexicographic order of their indices in range(len(list)).
    >>> [t for t in SubsetIterator(range(1,5),2)]
    [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]
    """
    def __init__(self, universe, subsetsize):
        self.universe = universe[:]
        self.subsetsize = subsetsize
    def __iter__(self):
        if (0 == self.subsetsize) or (len(self.universe) < self.subsetsize):
            return [[]].__iter__()
        if self.subsetsize == len(self.universe):
            return [self.universe].__iter__()
        
        self.a = range(self.subsetsize)    # a is used to store the indices of the subset
        self.a[-1] -= 1        # decrement the last index so first increment will return first subset
        return self
    def next(self):
        n = len(self.universe)
        r = self.subsetsize
        a = self.a
        if a[0] == n-r:
            raise StopIteration        
        a[-1] = a[-1] + 1
        if a[-1] > n-1:
            j = r-2
            while a[j] == n-1 - (r-1-j):
                j = j-1
            a[j] = a[j] + 1
            for i in range(j+1, r):
                a[i] = a[j] + i - j
        return [self.universe[x] for x in a]    

class TupleIterator:
    """Class for iterating over all lists a where a[j] in range(0,ranges[j]).
    
    >>> [t for t in TupleIterator([1,2,3])]
    [[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 1, 0], [0, 1, 1], [0, 1, 2]]
    """
    def __init__(self, ranges):
        self.ranges = ranges
    def __iter__(self):
        if 0 == len(self.ranges):
            return ().__iter__()
        
        self.a = [0] * len(self.ranges)
        self.a[-1] = -1        # set the last element to -1, so the first call to next returns [0,0,..,0]
        return self
    def next(self):
        a = self.a
        ranges = self.ranges
        r = len(a)
        j = r-1
        while a[j] == ranges[j] - 1:
            j = j - 1
            if j < 0:
                a = [0] * len(ranges)
                raise StopIteration
        a[j] = a[j] + 1
        a[j+1:] = [0] * (r-j-1)
        return a[:]            # we return a copy of the internal list (not necessary for this application)
    
    
class PartitionIterator:
    '''
    Iterates over all possible assignments of of 'whole' objects into 'numParts' parts.
    Returns a list of size numParts during each iteration.  Each element in this
    list specifies the portion of the whole allotted to that element.  The sum of
    all of the elements in the list will always be equal to whole.
    The optional maxInParts argument is a list that specifies the maximum portion
    that can be assigned to each part.  For example, if maxInParts[i]=k, the value of the
    i'th element of the list returned for each iteration will be at most k.
    
    >>> [partition for partition in PartitionIterator(3,3, [3,1,2], [1,0,0])]
    [[3, 0, 0], [2, 1, 0], [2, 0, 1], [1, 1, 1], [1, 0, 2]]
    '''
    
    def __init__(self, whole, numParts, maxInParts=None, minInParts=None):
        self.whole = whole
        self.numParts = numParts
        if maxInParts == None:
            self.maxInParts = [whole] * numParts
        else:
            self.maxInParts = maxInParts
        assert numParts <= len(self.maxInParts)
        
        if minInParts == None:
            self.minInParts = [0] * numParts
        else:
            self.minInParts = minInParts
        assert numParts <= len(self.minInParts)
    
    def __iter__(self):
        if self.numParts == 0:
            return [].__iter__()
        
        self.leadPart = min(self.whole, self.maxInParts[0])
        
        if self.leadPart < self.minInParts[0]:
            return [].__iter__()
        
        if self.numParts == 1:
            if (self.whole == self.leadPart) and (self.whole >= self.minInParts[0]):
                return [[self.leadPart]].__iter__()
            return [].__iter__()
        
        self.subIterator = PartitionIterator(self.whole - self.leadPart, self.numParts - 1, self.maxInParts[1:]).__iter__()        
        return self
    
    def next(self):

        while True:
            
            try:
                return [self.leadPart] + self.subIterator.next()
            except StopIteration:
                self.leadPart -= 1
                
                if self.minInParts[0] > self.leadPart:
                    raise StopIteration
                
                self.subIterator = PartitionIterator(self.whole - self.leadPart, self.numParts - 1, self.maxInParts[1:]).__iter__()                                
        
class SliceIterator:
    '''
    >>> l = range(7)
    >>> [s for s in SliceIterator(l, 2)]
    [[0, 1], [2, 3], [4, 5], [6]]
    '''
    
    def __init__(self, iterable, sliceLen=1):
        self._iterable = iterable
        self._sliceLen = sliceLen
        
    def __iter__(self):
        self._iterator = self._iterable.__iter__()
        return self
    
    def next(self):

        slice = []
        try:
            for _ in xrange(self._sliceLen):
                slice.append(self._iterator.next())
        except StopIteration:
            pass
        
        if 0 == len(slice):
            raise StopIteration
        
        return slice        

class PairIterator(object):
    '''
    Iterate through all pairs of items.
    
    >>> l = range(3)
    >>> [p for p in PairIterator(l)]
    [(0, 1), (0, 2), (1, 2)]
    '''
    
    
    def __init__(self, items):
        self.items = items
        
    def __iter__(self):
        self.headIndex = -1
        self._setHead()
        
        return self
        
    def next(self):
        while(True):
            try:
                self.tailIndex += 1
                return (self.head, self.items[self.tailIndex-1])
            except IndexError:
                self._setHead()
            
    def _setHead(self):
        if self.headIndex >= len(self.items)-1:
            raise StopIteration
        
        self.headIndex += 1
        self.tailIndex = self.headIndex + 1
        self.head = self.items[self.headIndex]

class PairSetIterator(object):
    '''
    >>> l = range(3)
    >>> [p for p in PairSetIterator(l)]
    [[(1, 2)], [(2, 0)], [(0, 1)]]
    '''
    
    def __init__(self, items):
        self.items = list(items)
        
    def __iter__(self):
        if 0 == len(self.items) % 2:
            return EvenPairSetIterator(self.items).__iter__()
      
        # The odd case is harder because there is always one item that isn't paired.
        self.unusedIndex = None
        self._setNextUnusedIndex()
      
        return self
  
    def next(self):
        while(True):
            try:
                return self.evenIterator.next()
            except StopIteration:
                self._setNextUnusedIndex()
                
    def _setNextUnusedIndex(self):
        if self.unusedIndex == len(self.items):
            raise StopIteration
        
        if None == self.unusedIndex:
            self.unusedIndex = 0
        else:
            self.unusedIndex += 1
            self.items.append(self.unused)
        
        self.unused = self.items.pop(0)
         
        self.evenIterator = EvenPairSetIterator(self.items).__iter__()
        
class EvenPairSetIterator(object):
    '''
    >>> l = range(4)
    >>> [p for p in EvenPairSetIterator(l)]
    [[(0, 1), (2, 3)], [(0, 2), (3, 1)], [(0, 3), (1, 2)]]
    '''
    
    def __init__(self, items):
        assert 0 == len(items) % 2
        self.items = list(items)
        
    def __iter__(self):
        if 2 >= len(self.items):
            return [[tuple(self.items)]].__iter__()
        
        self.head = self.items.pop(0)
        self.tailIndex = None
        self._setNextTail()
        
        return self
        
    def next(self):
        while(True):
            try:
                return [(self.head, self.tail)] + self.subIterator.next()
            except StopIteration:
                self._setNextTail()
                
    def _setNextTail(self):
        if self.tailIndex == len(self.items):
            raise StopIteration
        
        if None == self.tailIndex:
            self.tailIndex = 0
        else:
            self.tailIndex += 1
            self.items.append(self.tail)
        
        self.tail = self.items.pop(0)
        
        self.subIterator = EvenPairSetIterator(self.items)
        self.subIterator = EvenPairSetIterator(self.items).__iter__()
             
             
class ParallelPairIterator(object):
    '''
    >>> l = [(1,2), (2,3), (1,4), (4,5)]
    >>> [p for p in ParallelPairIterator(l)]
    [[(1, 2), (4, 5)], [(2, 3), (1, 4)], [(2, 3), (4, 5)]]
    '''
    
    def __init__(self, pairs):
        self.pairs = list(pairs)
        
    def __iter__(self):
        if 0 == len(self.pairs):
            return [[]].__iter__()
        
        head = self.pairs[0]
        
        # Select all other pairs that contain an element of head.
        self.myPairs = [pair for pair in self.pairs if (head[0] in pair) or (head[1] in pair)]
        
        self.head = None
        self.usedPairs = []
        self._nextHead()
        
        return self
        
    def next(self):
        while(True):
            try:
                return [self.head] + self.subIterator.next()
            except StopIteration:
                self._nextHead()
                
    def _nextHead(self):
        if 0 == len(self.myPairs):
            raise StopIteration
    
        if None != self.head:
            self.usedPairs.append(self.head)
             
        head = self.myPairs.pop(0)
        self.head = head
        self.pairs.remove(head)
        
        usedCommutingPairs = [0 for pair in self.usedPairs if (pair[0] not in head) and (pair[1] not in head)]
        if len(usedCommutingPairs) != 0:
            # At least one other previously used head pair commutes with the current head.
            # This means that the current head was in "remainingPairs" for that previous head,
            # and therefore there are no more maximal pair combinations with the current head.
            self._nextHead()
        
        remainingPairs = [pair for pair in self.pairs 
                               if (self.head[0] not in pair) and
                                  (self.head[1] not in pair)]
        self.subIterator = ParallelPairIterator(remainingPairs).__iter__()
      
      
def equal_slice_iterators(iterable, num_slices):
    '''
    Chops the sequence into equal size subsequences.
    The last subsequence may contain fewer elements.
    Note: Once the slices have been constructed,
    the original iterator should not be used anywhere else.
    
    >>> map(list, equal_slice_iterators((i for i in range(10)), 3))
    [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9]]
    '''
    
    # In the event that the iterable is an iterator or a generator
    # we will need independent copies for each slice, plus
    # one more to determine the length.
    iter_copies = itertools.tee(iterable, num_slices + 1)
    
    # We want to know the length of the iterable, but
    # we don't want to load the whole thing into memory.
    length = sum(1 for _ in iter_copies[-1])
    
    step = long(math.ceil(length / float(num_slices)))
    slices = [itertools.islice(iter_copies[i], i * step, (i+1) * step)
              for i in xrange(num_slices)]
    return slices
            
if __name__ == "__main__":
    import doctest
    doctest.testmod()