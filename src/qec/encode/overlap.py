'''
Created on 2011-08-29

@author: adam
'''
from adt.BitString import BitString
from copy import copy
from util.iteration import PairSetIterator, ParallelPairIterator
import numpy
import util.matrix

def overlapMatrix(stabilizers):
    O = numpy.dot(numpy.transpose(stabilizers), stabilizers)
    O[numpy.diag_indices(len(O))] = 0
    return O
    #return numpy.triu(O)
    
def useOverlap(A, c, t):
    #print 'before (c=', c, ', t=', t, '):'
    #print A
    At = numpy.transpose(A)
    overlap = numpy.multiply(At[c], At[t])
    #print 'overlap=', overlap
    At[t] ^= overlap
    At[c] &= overlap
    
    #print 'after: '
    #print numpy.transpose(At)
    
    return numpy.transpose(At)
    
def findOverlapCircuits(stabilizers):
    if True != util.matrix.gauss_jordan(stabilizers):
        print stabilizers
        raise Exception
    
    stabilizers = numpy.array(stabilizers)
    print stabilizers
    
    A = numpy.delete(stabilizers, range(len(stabilizers)), axis=1)
    
    print 'A='
    print A
    
    net = checkAllPairs(A, verbose=True)
    print net
        
def checkAllPairs(A, verbose=False, level=0):
    
    if level > 4:
        raise Exception
    
    overlapMtx = overlapMatrix(A)
    if verbose:
        print overlapMtx
    netMtx = overlapMtx - 1
    usefulPairs = getUsefulPairs(overlapMtx)
    if verbose:
        print 'useful pairs:', usefulPairs
    
    maxNet = 0
    maxCnots = ()
    if verbose:
        sizes = [0, 0, 0, 0, 0, 0, 0, 0]
        for pairs in ParallelPairIterator(usefulPairs):
            sizes[len(pairs)] += 1
            
        print 'total pairs to try:', sum(sizes), sizes
        
    for pairs in ParallelPairIterator(usefulPairs):
        #if verbose:
        #    print pairs
        net = 0
        newA = copy(A)
        for c,t in pairs:
            #print (c,t)
            net += netMtx[c,t]
            newA = useOverlap(newA, c, t)
            
        # Recursively check for more overlaps
        if 0 < len(pairs):
            subnet, subpairs = checkAllPairs(newA, False, 0)
           
            net += subnet
            cnots = (subpairs, pairs)
            if verbose:
                print 'net=', net, 'subnet=', subnet, 'cnots=', cnots
            if net > maxNet:
                maxNet = net
                maxCnots = cnots
    
    if verbose:
        print 'max net=', maxNet
        print 'cnots = ', maxCnots
    return maxNet, maxCnots 
    
def getUsefulPairs(overlapMtx):
    rowsU, colsU = numpy.triu_indices_from(overlapMtx)
    rowsL, colsL = numpy.tril_indices_from(overlapMtx)
    offDiag = zip(rowsU,colsU) + zip(rowsL,colsL)
    usefulPairs = [pair for pair in offDiag if overlapMtx[pair[0]][pair[1]] > 1]
    return usefulPairs

def pairsNotContaining(pairs, value):
    return [pair for pair in pairs if value not in pair]
    
# Steane [[7,1,3]]] X stabilizers
#steane7 = [
#           [0,0,0,1,1,1,1],
#           [1,0,1,0,1,0,1],
#           [0,1,1,0,0,1,1]
#           ]

steane7 = [
           [0,0,1,0,1,1,1],
           [1,0,0,1,1,0,1],
           [0,1,0,1,0,1,1]
           ]

golay = [
        [1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,1,1,1],
        [0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,1,0,0,0],
        [0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,1,0,0],
        [0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,1,0],
        [0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,1],
        [0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,1,0,1,1,1,0,0,1],
        [0,0,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,1,0,0,1,1],
        [0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,1,0,0,0,1,1,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,1,0,0,0,1,1],
        [0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,1,1,1,1,0],
        [0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,1,1,1,1]
        ]

if __name__ == "__main__":
    print overlapMatrix(golay)
    
    findOverlapCircuits(golay)