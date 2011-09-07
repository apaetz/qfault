'''
Created on 2011-08-29

This file contains functions for identifying and exploiting overlap between
stabilizer generators.  Overlap is loosely defined as the amount of similarity
between two or more stabilizer generators.  The idea is that when multiple
stabilizer generators share the same Pauli operators on two qubits, a CNOT gate
between those qubits can copy many Pauli operators at once and thereby reduce
the number of gates required to construct a stabilizer state. 

@author: adam
'''
from copy import copy
from util.iteration import ParallelPairIterator
import numpy
import util.matrix

def overlapMatrix(stabilizers):
    '''
    Returns a matrix O for which the entry O(i,j) contains the
    number of positions in which columns (qubits) i and j overlap.
    '''
    O = numpy.dot(numpy.transpose(stabilizers), stabilizers)
    O[numpy.diag_indices(len(O))] = 0
    return O
    #return numpy.triu(O)
    
def useOverlap(A, c, t):
    '''
    Consumes the overlap for qubits c and t.
    Overlapping entries in the target qubit (t) are
    removed, leaving non-overlapping entries intact.
    Overlapping entries in the control qubit (c) are
    left intact, and non-overlapping entries are removed
    (they cannot exist until *after* the overlap is used).
    '''
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
    '''
    Find overlap-based circuits for the given set
    of X stabilizer generators.
    TODO: gaussian elimination does not do what I want.  For
    now the stabilizers argument must be given in (I|A) format.
    '''
    
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
    '''
    Checks all "maximal" column pairings.  That is
    each column is paired with another column (one column
    may be idle if the total number is odd).  Overlap between
    each pair is exploited.  The remaining stabilizer matrix
    is analyzed by recursively checking all column pairings
    until no useful pairings remain.
    
    Returns the maximum reduction in CNOTs along with the
    schedule of overlap CNOTs.
    '''
    
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
    '''
    Returns a list of pairs (c,t) for which performing
    a CNOT from qubit c to qubit t would result in
    reducing the number of CNOT gates required overall.
    '''
    rowsU, colsU = numpy.triu_indices_from(overlapMtx)
    rowsL, colsL = numpy.tril_indices_from(overlapMtx)
    offDiag = zip(rowsU,colsU) + zip(rowsL,colsL)
    usefulPairs = [pair for pair in offDiag if overlapMtx[pair[0]][pair[1]] > 1]
    return usefulPairs
    
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