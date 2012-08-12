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
from qfault.util.iteration import ParallelPairIterator
import numpy

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
    
#    # afterA represents the entries of the control that must
#    # not exist until after the overlap is used.
#    afterA = numpy.zeros_like(At)
#    afterA[c] = At[c]
    
    beforeA = numpy.zeros_like(At)
    
    overlap = numpy.multiply(At[c], At[t])
    #print 'overlap=', overlap
    At[t] ^= overlap
    At[c] &= overlap
    
#    # For the after matrix, keep only the entries that are no
#    # longer in the control column of A.
#    afterA[c] ^= At[c]
#    
#    beforeA[c]
    
    #print 'after: '
    #print numpy.transpose(At)
    
    return numpy.transpose(At)#, numpy.transpose(afterA)
    
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
    maxCnots = []
    #maxAfter = numpy.zeros_like(A)
    maxA = copy(A)
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
#        afterA = numpy.zeros_like(A)
        for c,t in pairs:
            #print (c,t)
            net += netMtx[c,t]
            newA = useOverlap(newA, c, t)
            #afterA += after
            
        # Recursively check for more overlaps
        if 0 < len(pairs):
            subnet, subpairs, newA = checkAllPairs(newA, False, 0)
           
            net += subnet
            cnots = subpairs + [pairs]
            if verbose:
                print 'net=', net, 'subnet=', subnet, 'cnots=', cnots
                print 'A=', newA
            if net > maxNet:
                maxNet = net
                maxCnots = cnots
                maxA = newA
                #maxAfter = afterA
    
    if verbose:
        print 'max net=', maxNet
        print 'cnots = ', maxCnots
        print 'A=', maxA
    return maxNet, maxCnots, maxA 
    
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


def constructCircuitFromOverlapCnots(A, cnots, target):
    '''
    Construct an encoding circuit from the given overlap CNOT sequence.
    The matrix A is the redundancy matrix produced by checkAllPairs().  It
    includes all of the entries that must exit prior to the overlap CNOTs.
    '''
    
    
    At = numpy.transpose(A)
    rangeStart = numpy.zeros_like(At)
    
    targetT = numpy.transpose(target)
    rangeEnd = numpy.zeros_like(At)
    # TODO: compute upper bound on the number of rounds
    numpy.putmask(rangeEnd, targetT, 7)
    
    # First, figure out the earliest round in which the overlap CNOT may begin.
    # It can be calculated by looking at weights of the columns associated with each CNOT
    # control in the first round of overlap CNOTs.
    firstOverlapRound = 0
    for c,t in cnots[0]:
        weight = sum(At[c])
        firstOverlapRound = max(firstOverlapRound, weight)
        
    for r in range(len(cnots)):
        round = firstOverlapRound + r
        cnotsForRound = cnots[r]
        
        for c, t in cnotsForRound:
            
            # The overlapping control column bits must appear
            # before this round.
            mins = [min(end, round) for end in rangeEnd[c]]
            numpy.putmask(rangeEnd[c], At[c], mins)
            
            # The remainder of the control column cannot appear until
            # after this round.
            maxs = [max(start, round+1) for start in rangeStart[c]]
            numpy.putmask(rangeStart[c], targetT[c]^At[c], maxs)
            
            numpy.putmask(rangeStart[t], At[c], [round] * len(rangeStart[t]))
            
            # Execute the CNOT
            At[t] ^= At[c]
         
    rangeStart = numpy.transpose(rangeStart)
    rangeEnd = numpy.transpose(rangeEnd)   
    print numpy.transpose(At)
    print rangeStart
    print rangeEnd
    
    rows, cols = numpy.shape(rangeStart)
    ranges = [[0 for _ in range(cols)] for _ in range(rows)]
    for row in range(rows):
        for col in range(cols):
            ranges[row][col] = (rangeStart[row][col], rangeEnd[row][col])
        print ranges[row]
        
    
def constructCircuitFromRanges(target, overlapCnots, rangeStart, rangeEnd):
    pass
        
        
def executeCircuit(stabs, cnots):
    stabsT = numpy.transpose(stabs)
    for c, t in cnots:
        stabsT[t] ^= stabsT[c]
        
    return numpy.transpose(stabsT)
    
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

golayStart = [
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]
        ]

golayCircuit53 = [
(0,13), (4,17), (1,19), (7,14),
(8,13), (5,17), (2,12), (4,14), (10,15), (3,19), (0,20), (7,11), (9,18),
(13,21), (17,22), (8,12), (14,16), (3,15), (1,11), (10,20), (6,18),
(12,17), (19,13), (15,18), (20,22), (6,11), (9,14), (7,21),
(11,12), (13,16), (2,14), (17,15), (18,21), (22,19),
(14,20), (5,11), (10,12), (6,13), (7,15), (1,17), (2,18), (9,19), (8,22),
(1,14), (9,11), (5,13), (6,22),
(3,14), (9,17), (5,18), (0,11), (2,13),
(6,14)
                  ]

if __name__ == "__main__":
#    print overlapMatrix(golay)
#    findOverlapCircuits(golay)

    cnots= [[(2, 10), (6, 11)], [(1, 6), (8, 2), (3, 5), (4, 7), (9, 11)], [(0, 1), (2, 5), (3, 9), (6, 4), (7, 10), (11, 8)]]
    A = [
     [0,0,1,0,0,0,0,0,0,1,0,0],
     [1,0,0,0,0,0,0,0,1,0,0,0],
     [0,1,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,1,0,0,0,1,0,0,0],
     [0,0,0,1,0,0,1,0,0,0,0,0],
     [0,0,0,0,0,0,1,0,0,0,0,0],
     [1,0,0,0,0,0,0,1,0,0,0,0],
     [1,0,0,1,0,0,0,0,0,0,1,0],
     [0,1,1,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,1,0,0,0,0],
     [0,0,0,0,1,0,0,0,0,1,0,0]
     ]


    targetA = numpy.delete(golay, range(len(golay)), axis=1)
    
    #print cnots
    #constructCircuitFromOverlapCnots(A, cnots, targetA)
    
    print len(golayCircuit53)
    result = executeCircuit(golayStart, golayCircuit53)
    print golay - result
