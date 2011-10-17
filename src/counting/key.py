'''
Created on 2011-10-17

@author: adam
'''
from qec.error import xType, zType, dualType, PauliError
from qec.qecc import StabilizerCode
from util import listutils, bits

class DefaultErrorKeyGenerator(object):
    
    def getKey(self,e):
        return e
    
    def __repr__(self):
        return 'default_key'
    

class SyndromeKeyGenerator(object):
    '''
    Generates a key based on the error syndrome such that
    for any two errors e1, e2: key(e1) == key(e2) iff e1 and
    e2 are equivalent errors (modulo the stabilizers).
    
    The key is constructed from the error syndrome and some
    of the logical operators. It is given as a tuple of the
    form (syndrome, logical), where syndrome and logical are
    bit strings. A logical operator is used unless either
      1. It is in the stabilizer already.
      2. The corresponding dual logical operator is in the stabilizer.
    In case 1, the check-bit for the logical operator is already
    captured by the syndrome and need not be duplicated.  In case 2,
    errors that otherwise generate the same syndrome, but differ
    with respect to the logical operator are actually equivalent.
    This is because the corresponding dual operator can be applied
    (if necessary, and without affecting the logical state) to 
    ensure that the error always commutes with the logical operator.
    '''
    
    def __init__(self, code, paulis):
        self._paulis = paulis
        
        # For each code, determine which logical operators need to be included in the
        # key.
        logicals = tuple(code.logicalOperator(i,xType) for i in range(code.k)) + \
                   tuple(code.logicalOperator(i,zType) for i in range(code.k))

        
        
        # Ordered list of all parity checks, including logical operators.
        # Some of the stabilizer generators may also be logical operators (e.g,
        # a stabilizer state), so we need to eliminate duplicates.
        stabs = code.stabilizerGenerators()
        parityChecks = listutils.uniqify(logicals + stabs)
        
        stabs = set(stabs)
        nStabs = len(parityChecks) - len(logicals)
        
        # A parity check is active if:
        # 1. It is in the set of stabilizer generators
        # or
        # 2. The corresponding dual operator is not a stabilizer.
        k = code.k
        isActive = [(logicals[(2*i)%k] not in stabs) for i in range(len(logicals))] + \
                   [True] * nStabs
#        activeMeta = []
#        
#        for i in range(len(logicals[xType])):
#            if logicals[zType][i] not in stabs:
#                isActive[i] = 1
#            if not ((logicals[xType][i] in stabs) or (logicals[zType][i] in stabs)):
#                activeChecks += [logicals[dualType(t)][i] for t in paulis]
#                activeMeta += [(i,t) for t in paulis]

        self.code = code
        self.logicals = logicals
        self.parityChecks = parityChecks
        self.activeChecks = isActive
        self.activeBits = bits.listToBits(isActive)
        self.nStabs = nStabs
    
    def getKey(self, e):
        pc = StabilizerCode.Syndrome(e, self.parityChecks)
        key = pc & self.activeBits
        #print 'e=', blockErrors[block.name], 's=', key
        return key
    
#    def getError(self, key):
#        syndrome, logical = key
#        e = self.code.getSyndromeCorrection(syndrome)
#        for l,d in zip(self.logicals, self.dualLogicals):
#            if (logical & 1) == e.commutesWith(l):
#                e *= d
#            logical >>= 1
#        return e
    
    def decode(self, key):
        syndrome = key & ((1 << len(self.nStabs)) - 1)
        logicalChecks = key >> self.nStabs
        
        e = self.code.getSyndromeCorrection(syndrome)
        
        decoded = (0,0)
        k = self.code.k
        for i, check in enumerate(self.logicals):
            decoded[i/k] += (self.activeChecks[i] * ((logicalChecks & 1) ^ e.commutesWith(check))) << i%k
            logicalChecks >>= 1
    
        return PauliError(xbits=decoded[0], zbits=decoded[1])
    
    def __repr__(self):
        paulis = ''.join(str(p) for p in self._paulis)
        return 'syndrome_' + str(self.code) + paulis
    
class SyndromeKeyConverter(object):
    '''
    TODO: this should supply a method for converting from a syndrome key on one block
    to a syndrome key on another block (with a possibly different stabilizer state).
    '''
    
    def __init__(self, fromKeyGenerator, toKeyGenerator):
        pass
    
class MultiBlockSyndromeKeyGenerator(object):
    '''
    '''
    
    def __init__(self, blocks, paulis):
        self._blocks = blocks
        self._paulis = paulis
        
        self.generators = {block.name: SyndromeKeyGenerator(block.code, paulis) for block in blocks}
        
        if 1 == len(blocks):
            self.getKey = self.oneBlockKey
        else:
            self.getKey = self.tupleKey
        
    def oneBlockKey(self, blockErrors):
        block = self._blocks[0]
        e = blockErrors[block.name]
        return self.generators[block.name].getKey(e)

#    def concatenatedKey(self, blockErrors):
#        bits = self._bits
#        paulis = self._paulis
#        # Concatenate the errors on each block into a single bit string.
#        appendSyndrome = lambda s, block: (s << bits[block.name]) + block.code.getSyndrome(blockErrors[block.name], paulis)
#        key = reduce(appendSyndrome, self._blocks, 0)
#        return key

    def tupleKey(self, blockErrors):
        key = tuple(self.generators[block.name].getKey(blockErrors[block.name]) for block in self._blocks)
        #print 'blockErrors=', blockErrors, 'key=', key
        return key
    
    def __repr__(self):
        blocks = self._blocks[0].name.join('.' + block.name for block in self._blocks)
        paulis = ''.join(str(p) for p in self._paulis)
        return 'syndrome_' + blocks + paulis



if __name__ == '__main__':
    pass