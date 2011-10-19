'''
Created on 2011-10-17

@author: adam
'''
from qec.error import xType, zType, dualType, PauliError
from qec.qecc import StabilizerCode, Codeword
from util import listutils, bits
from util.cache import memoize

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
    
    def __init__(self, code):
                
        # Ordered list of all parity checks, including logical operators.
        # Some of the stabilizer generators may also be logical operators (e.g,
        # a stabilizer state), so we need to eliminate duplicates.
        parityChecks = code.stabilizerGenerators() + code.normalizerGenerators()
        
        nStabs = len(code.stabilizerGenerators())

        self.code = code
        self._parityChecks = parityChecks
        self.nStabs = nStabs
        
    def parityChecks(self):
        return self._parityChecks
    
    @memoize
    def getKey(self, e):
        key = StabilizerCode.Syndrome(e, self.parityChecks())
        #print 'e=', e, 'parityChecks=', self.parityChecks(), 'pc={0:b} active={1:b} key={2:b}'.format(pc, self.activeBits, key)
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
        normalizers = self.code.normalizerGenerators()
        nNorms = len(normalizers)
        logicalChecks = key & ((1 << nNorms) - 1)
        syndrome = key >> nNorms
        
        e = self.code.syndromeCorrection(syndrome)
        
        decoded = (0,0)
        for i, check in self.code.normalizerGenerators():
            qubit = i/2
            decoded[i % 2] += ((logicalChecks & 1) ^ e.commutesWith(check)) << qubit
            logicalChecks >>= 1
    
        return PauliError(xbits=decoded[0], zbits=decoded[1])
    
    def __repr__(self):
        return 'SyndromeKeyGenerator(' + str(self.code) + ')'
    
class MaskedKeyGenerator(object):
    
    def __init__(self, generator):
        self._generator = generator
        
    def mask(self):
        return 0
    
    def parityChecks(self):
        return self._generator.parityChecks()
    
    def getKey(self, e):
        return self._generator.getKey(e) & self.mask()
    
    def decode(self, key):
        return self._generator.decode(key)
    
class StabilizerStateKeyGenerator(MaskedKeyGenerator):
    
    def __init__(self, state):
        self._state = state
        self.lStabs = set(state.logicalStabilizers())
        code = state.getCode()
        
        # The hash mask eliminates parity checks of logical operators that are in the normalizer
        # of the code, but are not in the normalizer of the state because the corresponding dual
        # operator is now in the stabilizer.
        # TODO: The calculation of the mask here is a bit sticky because it assumes a particular
        # implementation of SyndromeKeyGenerator.getKey().  It would be nice to eliminate this
        # dependency.
        normMasks = [norm in self.lStabs for norm in code.normalizerGenerators()]
        self._mask = bits.listToBits(([1] * len(code.stabilizerGenerators())) + normMasks)
        
        super(StabilizerStateKeyGenerator, self).__init__(SyndromeKeyGenerator(code))
        
    def mask(self):
        return self._mask
    
    def __repr__(self):
        return 'StabilizerStateKeyGenerator(' + str(self._state) + ')'
    
class MultiBlockSyndromeKeyGenerator(object):
    '''
    '''
    
    def __init__(self, blocks):
        self._blocks = blocks
        
        def getGenerator(code):
            if isinstance(code, Codeword):
                return StabilizerStateKeyGenerator(code)
            return SyndromeKeyGenerator(code)
        
        self.generators = {block.name: getGenerator(block.code) for block in blocks}
        
        if 1 == len(blocks):
            self.getKey = self.oneBlockKey
        else:
            self.getKey = self.tupleKey
            
        self.keyLengths = [len(self.generators[block.name].parityChecks()) for block in blocks] 
        
    def parityChecks(self):
        checks = []
        for block in self._blocks:
            checks += self.generators[block.name].parityChecks()
        
        return checks
        
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
        gens = self.generators
        blocks = self._blocks
        keyLengths = self.keyLengths
        
        key = gens[blocks[0].name].getKey(blockErrors[blocks[0].name])
        for i,block in enumerate(blocks[1:]):
            key = (key << keyLengths[i-1]) + gens[block.name].getKey(blockErrors[block.name])
             
        #key = tuple(self.generators[block.name].getKey(blockErrors[block.name]) for block in self._blocks)
        #print 'blockErrors=', blockErrors, 'key=', key
        return key
    
    def __repr__(self):
        gens = tuple(self.generators[block.name] for block in self._blocks)
        return str(gens)



if __name__ == '__main__':
    pass