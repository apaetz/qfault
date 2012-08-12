'''
Created on 2011-10-17

Classes and functions for creating and manipulating syndrome keys.
The idea of this module is to abstract error syndrome processing away from the Component internals.
The outcome has been mixed because many components manipulate error syndromes in complicated ways.
More work is needed to come up with a satisfactory solution.

@author: adam
'''
from collections import OrderedDict
from qfault.counting.convolve import convolve_dict
from qfault.qec.error import xType, zType, dualType, PauliError, Pauli
from qfault.qec.qecc import StabilizerCode, Codeword
from qfault.util import listutils, bits
from qfault.util.cache import memoize
import logging

logger = logging.getLogger('counting.key')

rejectKey = 0

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
    
    @staticmethod
    def ParityChecks(stabilizerCode):
        try:
            # Get the underlying code, if dealing with a stabilizer state.
            stabilizerCode = stabilizerCode.getCode()
        except AttributeError:
            pass
        
        return stabilizerCode.stabilizerGenerators() + stabilizerCode.normalizerGenerators()
    
    def __init__(self, code, blockname):
                
        # Ordered list of all parity checks, including logical operators.
        # Some of the stabilizer generators may also be logical operators (e.g,
        # a stabilizer state), so we need to eliminate duplicates.
        parityChecks = self.ParityChecks(code)
        
        nStabs = len(code.stabilizerGenerators())

        self.blockname = blockname
        self.code = code
        self._parityChecks = parityChecks
        self.nStabs = nStabs
        
    def parityChecks(self):
        return self._parityChecks
    
    @memoize
    def getKey(self, e):
#        # e is given in little endian bit order (bit 0 = qubit 0), but
#        # Qecc objects require big endian.
#        X = bits.endianSwap(e[xType], self.code.n)
#        Z = bits.endianSwap(e[zType], self.code.n)
#        e = PauliError(xbits=X, zbits=Z)
        
        key = StabilizerCode.Syndrome(e, self.parityChecks())
        #print 'e=', e, 'parityChecks=', self.parityChecks(), 'key={0:b}'.format(key)
        return key
    
#    def keyMeta(self):
#        return SyndromeKeyMeta(self.parityChecks(), 1)
    
#    def getError(self, key):
#        syndrome, logical = key
#        e = self.code.getSyndromeCorrection(syndrome)
#        for l,d in zip(self.logicals, self.dualLogicals):
#            if (logical & 1) == e.commutesWith(l):
#                e *= d
#            logical >>= 1
#        return e
    
#    def decode(self, key):
#        normalizers = self.code.normalizerGenerators()
#        nNorms = len(normalizers)
#        logicalChecks = key & ((1 << nNorms) - 1)
#        syndrome = key >> nNorms
#        
#        e = self.code.syndromeCorrection(syndrome)
#        
#        decoded = (0,0)
#        for i, check in self.code.normalizerGenerators():
#            qubit = i/2
#            decoded[i % 2] += ((logicalChecks & 1) ^ e.commutesWith(check)) << qubit
#            logicalChecks >>= 1
#    
#        return PauliError(xbits=decoded[0], zbits=decoded[1])
    
    def __repr__(self):
        return 'SyndromeKeyGenerator(' + str(self.code) + ')'
    
    
def IntegerKey(value, nblocks=1):
    return tuple([value]*nblocks)
    
identity = lambda key: key
    
class KeyManipulator(object):
    
    def __init__(self, manipulator=identity):
        self._manipulator = manipulator
                
    def __call__(self, key):
        return self._manipulate(self._manipulator(key))
        
    def _manipulate(self, key):
        raise NotImplementedError
    
class IdentityManipulator(KeyManipulator):
    
    def _manipulate(self,key):
        return key
    
class KeyExtender(KeyManipulator):
    
    def __init__(self, manipulator, numBlocks, insertIndex):
        super(KeyExtender, self).__init__(manipulator)
        self._extension = tuple([0] * numBlocks)
        self._index = insertIndex
        
    def _manipulate(self, key):
        return key[:self._index] + self._extension + key[self._index:]
    
class KeyRemover(KeyManipulator):
    
    def __init__(self, manipulator, remove_indices):
        super(KeyRemover, self).__init__(manipulator)
        self._remove_indices = set(remove_indices)
        
    def _manipulate(self, key):
        key = listutils.remove_subsequence(key, self._remove_indices)
        return tuple(key)
    
class KeySplitter(KeyManipulator):
    
    def __init__(self, manipulator, splits):
        super(KeySplitter, self).__init__(manipulator)
        self._splits = splits
    

    def _manipulate(self, key):
        keys = [0] * (len(self._splits) + 1)
        lastSplit = 0
        
        for i,split in enumerate(self._splits):
            keys[i] = tuple(key[lastSplit:split])
            lastSplit = split
        keys[-1] = tuple(key[lastSplit:])
        
        return keys
    
class KeyPermuter(KeyManipulator):
    
    def __init__(self, manipulator, permutation):
        super(KeyPermuter, self).__init__(manipulator)
        self.permutation = permutation
        
    def _manipulate(self, key):
        return tuple(listutils.permute(key[:len(self.permutation)], self.permutation)) + key[len(self.permutation):]
    
class KeyConcatenator(KeyManipulator):
    
    def __init__(self, *manipulators):
        if 1 == len(manipulators):
            manipulator = manipulators[0]
        else:
            manipulator = MultiManipulatorAdapter(manipulators)
            
            
        super(KeyConcatenator, self).__init__(manipulator)

    def _manipulate(self, keys):
        return sum(keys, tuple())
    
class KeyMerger(KeyManipulator):
    '''
    Merges a key with multiple blocks into a key with a single block.
    '''
    def __init__(self, manipulator, keyLengths):
        super(KeyMerger, self).__init__(manipulator)
        self.keyLengths = keyLengths
        
    def _manipulate(self, key):
        merged = bits.concatenate(key[:len(self.keyLengths)], self.keyLengths, reverse=True)
        return (merged,) + key[len(self.keyLengths):]

class MultiManipulatorAdapter(object):
    
    def __init__(self, manipulators):
        self._manipulators = manipulators

    def __call__(self, keys):
        return tuple(self._manipulators[i](key) for i,key in enumerate(keys))
    
    
def keyForBlock(key, block):
    return (key[block],)

class KeyCopier(KeyManipulator):
    
    def __init__(self, manipulator, fromBlock, toBlock, mask=None):
        super(KeyCopier, self).__init__(manipulator)
        if None == mask:
            self._manipulate = self._manipulateNoMask
        else:
            self._mask = mask
            self._manipulate = self._manipulateMask
            
        self._toBlock = toBlock
        self._fromBlock = fromBlock
            
    def _manipulateMask(self, key):
        newKey = list(key)
        try:
            newKey[self._toBlock] ^= key[self._fromBlock] & self._mask
        except:
            pass
        return tuple(newKey)
    
    def _manipulateNoMask(self, key):
        newKey = list(key)
        newKey[self._toBlock] ^= key[self._fromBlock]
        return tuple(newKey)
    
class KeyMasker(KeyManipulator):
    
    def __init__(self, manipulator, mask, blocks=None):
        super(KeyMasker, self).__init__(manipulator)
        
        if None == blocks:
            self._manipulate = self._manipulateAll
        else:
            self.blocks = blocks
            self._manipulate = self._manipulateSelected
            
        self.mask = mask
        
    
    
    def _manipulateSelected(self, key):
        mask = self.mask
        newKey = list(key)
        # Mask only the selected blocks.
        for block in self.blocks:
            newKey[block] &= mask
        return tuple(newKey)
    
    def _manipulateAll(self, key):
        mask = self.mask
        newKey = tuple(k & mask for k in key)
        return newKey
    
class SyndromeKeyFilter(KeyManipulator):
        
    def __init__(self, code, manipulator):
        super(SyndromeKeyFilter, self).__init__(manipulator)
        self._nNorms = len(code.normalizerGenerators())
    
    def _manipulate(self, key):
        return ((key[0] >> self._nNorms) << self._nNorms,) + key[1:]
        
def convolveKeyCounts(counts1, counts2, key1Lengths, key2Lengths):
    
    # Return counts according to the larger of the two keys.
    if len(key2Lengths) > len(key1Lengths):
        keyLengths = key2Lengths
        k = len(key1Lengths)
    else:
        keyLengths = key1Lengths
        k = len(key2Lengths)
    
    if key1Lengths[:k] != key2Lengths[:k]:
        raise Exception('Incompatible key lengths {0}, {1}'.format(key1Lengths, key2Lengths))
    
    try:
        counts1 = {bits.concatenate(key, key1Lengths, reverse=True): count for key,count in counts1.iteritems()}
        counts2 = {bits.concatenate(key, key2Lengths, reverse=True): count for key,count in counts2.iteritems()}
    except Exception:
        raise
    
    counts = convolve_dict(counts1, counts2)
    
    return {bits.split(keybits, keyLengths, reverse=True): count for keybits,count in counts.iteritems()}
    
    
    

    
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
    
#    def keyMeta(self):
#        return self._generator.keyMeta()
    
class StabilizerStateKeyGenerator(MaskedKeyGenerator):
    
    def __init__(self, state, blockname):
        self._state = state
        #self.lStabs = set(state.logicalStabilizers())
        #code = state.getCode()
        
        skg = SyndromeKeyGenerator(state.getCode(), blockname)

        # The check mask eliminates parity checks of logical operators that are in the normalizer
        # of the code, but are not in the normalizer of the state because the corresponding dual
        # operator is now in the stabilizer.
        maskedChecks = set(state.stabilizerGenerators() + state.normalizerGenerators())
        masks = [check in maskedChecks for check in skg.parityChecks()]
        self._mask = bits.listToBits(masks)
        
        super(StabilizerStateKeyGenerator, self).__init__(skg)
        
    def mask(self):
        return self._mask
    
    def __repr__(self):
        return 'StabilizerStateKeyGenerator(' + str(self._state) + ')'
    
class MultiBlockSyndromeKeyGenerator(object):
    '''
    '''
    
    def __init__(self, blocks):
        self._blocks = blocks
        self._blocknames = [block.name for block in blocks]
        
        def getGenerator(code, blockname):
            if isinstance(code, Codeword):
                return StabilizerStateKeyGenerator(code, blockname)
            return SyndromeKeyGenerator(code, blockname)
        
        self.generators = {block.name: getGenerator(block.code, block.name) for block in blocks}
        
        pcSet = set([g.parityChecks() for g in self.generators.values()])
        if 1 < len(pcSet):
            raise Exception('Parity check mismatch. {0}'.format(pcSet))
        
        self._parityChecks = pcSet.pop()
        
#        # oneBlockKey is slightly more efficient.
#        if 1 == len(blocks):
#            self.getKey = self.oneBlockKey
#        else:
#            self.getKey = self.multiBlockKey
        
    def parityChecks(self):
        return self._parityChecks
    
#    def keyMeta(self):
#        return SyndromeKeyMeta(self._parityChecks, len(self._blocknames))
        
#    def oneBlockKey(self, blockErrors):
#        block = self._blocks[0]
#        e = blockErrors[block.name]
#        return self.generators[block.name].getKey(e)

    def getKey(self, blockErrors):
        gens = self.generators
        blocknames = self._blocknames
        
        key = tuple(gens[name].getKey(blockErrors.get(name, Pauli.I)) for name in blocknames)
        
        #print 'blockErrors=', blockErrors, 'key=', key
        return key
    
    def __repr__(self):
        gens = tuple(self.generators[block.name] for block in self._blocks)
        return str(gens)
    
    
class SyndromeKeyDecoder(object):
    
    def __init__(self, code):
        self._code = code
        self._normalizers = code.normalizerGenerators()
        
    def syndrome(self, key):
        return key >> len(self._normalizers)
        
    def decode(self, key):
#        logicalError = Pauli.I
        nNorms = len(self._normalizers)
        logicalChecks = key & ((1 << nNorms) - 1)
        syndrome = key >> nNorms
        
        e = self._code.syndromeCorrection(syndrome)
        
        commutations = [not e.commutesWith(check) for check in self._normalizers]
        commutations = bits.listToBits(commutations)
        decoded = logicalChecks ^ commutations
        
        logger.debug('correction=%s, normalizers=%s', e, self._normalizers)
        logger.debug('key=%s, decoded key=%s', key, decoded)
        return decoded
        
        
    def asPauli(self, key):
        nNorms = len(self._normalizers)
        logicalChecks = key & ((1 << nNorms) - 1)
        logicalChecks = bits.bitsToList(logicalChecks, nNorms)
        normalizerLookup = {norm: logicalChecks[i] for i,norm in enumerate(self._normalizers)}
        
        syndrome = key >> nNorms
        e = self._code.syndromeCorrection(syndrome)
        
        for logical in self._code.logicalOperators():
            Xl = logical[xType]
            Zl = logical[zType]
            if normalizerLookup[Xl]:
                # The error anti-commutes with the logical X operator.
                e *= Zl
            
            if normalizerLookup[Zl]:
                # The error anti-commutes with the logical Z operator.
                e *= Xl
                
        return e
            

if __name__ == '__main__':
    pass