'''
Created on 2011-10-17

@author: adam
'''
from collections import OrderedDict
from counting.convolve import convolveDict
from qec.error import xType, zType, dualType, PauliError, Pauli
from qec.qecc import StabilizerCode, Codeword
from util import listutils, bits
from util.cache import memoize
import logging
import operator

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
    
#def IntegerKeyMeta(nblocks=1):
#    return SyndromeKeyMeta([Pauli.I], nblocks=nblocks)
    
#class SyndromeKeyMeta(object):
#    
#    def __init__(self, parityChecks, nblocks=1):
#        self._parityChecks = tuple(parityChecks)
#        self.nblocks = nblocks
#
#    def length(self):
#        return len(self.parityChecks()) * self.nblocks
#        
#    def parityChecks(self):
#        return self._parityChecks
#    
#    def meta(self):
#        return self
#    
#    def __call__(self, key):
#        return key
#    
##    def blocknames(self):
##        return self._blocks.keys()
##    
##    def blockIndex(self, blockname):
##        return self._blocks[blockname]
##    
##    def blockRange(self, blockname):
##        blocklen = len(self.parityChecks())
##        start = self.blockIndex(blockname) * blocklen
##        end = start + blocklen
##        return range(start, end+1)
##    
##    def syndromeOf(self, key):
##        return key >> self._normalizerBits
##    
##    def syndromeMask(self):
##        return  ((1 << self.length()) - 1) ^ ((1 << self._normalizerBits) - 1)
#    
#    def __eq__(self, other):
#        try:
#            return (self.parityChecks() == other.parityChecks()) and (self.nblocks == other.nblocks) 
#        except:
#            return False
#        
#    def __hash__(self):
#        return self.parityChecks().__hash__() + self.nblocks
#    
#    def __str__(self):
#        return str(self.nblocks) + '*' + str(self.parityChecks())
#    
#    def __repr__(self):
#        return 'SyndromeKeyMeta(' + str(self.parityChecks()) + ',' + str(self.nblocks) + ')'
    
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
    
    def __init__(self, manipulator, start, end):
        super(KeyRemover, self).__init__(manipulator)
        self.start = start
        self.endplus1 = end + 1
        
    def _manipulate(self, key):
        return key[:self.start] + key[self.endplus1:]
    
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
    
class KeyConcatenator(KeyManipulator):
    
    def __init__(self, *manipulators):
        if 1 == len(manipulators):
            manipulator = manipulators[0]
        else:
            manipulator = MultiManipulatorAdapter(manipulators)
            
            
        super(KeyConcatenator, self).__init__(manipulator)

    def _manipulate(self, keys):
        return sum(keys, tuple())

class MultiManipulatorAdapter(object):
    
    def __init__(self, manipulators):
        self._manipulators = manipulators

    def __call__(self, keys):
        return tuple(self._manipulators[i](key) for i,key in enumerate(keys))
        
#def extendKeys(keys, keyMeta, blocksBefore=0, blocksAfter=0):
#    before = tuple([0] * blocksBefore)
#    after = tuple([0] * blocksAfter)
#    keymap = {key: before + key + after for key in keys}
#    keyMeta = SyndromeKeyMeta(keyMeta.parityChecks(), keyMeta.nblocks + blocksBefore + blocksAfter)
#    
#    return keymap, keyMeta
#
#def keyExtender(keyMeta, blocksBefore=0, blocksAfter=0):
#    before = tuple([0] * blocksBefore)
#    after = tuple([0] * blocksAfter)
#    def extendKey(key):
#        return before + key + after
#    
#    keyMeta = SyndromeKeyMeta(keyMeta.parityChecks(), keyMeta.nblocks + blocksBefore + blocksAfter)
#    
#    return extendKey, keyMeta
#
#def keySplitter(keyMeta, splitBlock):
#    
#    meta1 = SyndromeKeyMeta(keyMeta.parityChecks(), splitBlock)
#    meta2 = SyndromeKeyMeta(keyMeta.parityChecks(), keyMeta.nblocks - splitBlock)
#    
#    def split(key):
#        return key[:splitBlock], key[splitBlock:]
#
#    return split, meta1, meta2
#
#def keyConcatenator(keyMeta1, keyMeta2):
#    
#    if keyMeta1.parityChecks() != keyMeta2.parityChecks():
#        raise Exception('Incompatible keys {0}, {1}'.format(keyMeta1, keyMeta2))
#    
#    meta = SyndromeKeyMeta(keyMeta1.parityChecks(), keyMeta1.nblocks + keyMeta2.nblocks)
#    
#    def cat(key1, key2):
#        return key1 + key2
#    
#    return cat, meta

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
        newKey[self._toBlock] ^= key[self._fromBlock] & self._mask
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
    
#def keyCopier(keyMeta, fromBlock, toBlock, mask=None):
#    '''
#    Copy keys (with the given key metatadata) from one block to another block.  The optional
#    mask specifies which bits of fromBlock are copied.  By default, all bits are copied.
#    '''
#            
#    if None == mask:
#        blocklen = len(keyMeta.parityChecks())
#        # Select all of the bits of fromBlock
#        mask = (1 << blocklen) - 1
#    
#    def newKey(key):
#        newKey = list(key)
#        newKey[toBlock] ^= key[fromBlock] & mask
#        return tuple(newKey)
#    
#    return newKey

#def copyKeys(keys, fromBlock, toBlock, mask=None):
#    '''
#    Copy keys from one block onto another block.  The optional
#    mask specifies which bits of fromBlock are copied.  By default, all bits are copied.
#    '''
#            
#    if None == mask:
#        blocklen = len(keyMeta.parityChecks())
#        # Select all of the bits of fromBlock
#        mask = (1 << blocklen) - 1
#    
#    def newKey(key):
#        newKey = list(key)
#        newKey[toBlock] ^= key[fromBlock] & mask
#        return tuple(newKey)
#    
#    return {key: newKey(key) for key in keys}
        
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
    
    counts = convolveDict(counts1, counts2)
    
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
        
    def decode(self, key):
#        logicalError = Pauli.I
        nNorms = len(self._normalizers)
        logicalChecks = key & ((1 << nNorms) - 1)
        syndrome = key >> nNorms
        
        e = self._code.syndromeCorrection(syndrome)
        
#        # TODO: this algorithm is neither robust nor readable.
#        decoded = [0,0]
#        for i, check in enumerate(reversed(self._normalizers)):
#            qubit = i/2
#            decoded[i % 2] += ((logicalChecks & 1) ^ (not e.commutesWith(check))) << qubit
#            logicalChecks >>= 1

#        logicalError += PauliError(xbits=decoded[0], zbits=decoded[1])
        
        decoded = 0
        commutations = [not e.commutesWith(check) for check in self._normalizers]
        commutations = bits.listToBits(commutations)
        decoded = logicalChecks ^ commutations
        
        logger.debug('correction=%s, normalizers=%s', e, self._normalizers)
        logger.debug('key=%s, decoded key=%s', key, decoded)
        return decoded
        

if __name__ == '__main__':
    pass