'''
Created on 2011-10-25

@author: adam
'''
from counting.component.base import CountableComponent, Component,\
    ParallelComponent
from counting.key import KeyCopier, KeyMasker, IdentityManipulator,\
    SyndromeKeyGenerator
from counting.location import Locations
from qec.error import zType, xType, Pauli
from qec.qecc import StabilizerState, StabilizerCode
from util import counterUtils, bits
import logging
from counting.countErrors import mapCounts
from util.counterUtils import locrest
from counting.block import Block

logger = logging.getLogger('component')

        
class TransCnot(CountableComponent):
    '''
    Transversal Controlled-NOT.
    '''
    
    ctrlName = 'ctrl'
    targName = 'targ'
    
    def __init__(self, kGood, ctrlCode, targCode, blockorder=[ctrlName, targName]):
        n = ctrlCode.blockLength()
        if n != targCode.blockLength():
            raise Exception('Control ({0}) and target ({1}) blocklengths do not match.'.format(n, targCode.blockLength()))
        
        nickname='transCNOT.'+str(n)
        locs = Locations([counterUtils.loccnot(self.ctrlName, i, self.targName, i) for i in range(n)], nickname)
        
        super(TransCnot, self).__init__(kGood, locs)
        self.blockorder = blockorder
        self.codes = {self.ctrlName: ctrlCode, self.targName: targCode}
        
    def inBlocks(self):
        return tuple(Block(name, self.codes[name]) for name in self.blockorder)
    
    def outBlocks(self):
        outCodes = self.codes
        # Input codes may actually be codewords.  But the output may be
        # entangled, so the output codes are just the underlying code.
        for block, code in outCodes.iteritems():
            try:
                outCodes[block] = code.getCode()
            except AttributeError:
                pass
            
        return tuple(Block(name, outCodes[name]) for name in self.blockorder)
        
    def _logicalChecks(self, codes):
        code = TransCnotCode(codes)
        return code
    
    def propagateOperator(self, op):
        ctrl = self.ctrlName
        targ = self.targName
        newOp = {ctrl: op[ctrl] ^ op[targ][zType],
                 targ: op[targ] ^ op[ctrl][xType]}
        
        return newOp
    
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        
        # TODO: this assumes that the parity checks for both blocks are identical.
        # Need to explicitly check this condition?
        parityChecks = SyndromeKeyGenerator(self.outBlocks()[0].getCode(), None).parityChecks()
        
        # On the control input X errors propagate through to the target
        # block.
        fromCtrlMask = bits.listToBits((0 == check[xType]) for check in parityChecks)
        
        # On the target input Z errors propagate through to the control
        # block.
        fromTargMask = bits.listToBits((0 == check[zType]) for check in parityChecks)
        
        ctrlNum = self.blockorder.index(self.ctrlName)
        targNum = not ctrlNum
        
        ctrlCopier = KeyCopier(subPropagator, ctrlNum, targNum, mask=fromCtrlMask)
        targCopier = KeyCopier(ctrlCopier, targNum, ctrlNum, mask=fromTargMask)
            
        return targCopier
        
class TransCnotCode(StabilizerCode):
    
    def __init__(self, codes):
        self._codes = codes
        
        lengths = [code.n for code in codes.values()]
        if not all(lengths[0] == n for n in lengths):
            raise Exception('Codes are not all of the same length.')
        
        subblockLength = lengths[0]
        
        n = sum(code.n for code in codes.values())
        k = sum(code.k for code in codes.values())
        d = min(code.d for code in codes.values())
        name = ''.join(str(code) for code in codes.values())
        super(TransCnotCode, self).__init__(name, n, k, d)
        
        self.subblockLength = subblockLength

    def stabilizerGenerators(self):
        stabs  = {name: self._codes[name].stabilizerGenerators() 
                  for name in [TransCnot.ctrlName, TransCnot.targName]}
        
        # We assume that the transversal CNOT operation is a valid operation
        # in the code, and therefore the stabilizer generators do not change.
        # Additionally, we assume that the result of the CNOT leaves the
        # two blocks unentangled (as in, e.g., error correction).
        # Thus, we need only extend the operators into the larger space.
        
        I = Pauli.I ** self.subblockLength
        
        block0Stabs = [I + stab for stab in stabs[self._codes.keys()[0]]]
        block1Stabs = [stab + I for stab in stabs[self._codes.keys()[1]]]
        
        return tuple(block0Stabs + block1Stabs)
    
    def normalizerGenerators(self):
        ctrlNorms = self._codes[TransCnot.ctrlName].normalizerGenerators()
        targNorms = self._codes[TransCnot.targName].normalizerGenerators()
        
        return self.PropagateOperators(ctrlNorms, targNorms, self.subblockLength, self._codes.keys())
    
    @staticmethod
    def PropagateOperators(ctrlOps, targOps, blocklength, blockorder):
        newOps = []
        
        I = Pauli.I ** blocklength
        
        # For now, assume that the block ordering is [ctrl, targ].
        for stab in ctrlOps:
            if 0 == stab[zType]:
                # This is an X stabilizer. (X -> XX)
                newOps.append(stab+stab)
            else:
                # This is a Z stabilizer. (Z -> ZI)
                newOps.append(stab+I)
                
        for stab in targOps:
            if 0 == stab[zType]:
                # This is an X stabilizer. (X -> IX)
                newOps.append(I+stab)
            else:
                # This is a Z stabilizer. (Z -> ZZ)
                newOps.append(stab+stab)
                
        # Check our original assumption.  Block ordering is
        # big endian, so blockorder[0] is the MSBs.
        if TransCnot.ctrlName != blockorder[0]:
            # Swap the blocks around.
            shift = blocklength
            mask = (1 << shift) - 1
            newOps = [(stab >> shift) ^ ((stab & mask) << shift) for stab in newOps]
                
        return newOps
    
class TransMeas(CountableComponent):
    
    def __init__(self, kGood, code, basis, blockname='0'):
        n = code.blockLength()
        nickname = 'transMeas' + str(basis) + str(n)
        if Pauli.X == basis:
            loc = counterUtils.locXmeas
            self._basisType = xType
        elif Pauli.Z == basis:
            loc = counterUtils.locZmeas
            self._basisType = zType
        else:
            raise Exception('{0}-basis measurement is not supported'.format(basis))
        
        locs = Locations([loc(blockname, i) for i in range(n)], nickname)
        super(TransMeas, self).__init__(kGood, locs)
        
        self._basis = basis
        self._block = Block(blockname, code)
        
    def inBlocks(self):
        return (self._block,)
        
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        code = self.inBlocks()[0].getCode()
        parityChecks = SyndromeKeyGenerator(code, None).parityChecks()
        
        # Eliminate bits corresponding to checks of the same type as the
        # measurement.  These errors/syndromes cannot be detected.
        mask = bits.listToBits((0 != check[self._basisType]) for check in parityChecks)
        
        blocks = range(len(self.inBlocks()))
        return KeyMasker(subPropagator, mask, blocks=blocks)
        
class TransRest(CountableComponent):
    
    def __init__(self, kGood, code, blockname='0'):
        nickname = 'transRest' + str(code.n)
        locs = Locations([locrest(blockname, i) for i in range(code.n)], nickname)
        
        super(TransRest, self).__init__(kGood, locs)
        
        self._block = Block(blockname, code)
        
    def inBlocks(self):
        return (self._block,)
  
