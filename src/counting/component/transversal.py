'''
Created on 2011-10-25

@author: adam
'''
from counting.component.base import CountableComponent, Component
from counting.key import keyExtender, keyCopier
from counting.location import Locations
from qec.error import zType, xType, Pauli
from qec.qecc import StabilizerState, StabilizerCode
from util import counterUtils, bits
import logging

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
        
        if isinstance(ctrlCode, StabilizerState): ctrlCode = ctrlCode.getCode()
        if isinstance(targCode, StabilizerState): targCode = targCode.getCode()
        
        codes = {self.ctrlName: ctrlCode, self.targName: targCode}
        
        super(TransCnot, self).__init__(locs, blockorder, codes, kGood, nickname)
        self.blockorder = blockorder
        
    def _logicalChecks(self, codes):
        code = TransCnotCode(codes)
        return code
    
    def propagateOperator(self, op):
        ctrl = self.ctrlName
        targ = self.targName
        newOp = {ctrl: op[ctrl] ^ op[targ][zType],
                 targ: op[targ] ^ op[ctrl][xType]}
        
        return newOp
    
    def keyPropagator(self, keyMeta, blockname):
        parityChecks = keyMeta.parityChecks()
        
        blocknum = self.blockorder.index(blockname)
        
        # Extend the single input block to both blocks
        before = 1 * (1 == blocknum)
        after = 1 * (0 == blocknum)
        extender, keyMeta = keyExtender(keyMeta, blocksBefore=before, blocksAfter=after)
        
        if self.ctrlName == blockname:
            # We're on the control input.  X errors propagate through to the target
            # block.
            mask = bits.listToBits((0 == check[xType]) for check in parityChecks)
        else:
            # We're on the target input.  Z errors propagate through to the control
            # block.
            mask = bits.listToBits((0 == check[zType]) for check in parityChecks)
        
        copier = keyCopier(keyMeta, blocknum, not blocknum, mask=mask)
    
        def propagator(key):
            return copier(extender(key))
        
        return propagator, keyMeta
        
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
    
class CnotConvolver(Component):
    
    ctrlName = 'ctrl'
    targName = 'targ'
    cnotName = 'cnot'
    
    def __init__(self, kGood, kGoodCnot, ctrlInput, targInput, ctrlName=ctrlName, targName=targName):
        
        # Construct a transversal CNOT component from the two input codes.
        ctrlCode = ctrlInput.outBlocks()[0].getCode()
        targCode = targInput.outBlocks()[0].getCode()
        cnot = TransCnot(kGoodCnot, ctrlCode, targCode)
        
        super(CnotConvolver, self).__init__(kGood, subcomponents={ctrlName: ctrlInput, 
                                                                  targName: targInput,
                                                                  self.cnotName: cnot})
        self.ctrlName = ctrlName
        self.targName = targName
        
    def outBlocks(self):
        self.subcomponents()[self.cnotName].outBlocks()
    
    def _convolve(self, results, pauli):
        
        cnot = self.subcomponents()[self.cnotName]
        ctrlResult = results[self.ctrlName]
        targResult = results[self.targName]
        cnotResult = results[self.cnotName]
        
        # First, propagate the input results through the CNOT    
        ctrlResult.counts, ctrlResult.keyMeta = cnot.propagateCounts(ctrlResult.counts,
                                                             ctrlResult.keyMeta,
                                                             cnot.ctrlName)
        targResult.counts, targResult.keyMeta = cnot.propagateCounts(targResult.counts,
                                                             targResult.keyMeta,
                                                             cnot.targName)
        
        ctrlResult.blocks = cnotResult.blocks
        targResult.blocks = cnotResult.blocks
                    
        # Now convolve.
        return super(CnotConvolver, self)._convolve(results, pauli)            

    
class TransMeas(CountableComponent):
    
    def __init__(self, kGood, code, basis, blockname='0'):
        n = code.blockLength()
        nickname = 'transMeas' + str(basis) + str(n)
        if Pauli.X == basis:
            loc = counterUtils.locXmeas
        elif Pauli.Z == basis:
            loc = counterUtils.locZmeas
        else:
            raise Exception('{0}-basis measurement is not supported'.format(basis))
        
        locs = Locations([loc(blockname, i) for i in range(n)], nickname)
        super(TransMeas, self).__init__(locs, [blockname], {blockname: code}, kGood, nickname)
