'''
Created on 2011-10-25

Basic transversal components such as: CNOT, rest, measurement.

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
        
#    def _logicalChecks(self, codes):
#        code = TransCnotCode(codes)
#        return code
    
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
        
           
class TransMeas(CountableComponent):
    '''
    Transversal measurement in either the 'X' or 'Z' basis.
    '''
    
    def __init__(self, kGood, code, basis, blockname=''):
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
    '''
    Transversal rest.
    '''
    
    def __init__(self, kGood, code, blockname=''):
        nickname = 'transRest' + str(code.n)
        locs = Locations([locrest(blockname, i) for i in range(code.n)], nickname)
        
        super(TransRest, self).__init__(kGood, locs)
        
        self._block = Block(blockname, code)
        
    def inBlocks(self):
        return (self._block,)
  
