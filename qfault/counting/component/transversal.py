'''
Created on 2011-10-25

Basic transversal components such as: CNOT, rest, measurement.

@author: adam
'''
from copy import copy
from qfault.circuit import location
from qfault.circuit.block import Block
from qfault.circuit.location import Locations
from qfault.counting.component.base import CountableComponent, Component, \
    ParallelComponent
from qfault.counting.count_locations import map_counts
from qfault.counting.key import KeyCopier, KeyMasker, IdentityManipulator, \
    SyndromeKeyGenerator
from qfault.qec.error import zType, xType, Pauli
from qfault.qec.qecc import StabilizerState, StabilizerCode
from qfault.util import counterUtils, bits
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
        locs = Locations([location.cnot(self.ctrlName, i, self.targName, i) for i in range(n)], nickname)
        
        super(TransCnot, self).__init__(kGood, locs)
        self.blockorder = blockorder
        self.codes = {self.ctrlName: ctrlCode, self.targName: targCode}
        
    def inBlocks(self):
        return tuple(Block(name, self.codes[name]) for name in self.blockorder)
    
    def outBlocks(self):
        outCodes = copy(self.codes)
        # Input codes may actually be codewords.  But the output may be
        # entangled, so the output codes are just the underlying code.
        for block, code in outCodes.iteritems():
            try:
                outCodes[block] = code.get_code()
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
        code = self.outBlocks()[0].get_code()
        parityChecks = SyndromeKeyGenerator(code).parityChecks()
        identity = Pauli.I ** code.blockLength()
        
        # On the control input X errors propagate through to the target
        # block.
        fromCtrlMask = bits.listToBits((identity == check.partial(Pauli.X)) 
                                       for check in parityChecks)
        
        # On the target input Z errors propagate through to the control
        # block.
        fromTargMask = bits.listToBits((identity == check.partial(Pauli.Z)) 
                                       for check in parityChecks)
        
        ctrlNum = self.blockorder.index(self.ctrlName)
        targNum = ctrlNum ^ 1
        
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
        locs = Locations([location.meas(basis, blockname, i) for i in range(n)], nickname)
        super(TransMeas, self).__init__(kGood, locs)
        
        self._basis = basis
        self._block = Block(blockname, code)
        
    def inBlocks(self):
        return (self._block,)
        
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        code = self.inBlocks()[0].get_code()
        parityChecks = SyndromeKeyGenerator(code).parityChecks()
        identity = Pauli.I ** code.blockLength()
        
        # Eliminate bits corresponding to checks of the same type as the
        # measurement.  These errors/syndromes cannot be detected.
        mask = bits.listToBits(identity != check.partial(self._basis) 
                               for check in parityChecks)
        
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
  
