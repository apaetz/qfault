'''
Created on 2012-04-02

@author: adam
'''

from qfault.counting.block import Block
from qfault.counting.component import transversal
from qfault.counting.component.base import ParallelComponent, SequentialComponent, Prep, Filter
from qfault.counting.component.bell import BellPair, BellMeas
from qfault.counting.component.block import BlockDiscard,\
    BlockCombine
from qfault.counting.component.teleport import Teleport, TeleportED
from qfault.counting.key import IdentityManipulator, SyndromeKeyGenerator,\
    KeyManipulator
from qfault.qec import ed422
from qfault.qec.error import Pauli, xType
from qfault.util import bits, listutils
import logging

logger = logging.getLogger('scheme.fibonacci')

class SubblockTeleport(SequentialComponent):
    
    def __init__(self, kGood, bell_pair, block_to_teleport=0, teleport_output_block=2, postselect=False):
        code = bell_pair.outBlocks()[0].getCode()
        
        cnot = transversal.TransCnot(kGood, code, code)
        discard = BlockDiscard(cnot.outBlocks(), [block_to_teleport ^ 1])
        bm = BellMeas(kGood, code, kGood, kGood, kGood)
        if postselect:
            if 2 != teleport_output_block:
                raise RuntimeError('teleportation with postselection supports output block 2 only')
            teleport = TeleportED(kGood, bell_pair, bm)
        else:
            teleport = Teleport(kGood, bell_pair, bm, output_block=teleport_output_block)
        
        super(SubblockTeleport, self).__init__(kGood, subcomponents=[bell_pair, cnot, discard, teleport])
        
    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
        r = super(SubblockTeleport, self).count(noiseModels, pauli, inputResult, kMax)
        return r
        
class BlockTeleport(SequentialComponent):
    
    def __init__(self, kGood, bell_pair, meas_basis, block_to_meas):
        code = bell_pair.outBlocks()[0].getCode()
        cnot = transversal.TransCnot({p: 1 for p in kGood.keys()}, code, code)
        
        bpc = SequentialComponent(kGood, subcomponents=[bell_pair, cnot])
        
        discard1 = BlockDiscard(cnot.outBlocks(), [0])
        if 1 == block_to_meas:
            discard2 = BlockDiscard(cnot.outBlocks(), [1])
        else:
            discard2 = BlockDiscard(cnot.outBlocks(), [0])

        bp_cnot1 = SequentialComponent(kGood, [bpc, discard1])
        bp_cnot2 = SequentialComponent(kGood, [bpc, discard2])
        bp_cnot_parallel = ParallelComponent(kGood, bp_cnot1, bp_cnot2)
            
        discard3 = BlockDiscard(cnot.outBlocks(), [block_to_meas - 1])
        
        if xType == meas_basis:
            m = Pauli.X
        else:
            m = Pauli.Z
        meas = transversal.TransMeas(kGood, code, m)
        
        super(BlockTeleport, self).__init__(kGood, subcomponents=[bp_cnot_parallel, cnot, discard3, meas])
        
    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
        r = super(BlockTeleport, self).count(noiseModels, pauli, inputResult, kMax)
        return r
        
        
class CnotGadget(SequentialComponent):
    
    def __init__(self, kGood, j):
        bp1 = BP1Max(kGood)
        bpj = BP1LevelJ(kGood, bp1, j)
        code = bpj.outBlocks()[0].getCode()
        
        bp_in = ParallelComponent(kGood, bpj, bpj)
        cnot = transversal.TransCnot(kGood, code, code)

        super(CnotGadget, self).__init__(kGood, subcomponents=[bp_in, cnot])
        
class CnotExRec(SequentialComponent):
    
    def __init__(self, kGood, m, j, output_block):
        gadget = CnotGadget(kGood, j)
        discard = BlockDiscard(gadget.outBlocks(), [output_block ^ 1])
        
        code = gadget.outBlocks()[output_block].getCode()
        
        bp1 = BP1Max(kGood)
        bpj = BP1LevelJ(kGood, bp1, j)
        
        bell_pair = ParallelComponent(kGood, bpj, bpj)
        bm = BellMeas(kGood, code)
        if xType == m:
            output_block = 1
        else:
            output_block = 0
        teleport = Teleport(kGood, bell_pair, bm, output_block=output_block)
        
        super(CnotExRec, self).__init__(kGood, subcomponents=[gadget, discard, teleport])

        
        
class BP1Max(SequentialComponent):
    '''
    Component that outputs just a single block of the 1-BP.  This block
    contains syndrome counts representing the maximum over both 1-BP blocks. 
    '''
    
    def __init__(self, kGood):
        
        code = ed422.ED412Code(gaugeType=None)
        
        prepZ = Prep(kGood, ed422.prepare(Pauli.Z, Pauli.X), code)
        prepX = Prep(kGood, ed422.prepare(Pauli.X, Pauli.Z), code)
        
        bp = BellPair(kGood, prepX, prepZ, kGood)
        
        bell_meas = BellMeas(kGood, code, kGoodCnot=kGood, kGoodMeasX=kGood, kGoodMeasZ=kGood)
        teleport = TeleportED(kGood, bp, bell_meas)
        parallel_teleport = ParallelComponent(kGood, teleport, teleport)
        combine = BlockCombine(parallel_teleport.outBlocks())
        
        # Temp: remove this!
#        discard = BlockDiscard(bp.outBlocks(), [1])
#        super(BP1Max, self).__init__(kGood, subcomponents=[bp, discard])
#        super(BP1Max, self).__init__(kGood, subcomponents=[prepZ])
        super(BP1Max, self).__init__(kGood, subcomponents=[bp, parallel_teleport, combine])
        
#    def outBlocks(self):
#        return (self[-1].outBlocks()[0],)
#        
#    @fetchable
#    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
#        result = super(BP1Max, self).count(noiseModels, pauli, inputResult, kMax)
#        
#        discard1 = BlockDiscard(self[-1].outBlocks(), [0])
#        discard2 = BlockDiscard(self[-1].outBlocks(), [1])
#        
#        # Separate into two single blocks
#        counts1 = discard2.count(noiseModels, Pauli.Y, result).counts
#        counts2 = discard1.count(noiseModels, Pauli.Y, result).counts
#        
#        # Take the maximum count for each fault order k and each syndrome.
#        counts = count_errors.maxCount(counts1, counts2)
#        
#        return CountResult(counts, self.outBlocks() + result.blocks[2:])

class GaugeFilter(Filter):
    
    def __init__(self, old_code, new_code):
        self._new_code = new_code
        self._old_code = old_code
        self._parity_checks = SyndromeKeyGenerator(self._old_code, '').parityChecks()
        self._gauge_operators = []
        for ops in new_code.gaugeOperators():
            self._gauge_operators += list(ops.values())
        
        super(GaugeFilter, self).__init__()
    
    def inBlocks(self):
        return (Block('', self._old_code),)
    
    def outBlocks(self):
        return (Block('', self._new_code),)
    
    def propagateCounts(self, inputResult):
        result = super(GaugeFilter, self).propagateCounts(inputResult)
        result.blocks[0].code = self._new_code
        return result
    
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        return self.GaugelessPropagator(subPropagator, self._parity_checks, self._gauge_operators)
        
    class GaugelessPropagator(KeyManipulator):
        
        def __init__(self, manipulator, parity_checks, gauge_operators):
            super(GaugeFilter.GaugelessPropagator, self).__init__()
            gauge_operators = set(gauge_operators)
            self._gauge_bits = [i for i in range(len(parity_checks)) if parity_checks[i] in gauge_operators]
            self._n = len(parity_checks)
            
        def _manipulate(self, key):
            key_bit_list = listutils.remove_subsequence(bits.bitsToList(key[0], self._n), self._gauge_bits)
            return (bits.listToBits(key_bit_list),) + key[1:]
        
class BP1LevelJ(SequentialComponent):
    
    def __init__(self, kGood, bp1, j):
        if 2 == j:
            # Level-2 sub-block teleportation is just level-1 teleportation with postselection.
            # So we can count it precisely.
            bp = ParallelComponent(kGood, bp1, bp1)    
            block_to_teleport = j % 2
            sub = SubblockTeleport(kGood, bp, block_to_teleport=block_to_teleport, postselect=True)
        else:
            # For j=1, this is obviously correct.  For j > 2, the sub-block teleportation outputs
            # the bottom half of a level-(j-1) BP.  But logical corrections are based on level-(j-1)
            # (which is greater than one), so there is nothing else to do here.
            sub = bp1                
            
        super(BP1LevelJ, self).__init__(kGood, subcomponents=[sub])