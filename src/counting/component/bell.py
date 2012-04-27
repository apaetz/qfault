'''
Created on 2011-10-25

Components for preparing Bell states and making Bell measurements.

@author: adam
'''
from counting.component.base import ParallelComponent,\
    SequentialComponent, Filter
from counting.component.transversal import TransCnot, TransMeas
from qec.error import Pauli, xType, zType
import logging
from counting.key import SyndromeKeyGenerator, MultiBlockSyndromeKeyGenerator,\
    KeyManipulator, IdentityManipulator
from counting.block import Block
import counting

logger = logging.getLogger('component')

class BellPair(SequentialComponent):
    '''
    Prepares the encoded Bell state |+>|+>.
    '''

    def __init__(self, kGood, plus, zero, kGoodCnot):
        # Construct a transversal CNOT component from the two input codes.
        ctrlCode = plus.outBlocks()[0].getCode()
        targCode = zero.outBlocks()[0].getCode()
        cnot = TransCnot(kGoodCnot, ctrlCode, targCode)
        
        prep = ParallelComponent(kGood, plus, zero)
        bf = BellFilter(cnot.outBlocks()[0].getCode())

        # Order is important here.  The preps must be in front of the CNOT.
        super(BellPair, self).__init__(kGood, subcomponents=(prep, cnot, bf))
        
        
class BellFilter(Filter):
    '''
    Logical XX and ZZ errors act trivially on Bell pairs.  This component
    maps logical XX and ZZ syndrome counts to trivial syndrome counts.
    '''
    
    def __init__(self, code):
        super(BellFilter, self).__init__()
        
        self.code = code
        logicals = code.logicalOperators()[0]
        blocks = self.inBlocks()
        generator = MultiBlockSyndromeKeyGenerator(blocks)
        xx = {blocks[0].name: logicals[xType], blocks[1].name: logicals[xType]}
        zz = {blocks[0].name: logicals[zType], blocks[1].name: logicals[zType]}
        xxKey = generator.getKey(xx)
        zzKey = generator.getKey(zz)
        
        self.keys = (xxKey, zzKey)
        self.trivialKey = generator.getKey({str(blocks[0]): Pauli.I, str(blocks[1]): Pauli.I})
        
    def inBlocks(self):
        return (Block("A", self.code), Block("B", self.code))
        
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        return self.LogicalPropagator(subPropagator, self.keys, self.trivialKey)
        
    class LogicalPropagator(KeyManipulator):
        
        def __init__(self, manipulator, logicalKeys, trivialKey):
            self.logicalKeys = logicalKeys
            self.trivialKey = trivialKey
            super(BellFilter.LogicalPropagator, self).__init__(manipulator)
            
        def _manipulate(self, key):
            if key in self.logicalKeys:
                return self.trivialKey
            return key
    
class BellMeas(SequentialComponent):
    '''
    Performs a transversal measurement in the Bell basis.
    Output block 0 is the result of X-basis measurement.
    Output block 1 is the result of Z-basis measurement.
    '''
    
    def __init__(self, kGood, code, kGoodMeasX=None, kGoodMeasZ=None, kGoodCnot=None):
        if None == kGoodMeasX: kGoodMeasX = kGood
        if None == kGoodMeasZ: kGoodMeasZ = kGood
        if None == kGoodCnot: kGoodCnot = kGood
        
        measX = TransMeas(kGoodMeasX, code, Pauli.X)
        measZ = TransMeas(kGoodMeasZ, code, Pauli.Z)
        meas = ParallelComponent(kGood, measX, measZ)
        
        subs = (TransCnot(kGoodCnot, code, code), meas)
        
        super(BellMeas, self).__init__(kGood, subcomponents=subs)
        self.code = code
