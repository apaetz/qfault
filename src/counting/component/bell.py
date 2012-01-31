'''
Created on 2011-10-25

Components for preparing Bell states and making Bell measurements.

@author: adam
'''
from counting.component.base import ParallelComponent,\
    SequentialComponent
from counting.component.transversal import TransCnot, TransMeas
from qec.error import Pauli
import logging

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

        # Order is important here.  The preps must be in front of the CNOT.
        super(BellPair, self).__init__(kGood, subcomponents=(prep, cnot))
    
    
class BellMeas(SequentialComponent):
    '''
    Performs a transversal measurement in the Bell basis.
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
