'''
Created on 2011-10-25

@author: adam
'''
from counting.block import Block
from counting.component.base import Component, ConcatenatedComponent
from counting.component.transversal import CnotConvolver, TransCnot, TransMeas
from qec.error import Pauli
import logging

logger = logging.getLogger('component')

class BellPair(CnotConvolver):
    
    ctrlName = '|+>'
    targName = '|0>'

    def __init__(self, kGood, plus, zero, kGoodCnot):
        super(BellPair, self).__init__(kGood, kGoodCnot, plus, zero, self.ctrlName, self.targName)
    
    
class BellMeas(Component):

    cnotName = 'cnot'
    measName = 'measX,measZ'
    
    def __init__(self, kGood, code, kGoodMeasX=None, kGoodMeasZ=None, kGoodCnot=None):
        if None == kGoodMeasX: kGoodMeasX = kGood
        if None == kGoodMeasZ: kGoodMeasZ = kGood
        if None == kGoodCnot: kGoodCnot = kGood
        
        measX = TransMeas(kGoodMeasX, code, Pauli.X)
        measZ = TransMeas(kGoodMeasZ, code, Pauli.Z)
        meas = ConcatenatedComponent(kGood, measX, measZ)
        
        subs = {self.cnotName: TransCnot(kGoodCnot, code, code),
                self.measName: meas}
        
        super(BellMeas, self).__init__(kGood, subcomponents=subs)
        self.code = code
        
    def inBlocks(self):
        return self[self.cnotName].inBlocks()
        
    def outBlocks(self):
        measX = Block(self.measXName, self.code)
        measZ = Block(self.measZName, self.code)
        return (measX, measZ)

    def keyPropagator(self, keyMeta):
        cnot = self.subcomponents()[self.cnotName]
        return cnot.keyPropagator(keyMeta)
        
    def _convolve(self, results, noiseModels, pauli):
        
        cnotResult = results[self.cnotName]
        meas = self[self.measName]
        
        # Propagate the CNOT counts through the measurements.
        cnotResult.counts, cnotResult.keyMeta = meas.propagateCounts(cnotResult.counts, cnotResult.keyMeta)
        
#        measX = results[self.measXName]
#        measZ = results[self.measZName]
#        
#        measX.counts, measX.keyMeta = extendCounts(measX.counts, measX.keyMeta, blocksAfter=1)
#        measZ.counts, measZ.keyMeta = extendCounts(measZ.counts, measZ.keyMeta, blocksBefore=1)
#            
#        measX.blocks = cnot.blocks
#        measZ.blocks = cnot.blocks
            
        # Now convolve.
        return super(BellMeas, self)._convolve(results, noiseModels, pauli)
