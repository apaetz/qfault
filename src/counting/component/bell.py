'''
Created on 2011-10-25

@author: adam
'''
from counting.block import Block
from counting.component.base import Component
from counting.component.transversal import CnotConvolver, TransCnot, TransMeas
from counting.countErrors import extendCounts
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
    measXName = 'measX'
    measZName = 'measZ'
    
    def __init__(self, kGood, code, kGoodMeasX=None, kGoodMeasZ=None, kGoodCnot=None):
        if None == kGoodMeasX: kGoodMeasX = kGood
        if None == kGoodMeasZ: kGoodMeasZ = kGood
        if None == kGoodCnot: kGoodCnot = kGood
        
        subs = {self.cnotName: TransCnot(kGoodCnot, code, code),
                self.measXName: TransMeas(kGoodMeasX, code, Pauli.X, self.measXName),
                self.measZName: TransMeas(kGoodMeasZ, code, Pauli.Z, self.measZName)}
        
        super(BellMeas, self).__init__(kGood, subcomponents=subs)
        self.code = code
        
    def outBlocks(self):
        measX = Block(self.measXName, self.code)
        measZ = Block(self.measZName, self.code)
        return (measX, measZ)

    def keyPropagator(self, keyMeta, blockname):
        cnot = self.subcomponents()[self.cnotName]
        namemap = {self.measXName: cnot.ctrlName, self.measZName: cnot.targName}
        return cnot.keyPropagator(keyMeta, namemap[blockname])
        
    def _convolve(self, results, noiseModels, pauli):
        
        cnot = results[self.cnotName]
        measX = results[self.measXName]
        measZ = results[self.measZName]
        
        measX.counts, measX.keyMeta = extendCounts(measX.counts, measX.keyMeta, blocksAfter=1)
        measZ.counts, measZ.keyMeta = extendCounts(measZ.counts, measZ.keyMeta, blocksBefore=1)
            
        measX.blocks = cnot.blocks
        measZ.blocks = cnot.blocks
            
        # Now convolve.
        return super(BellMeas, self)._convolve(results, noiseModels, pauli)
