'''
Created on 2011-11-14

@author: adam
'''
from counting.component.adapter import InputAdapter
from counting.component.base import Prep, Concatenator, FixedOutput
from counting.component.bell import BellPair, BellMeas
from counting.component.ec import TECDecodeAdapter, ConcatenatedTEC
from counting.component.exrec import ExRecForward
from counting.component.teleport import TeleportED
from counting.component.transversal import TransCnot
from qec import ed422, error
from qec.error import Pauli, PauliError
from qec.qecc import StabilizerState
from scheme import Scheme
from counting.key import MultiBlockSyndromeKeyGenerator
import logging

class KnillScheme(Scheme):
    '''
    classdocs
    '''


    def __init__(self, kPrep, kCnot, kEC, kExRec):
        '''
        Constructor
        '''
        
        prepZ = Prep(kPrep, 
                      ed422.prepare(Pauli.Z, Pauli.X), 
                      StabilizerState(ed422.ED412Code(gaugeType=error.xType), [error.zType]))
        prepX = Prep(kPrep, 
                          ed422.prepare(Pauli.X, Pauli.Z), 
                          StabilizerState(ed422.ED412Code(gaugeType=error.xType), [error.xType]))
        
        bellPair = BellPair(kEC, prepX, prepZ, kCnot)
        bellMeas = BellMeas(kEC, ed422.ED412Code(), kGoodCnot=kCnot, kGoodMeasX=kCnot, kGoodMeasZ=kCnot)
        
        ed = TeleportED(kEC, bellPair, bellMeas)
        
        self.ed = ed
        self.kExRec = kExRec
        self.kCnot = kCnot
        
    def getExRecs(self):
        code = ed422.ED412Code(gaugeType=error.xType)
        keyGen = MultiBlockSyndromeKeyGenerator(self.ed.inBlocks())
        blockname = self.ed.inBlocks()[0].name
        
        inKeys = set()
        for eX in range(1 << 4):
            for eZ in range(1 << 4):
                pauliError = PauliError(xbits=eX, zbits=eZ)
                inKeys.add(keyGen.getKey({blockname: pauliError}))
        
        exrecs = {}
        
        # Generate CNOT exRecs for all possible configurations of the leading EDs.
        for haveLEDa, haveLEDb in [#(False, False),
                                   #(False, True),
                                   #(True, False),
                                   (True, True)]: 
            if haveLEDa:
                ledsA = [InputAdapter(self.ed, keyGen.getKey({blockname: Pauli.I}))]
            else:
                ledsA = [FixedOutput(code, [{key:1}], keyGen.keyMeta()) for key in inKeys]
                
            if haveLEDb:
                ledsB = [InputAdapter(self.ed, keyGen.getKey({blockname: Pauli.I}))]
            else:
                ledsB = [FixedOutput(code, [{key:1}], keyGen.keyMeta()) for key in inKeys]
                
            leds = []
            for ledA in ledsA:
                for ledB in ledsB:
                    led = Concatenator(self.kExRec, ledA, ledB)
                    leds.append(led)
                    
            cnot = TransCnot(self.kCnot, code, code)
            ted = TECDecodeAdapter(self.ed)
            ted = ConcatenatedTEC(self.kExRec, ted, ted)
        
            exrecs['cnot'] = exrecs.get('cnot', []) + [ExRecForward(self.kExRec, led, cnot, ted) for led in leds]

        return exrecs
    
if __name__ == '__main__':
    from counting import countParallel
    countParallel.setPool(countParallel.DummyPool())
    
    logging.getLogger('counting.threshold').setLevel(logging.DEBUG)
    
    kPrep = {Pauli.Y: 2}
    kCnot = {Pauli.Y: 2}
    kEC = {Pauli.Y: 4}
    kExRec = {Pauli.Y: 11}
    scheme = KnillScheme(kPrep, kCnot, kEC, kExRec)
    scheme.count()