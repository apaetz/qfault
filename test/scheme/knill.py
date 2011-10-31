'''
Created on 2011-10-30

@author: adam
'''
from counting.component import base
from counting.component.base import Prep, InputAdapter, ConcatenatedComponent
from counting.component.bell import BellPair, BellMeas
from counting.component.ec import TECDecodeAdapter, ConcatenatedTEC
from counting.component.exrec import ExRec
from counting.component.teleport import TeleportED
from counting.component.transversal import TransCnot
from qec import ed422, qecc, error, error
from qec.error import Pauli
from qec.qecc import StabilizerState
from settings.noise import CountingNoiseModelX, CountingNoiseModelZ
import unittest


class TestKnill(unittest.TestCase):

    def makeED(self, kGood):
        
        prepZ = Prep(kGood, 
                          ed422.prepare(Pauli.Z, Pauli.X), 
                          StabilizerState(ed422.ED412Code(gaugeType=error.xType), [error.zType]))
        prepX = Prep(kGood, 
                          ed422.prepare(Pauli.X, Pauli.Z), 
                          StabilizerState(ed422.ED412Code(gaugeType=error.xType), [error.xType]))
        
        bellPair = BellPair(kGood, prepX, prepZ, kGood)
        bellMeas = BellMeas(kGood, ed422.ED412Code(), kGood, kGood, kGood)
        
        teleportED = TeleportED(kGood, bellPair, bellMeas)
        
        return teleportED


    def testIdentityExRec(self):
        kGood = {pauli: 1 for pauli in [Pauli.X, Pauli.Z, Pauli.Y]}
    
    
        noises = {Pauli.X: CountingNoiseModelX(),
                  Pauli.Z: CountingNoiseModelZ(),
                  Pauli.Y: None,
                  }
        
        pauli = Pauli.X
        code = ed422.ED412Code()
    
        id = base.Empty(code)
        
        ed = self.makeED(kGood)
        led = InputAdapter(ed, (0,))
        ted = TECDecodeAdapter(ed)
        
        exRec = ExRec(kGood, led, id, ted)
        
        result = exRec.count(noises, pauli)
        print result.counts



if __name__ == "__main__":
    
    import util.cache
    util.cache.enableFetch(False)
    util.cache.enableMemo(False)
    
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()