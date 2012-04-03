'''
Created on 2012-04-02

@author: adam
'''

from counting.component.base import ParallelComponent, \
    SequentialComponent, Prep
from counting.component.bell import BellPair, BellMeas
from counting.component.teleport import Teleport
from counting.countErrors import error
from qec import ed422
from qec.error import Pauli
from qec.qecc import StabilizerState
from scheme import Scheme

class FibonacciScheme(Scheme):
    '''
    classdocs
    '''


    def __init__(self, kBPT):
        '''
        Constructor
        '''
        self.code = ed422.ED412Code(gaugeType=error.xType)
        
        prepZ = Prep(kBPT, ed422.prepare(Pauli.Z, Pauli.X), StabilizerState(self.code, [error.zType]))
        prepX = Prep(kBPT, ed422.prepare(Pauli.X, Pauli.Z), StabilizerState(self.code, [error.xType]))
        
        bellPair = BellPair(kBPT, prepX, prepZ, kBPT)
        bellMeas = BellMeas(kBPT, self.code, kGoodCnot=kBPT, kGoodMeasX=kBPT, kGoodMeasZ=kBPT)
        teleport = Teleport(kBPT, bellPair, bellMeas)
        
        self.bpt = BellPairTeleport(kBPT, bellPair, teleport)
        
    def count(self):
        return self.bpt.count(self.defaultNoiseModels, Pauli.Y)
    
    
class BellPairTeleport(SequentialComponent):
    
    def __init__(self, kGood, bellPair, teleport):
        parallelTeleport = ParallelComponent(kGood, teleport, teleport)
        super(SequentialComponent, self).__init__(kGood, subcomponents=[bellPair, parallelTeleport])
        
        
if __name__ == '__main__':
    from counting import countParallel
    countParallel.setPool(countParallel.DummyPool())
    
    kBPT = {Pauli.Y: 2}
    scheme = FibonacciScheme(kBPT)
    
    count = scheme.count()
    print count.counts