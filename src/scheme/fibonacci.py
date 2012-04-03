'''
Created on 2012-04-02

@author: adam
'''

from counting.block import Block
from counting.component.base import ParallelComponent, SequentialComponent, Prep,\
    ConcatenationFilter, BlockPermutation
from counting.component.bell import BellPair, BellMeas
from counting.component.teleport import Teleport
from counting.countErrors import error, mapCounts
from counting.key import IdentityManipulator, KeyMerger, SyndromeKeyGenerator
from qec import ed422
from qec.error import Pauli
from qec.qecc import StabilizerState, ConcatenatedCode
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
        self.bp2 = BellPair(kBPT, PrepLevel2(kBPT, self.bpt, Pauli.X), PrepLevel2(kBPT, self.bpt, Pauli.Z), kBPT)
        self.bpt2 = BellPairTeleport(kBPT, self.bp2, teleport)
        
    def count(self):
        return self.bpt2.count(self.defaultNoiseModels, Pauli.Y)
    
    
class BellPairTeleport(SequentialComponent):
    
    def __init__(self, kGood, bellPair, teleport):
        parallelTeleport = ParallelComponent(kGood, teleport, teleport)
        super(SequentialComponent, self).__init__(kGood, subcomponents=[bellPair, parallelTeleport])
        

class PrepLevel2(SequentialComponent):
    
    def __init__(self, kGood, bellPair, eigenstate=Pauli.Z):
        twinBP = ParallelComponent(kGood, bellPair, bellPair)
        subs = [twinBP]
        
        if (Pauli.Z == eigenstate):
            # |0> is prepared by permuting the second and third logical blocks.
            perm = BlockPermutation(twinBP.outBlocks(), [0,2,1,3])
            subs.append(perm)
            
        code = subs[-1].outBlocks()[0].getCode()
        subs.append(ConcatenationFilter(code, code))
            
        super(SequentialComponent, self).__init__(kGood, subcomponents=subs)

#    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
#        result = super(PrepLevel2, self).count(noiseModels, pauli, inputResult, kMax)
#        result = self.propagateCounts(result)
#        return result
#
#    def outBlocks(self):
#        '''
#        Output a single level-2 block rather than four level-1 blocks.
#        '''
#        subblocks = super(PrepLevel2, self).outBlocks()
#        code = subblocks[0].getCode()
#        catCode = ConcatenatedCode(code, code)
#        return tuple([Block('2-Prep', catCode)])
#
#
#    def keyPropagator(self, subPropagator=IdentityManipulator()):
#        '''
#        Convert the four level-1 blocks into a single level-2 block.
#        '''
#        subPropagator = super(PrepLevel2, self).keyPropagator(subPropagator)
#        subblocks = super(PrepLevel2, self).outBlocks()
#        keyLengths = [len(SyndromeKeyGenerator(block.getCode(), None).parityChecks()) for block in subblocks]
#
#        return KeyMerger(subPropagator, keyLengths)

    
        
        
if __name__ == '__main__':
    import util.cache
    util.cache.enableFetch(False)
#    util.cache.enableMemo(False)
    
    from counting import countParallel
    countParallel.setPool(countParallel.DummyPool())
    
    kBPT = {Pauli.Y: 2}
    scheme = FibonacciScheme(kBPT)
    
    count = scheme.count()
    print count.counts
    print count.blocks