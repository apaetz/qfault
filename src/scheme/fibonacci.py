'''
Created on 2012-04-02

@author: adam
'''

from counting.block import Block
from counting.component import transversal
from counting.component.base import ParallelComponent, SequentialComponent, Prep, \
    ConcatenationFilter
from counting.component.bell import BellPair, BellMeas
from counting.component.block import BlockDiscard, BlockPermutation
from counting.component.teleport import Teleport, TeleportED
from counting.countErrors import error, mapCounts, maxCount
from counting.key import IdentityManipulator, KeyMerger, SyndromeKeyGenerator
from qec import ed422
from qec.error import Pauli
from qec.qecc import StabilizerState, ConcatenatedCode
from scheme import Scheme
from counting.component.adapter import IdealDecoder

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
        teleport = TeleportED(kBPT, bellPair, bellMeas)
        
        self.bp = bellPair
        self.bpt = BellPairTeleport(kBPT, bellPair, teleport)
        self.bp2 = BellPair(kBPT, PrepLevel2(kBPT, self.bpt, Pauli.X), PrepLevel2(kBPT, self.bpt, Pauli.Z), kBPT)
        self.bpt2 = BellPairTeleport(kBPT, self.bp2, teleport)
        
        permutations = ([0,2,1,3], [0,3,1,2], [1,2,0,3], [1,3,0,2])
    
        self.sbtList = [SubblockTeleport(kBPT, self.bpt, perm) for perm in permutations]
        
    def count(self):
#        return self.bpt2.count(self.defaultNoiseModels, Pauli.Y)

#        results = [sbt.count(self.defaultNoiseModels, Pauli.Y) for sbt in self.sbtList]
#        countss = [result.counts for result in results]
#        print countss
#        countsMax = maxCount(*countss)
#        print countsMax
#        
#        return countsMax
#    
        bptDecode = BPTLevelOneSingleBlock(self.bpt.kGood, self.bpt)
        result = bptDecode.count(self.defaultNoiseModels, Pauli.Y)
        print result.counts
        
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
        
class SubblockTeleport(SequentialComponent):
    
    def __init__(self, kGood, bellPair, permutation):
        code = bellPair.outBlocks()[0].getCode()
        
        twinBP = ParallelComponent(kGood, bellPair, bellPair)
        perm = BlockPermutation(twinBP.outBlocks(), permutation)
        discard2 = BlockDiscard(perm.outBlocks(), 2)
        cnot = transversal.TransCnot(kGood, code, code)
        discard1 = BlockDiscard(cnot.outBlocks(), 1)
        bm = BellMeas(kGood, code, kGood, kGood, kGood)
        teleport = Teleport(kGood, bellPair, bm)
        
        super(SubblockTeleport, self).__init__(kGood, subcomponents=[twinBP, perm, discard2, cnot, discard1, teleport])    
        
        
class BPTLevelOneSingleBlock(SequentialComponent):
    
    def __init__(self, kGood, bpt):
        discard = BlockDiscard(bpt.outBlocks(), 1)
        decode = IdealDecoder(discard.outBlocks()[0].getCode())
        
        super(BPTLevelOneSingleBlock, self).__init__(kGood, subcomponents=[bpt, discard, decode])
        
if __name__ == '__main__':
    import util.cache
    util.cache.enableFetch(False)
#    util.cache.enableMemo(False)
    
    from counting import countParallel
    countParallel.setPool(countParallel.DummyPool())
    
    kBPT = {Pauli.Y: 2}
    scheme = FibonacciScheme(kBPT)
    
    count = scheme.count()
#    print count.counts
#    print count.blocks