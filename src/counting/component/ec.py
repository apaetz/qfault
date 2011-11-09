'''
Created on 2011-10-27

@author: adam
'''
from counting.component.base import Component, ConcatenatedComponent
from counting.convolve import convolveDict
from counting.countParallel import convolve
from counting.key import SyndromeKeyGenerator, SyndromeKeyDecoder, \
    SyndromeKeyMeta, KeyConcatenator
from counting.result import CountResult
from util import listutils, bits
from util.cache import fetchable
from qec.qecc import QeccNone
from counting.block import Block
import operator
from qec.error import Pauli, xType, zType
from counting import key

class ConcatenatedTEC(ConcatenatedComponent):
    # TODO: Using ConcatenatedComponent directly doesn't work because the counts
    # returned by TECDecodeAdapter are non-standard.
    
    def __init__(self, kGood, tec1, tec2):
        super(ConcatenatedTEC, self).__init__(kGood, tec1, tec2)
        
    def _convolve(self, results, noiseModels, pauli):
        concatenator = KeyConcatenator(results[0].keyMeta, results[1].keyMeta)
        def keyOp(key1, key2):
            return concatenator((key1, key2))
        
        def countMul(counts1, counts2):
            return convolveDict(counts1, counts2, keyOp=operator.add)
        
        def mergeOp(*tecCounts):
            if 1 == len(tecCounts): return tecCounts[0]
            
            merged = {}
            for counts in tecCounts:
                for key in counts:
                    merged[key] = listutils.addDicts(merged.get(key, {}), counts[key])
                    
            return merged
        
        convolved = convolve(results[0].counts, 
                             results[1].counts, 
                             kMax=self.kGood[pauli], 
                             extraArgs=[keyOp, countMul, listutils.addDicts, {}],
                             splitListsInto=[1,1],
                             listMergeOp=mergeOp)
        
        convRejected = convolve(results[0].rejected, 
                             results[1].rejected, 
                             kMax=self.kGood[pauli],
                             extraArgs=[keyOp, countMul, listutils.addDicts, {}],
                             splitListsInto=[1,1],
                             listMergeOp=mergeOp)
        
        
        return CountResult(convolved, concatenator.meta(), self.outBlocks(), rejectedCounts=convRejected)
        
#    def _convolveTEC(self, tecCounts1, tecCounts2, meta):
#        # TEC counts are indexed by [key][pauli]
#        # which gives the weighted count for the
#        # logical error 'pauli', given an input
#        # error of 'key'.
#        
#        concatenator = KeyConcatenator(meta)
#        
#        convolved = {}
#        for key1, counts1 in results1.iteritems():
#            for key2, counts2 in results2.iteritems():
#                key = concatenator((key1,key2))
#                counts = convolve.convolveDict(counts1, counts2)
#                convolved[key] = listutils.addDicts(convolved.get(key, {}), counts)
#                
#        return convolved

class TECDecodeAdapter(Component):
    
    def __init__(self, tec):
        self._tec = tec
        
    def count(self, noiseModels, pauli):
        # This must return a count result for which the counts
        # are indexable by k (the number of faults) and then
        # an input key.  The value of [k][key] is a dictionary
        # indexed by (logical) Pauli error.  Finally, the
        # value [k][key][pauli] is a weighted count of the
        # ways to produce the logical error 'pauli' given
        # k failures in the TEC when the input to the TEC
        # is 'key'.
        
        counts, rejected, meta, blocks = self.lookupTable(self._tec, noiseModels, pauli)
        return CountResult(counts, meta, blocks, rejectedCounts=rejected)
        
    def outBlocks(self):
        code = QeccNone(1)
        return tuple(Block(block.name, code) for block in self._tec.outBlocks())
        
    @staticmethod
    @fetchable
    def lookupTable(tec, noiseModels, pauli):
        code = tec.inBlocks()[0].getCode()
        keyMeta = SyndromeKeyGenerator(code, '').keyMeta()
        nchecks = len(keyMeta.parityChecks())
        
        if Pauli.X == pauli:
            keyMask = ((0 != check[zType]) for check in keyMeta.parityChecks())
        elif Pauli.Z == pauli:
            keyMask = ((0 != check[xType]) for check in keyMeta.parityChecks())
        else:
            keyMask = (1 for _ in keyMeta.parityChecks())
            
        keyMask = bits.listToBits(keyMask)
        keys = set(key & keyMask for key in xrange(1 << nchecks))
        
        decoder = SyndromeKeyDecoder(code)
        countLookup = [{} for _ in range(tec.kGood[pauli] + 1)]
        rejectLookup = [{} for _ in range(tec.kGood[pauli] + 1)]
        for key in keys:
            inCount = [{(key,): 1}]
            inResult = CountResult(inCount, keyMeta, [None])
            outResult = tec.count(noiseModels, pauli, inResult)
            outMeta = outResult.keyMeta
            outBlocks = outResult.blocks
            
            dCounts = TECDecodeAdapter.decodeCounts(outResult.counts, decoder)
            for k in range(len(countLookup)):
                countLookup[k][(key,)] = dCounts[k]
                rejectLookup[k][(key,)] = outResult.rejected[k]
            
        return countLookup, rejectLookup, outMeta, outBlocks
            
    @staticmethod
    def decodeCounts(counts, decoder):
        decodeCounts = []
        for countK in counts:
            decodeCountK = {}
            for key, val in countK.iteritems():
                decoded = decoder.decode(key)
                decodeCountK[decoded] = decodeCountK.get(decoded, 0) + val
                    
            decodeCounts.append(decodeCountK)
            
        return decodeCounts