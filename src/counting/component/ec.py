'''
Created on 2011-10-27

@author: adam
'''
from copy import copy
from counting import key
from counting.block import Block
from counting.component.adapter import ComponentAdapter, IdealDecoder
from counting.component.base import Component, ParallelComponent, \
    SequentialComponent
from counting.convolve import convolveDict
from counting.countErrors import mapCounts
from counting.countParallel import convolve
from counting.key import SyndromeKeyGenerator, SyndromeKeyDecoder, \
    KeyConcatenator, KeyMasker, IdentityManipulator
from counting.result import CountResult
from qec.error import Pauli, xType, zType
from qec.qecc import QeccNone
from util import listutils, bits
from util.cache import fetchable
from util.patterns import DelegatingAdapter
import logging
import operator

class ConcatenatedTEC(ParallelComponent):
    # TODO: Using ParallelComponent directly doesn't work because the counts
    # returned by DecodeAdapter are non-standard.
    
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
        
        
        return CountResult(convolved, self.outBlocks(), rejectedCounts=convRejected)
        
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
    
class TECAdapter(ComponentAdapter):
    
    def __init__(self, tec):
        super(TECAdapter, self).__init__(tec)
        self.tec = tec
    
    @fetchable
    def lookupTable(self, noiseModels, pauli):
        code = self.tec.inBlocks()[0].getCode()
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
        
        countLookup = [{} for _ in range(self.tec.kGood[pauli] + 1)]
        rejectLookup = [{} for _ in range(self.tec.kGood[pauli] + 1)]
        for key in keys:
            inCount = [{(key,): 1}]
            inResult = CountResult(inCount, [None])
            outResult = self.tec.count(noiseModels, pauli, inResult)
            outMeta = outResult.keyMeta
            outBlocks = outResult.blocks
            
            counts = self._processCounts(outResult.counts)
            for k in range(len(countLookup)):
                countLookup[k][(key,)] = counts[k]
                rejectLookup[k][(key,)] = outResult.rejected[k]
            
        return countLookup, rejectLookup, outMeta, outBlocks
        
    def count(self, noiseModels, pauli):
        # This must return a count result for which the counts
        # are indexable by k (the number of faults) and then
        # an input key.  The value of [k][key] is a dictionary
        # indexed by (logical) Pauli error.  Finally, the
        # value [k][key][pauli] is a weighted count of the
        # ways to produce the logical error 'pauli' given
        # k failures in the TEC when the input to the TEC
        # is 'key'.
        
        counts, rejected, meta, blocks = self.lookupTable(noiseModels, pauli)
        return CountResult(counts, blocks, rejectedCounts=rejected)
        
    def outBlocks(self):
        code = QeccNone(1)
        return tuple(Block(block.name, code) for block in self.tec.outBlocks())
    
    def _processCounts(self, counts):
        return counts
    
    
class DecodeAdapter(SequentialComponent):
    
    def __init__(self, tec):
        outBlocks = tec.outBlocks()
        idealDecoders = [IdealDecoder(block.getCode()) for block in outBlocks]
        decoder = ParallelComponent({}, *idealDecoders)
        super(DecodeAdapter, self).__init__(tec.kGood, (tec, decoder))
    
#    @staticmethod
#    def decodeCounts(counts, decoder):
#        decodeCounts = []
#        for countK in counts:
#            decodeCountK = {}
#            for key, val in countK.iteritems():
#                decoded = decoder.decode(key)
#                decodeCountK[decoded] = decodeCountK.get(decoded, 0) + val
#                    
#            decodeCounts.append(decodeCountK)
#            
#        return decodeCounts
#    
#    def _processCounts(self, counts):
#        code = self.tec.inBlocks()[0].getCode()
#        decoder = SyndromeKeyDecoder(code)
#        return self.decodeCounts(counts, decoder)
            
            
class LECSyndromeAdapter(ComponentAdapter):
    
    def __init__(self, lec):
        super(LECSyndromeAdapter, self).__init__(lec)
        self.lec = lec
    
    def logicalMasker(self, keyMeta):        
        # Strip off the logical syndrome information.
        block = self.lec.outBlocks()[0]
        stabilizers = set(block.getCode().stabilizerGenerators())   
        parityChecks = keyMeta.parityChecks() 
        syndromeBits = [(check in stabilizers) for check in parityChecks]
        syndromeMask = bits.listToBits(syndromeBits)
        
        masker = KeyMasker(keyMeta, syndromeMask)
        
        return masker
    
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        return self.logicalMasker(self.lec.keyPropagator(subPropagator))
    
    def _postCount(self, result, noiseModels, pauli):
        result = self.lec._postCount(result, noiseModels, pauli)
        self._log(logging.DEBUG, 'Stripping logical syndrome information.')
        masker = self.logicalMasker(result.keyMeta)
        result.counts = mapCounts(result.counts, masker)
        result.keyMeta = masker.meta()
        return result
