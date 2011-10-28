'''
Created on 2011-10-27

@author: adam
'''
from counting.key import SyndromeKeyGenerator, SyndromeKeyDecoder,\
    SyndromeKeyMeta
from counting.result import CountResult
from util.cache import fetchable
from counting.component.base import Component



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
        
        counts, meta, blocks = self.lookupTable(self._tec, noiseModels, pauli)
        return CountResult(counts, meta, blocks)
        
    @staticmethod
    @fetchable
    def lookupTable(tec, noiseModels, pauli):
        code = tec.inBlocks()[0].getCode()
        keyMeta = SyndromeKeyGenerator(code, '').keyMeta()
        nchecks = len(keyMeta.parityChecks())
        decoder = SyndromeKeyDecoder(code)
        lookup = [{} for _ in range(tec.kGood[pauli] + 1)]
        for key in xrange(1 << nchecks):
            inCount = [{(key,): 1}]
            inResult = CountResult(inCount, keyMeta, None)
            outResult = tec.count(noiseModels, pauli, inResult)
            outMeta = outResult.keyMeta
            outBlocks = outResult.blocks
            
            dCounts = TECDecodeAdapter.decodeCounts(outResult.counts, decoder)
            for k in range(len(lookup)):
                lookup[k][(key,)] = dCounts[k]
            
        return lookup, outMeta, outBlocks
            
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