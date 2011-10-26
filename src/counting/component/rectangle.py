'''
Created on 2011-10-26

@author: adam
'''
from counting.component.base import InputDependentComponent
from counting.component.transversal import TransCnot
from counting.key import keySplitter, keyForBlock, keyConcatenator,\
    SyndromeKeyMeta, SyndromeKeyGenerator, SyndromeKeyDecoder
from counting.countErrors import mapCounts
from counting.result import CountResult
from util.cache import fetchable, memoizeFetchable
from counting.countParallel import convolve
from util import listutils

class CnotRectangle(InputDependentComponent):
    '''
    classdocs
    '''

    cnotName = 'cnot'
    tecName = 'TEC'

    def __init__(self, kGood, kGoodCnot, tec):
        # Construct a transversal CNOT component from the two input codes.
        code = tec.outBlocks()[0].getCode()
        cnot = TransCnot(kGoodCnot, code, code)

        subcomponents = {
                         self.cnotName: cnot,
                         self.tecName: tec
                        }
        super(CnotRectangle, self).__init__(kGood, subcomponents=subcomponents)
        
#    def _propagateInput(self, inputResult):
#        cnot = self[self.cnotName]
#        ucTeleport = self[self.ucTeleportName]
#        
#        cnotPropagator, cnotMeta = cnot.keyPropagator(inputResult.keyMeta)
#        ucTeleportPropagator, ucMeta = ucTeleport.keyPropagator(cnotMeta)
#        ucConcatenator, outMeta = keyConcatenator(ucMeta, ucMeta)
#        
#        def propagate(key):
#            cnotKey = cnotPropagator(key)
#            key0 = keyForBlock(cnotKey, 0, cnotMeta)
#            key1 = keyForBlock(cnotKey, 1, cnotMeta)
#            key0 = ucTeleportPropagator(key0)
#            key1 = ucTeleportPropagator(key1)
#            return ucConcatenator(key0, key1)
#        
#        propagated = mapCounts(inputResult.counts, propagate)
#        
#        return propagated, outMeta

    def _propagateInput(self, inputResult):
        cnot = self[self.cnotName]
        counts, keyMeta = cnot.propagateCounts(inputResult.counts, inputResult.keyMeta)
        
        return counts, keyMeta

    def _count(self, noiseModels, pauli):
        # TODO: also initialize TEC decoder here?
        return {self.cnotName: self[self.cnotName].count(noiseModels, pauli)}
        
    def _convolve(self, results, noiseModels, pauli, inputResult):
              
        # Convolve the input with the transversal CNOT.
        convolved = super(CnotRectangle, self)._convolve(results, noiseModels, pauli, inputResult)
        
        # Split the keys into two parts; one for each TEC
        cnotResult = results[self.cnotName]
        splitter, _, _ = keySplitter(cnotResult.keyMeta, 1)
        convolved.counts = mapCounts(convolved.counts, splitter)
        
        decodeCounter = TECDecoder(self[self.tecName], noiseModels, pauli)
        
        decodeCounts = []
        for k,countsK in enumerate(convolved.counts):
            for splitKey, count in countsK.iteritems():
                # Use the decodeCounter to determine the counts for 
                # each TEC given the input key.
                decodeCounts0 = decodeCounter.getCounts(splitKey[0])
                decodeCounts1 = decodeCounter.getCounts(splitKey[1])
                
                # Combine the two counts to get a single count over two blocks.
                # The result is a list of counts indexed by k, the number of failures.
                # The value convolved[k] is a dictionary indexed by 
                # 2-qubit (logical) errors. 
                convolved = convolve(decodeCounts0, decodeCounts1, kMax=self.kGood[pauli] - k)
                
                # TODO multiply by count
            
                for j in len(convolved):
                    decodeCounts[j+k] = listutils.addDicts(decodeCounts[j+k], convolved)
               
        return CountResult(decodeCounts, None, None)
        
class TECDecoder(object):
    # TODO: this could probably be a component.
    
    def __init__(self, tec, noiseModels, pauli):
        self._tec = tec
        self._lookup = self.lookupTable(tec, noiseModels, pauli)
        
    @staticmethod
    @fetchable
    def lookupTable(tec, noiseModels, pauli):
        code = tec.inBlocks()[0].getCode()
        keyMeta = SyndromeKeyGenerator(code, '').keyMeta()
        nchecks = len(keyMeta.parityChecks())
        decoder = SyndromeKeyDecoder(code)
        lookup = {}
        for key in xrange(1 << nchecks):
            inCount = [{(key,): 1}]
            inResult = CountResult(inCount, keyMeta, None)
            outResult = tec.count(noiseModels, pauli, inResult)
            dCounts = TECDecoder.decodeCounts(outResult.counts, decoder)
            lookup[key] = dCounts
            
        return lookup
            
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
    
    def getCounts(self, key):
        return self._lookup[key]
        
        