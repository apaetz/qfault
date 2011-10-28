'''
Created on 2011-10-25

@author: adam
'''
from counting.component.base import Component
from counting.component.transversal import TransCnot
from counting.result import CountResult
from util import listutils
from counting.countParallel import convolve
from counting.component.ec import TECDecodeAdapter

class ExRec(Component):
    
    lecName = 'LEC'
    gaName = 'Ga'
    tecName = 'TEC'
    
    def __init__(self, kGood, lec, gadget, tec):
        tec = TECDecodeAdapter(tec)
        subs = {self.lecName: lec, self.gaName: gadget, self.tecName: tec}
        super(ExRec, self).__init__(kGood, subcomponents=subs)
        
    def _convolve(self, results, noiseModels, pauli):
        # The LEC should have already stripped all counts of the
        # logical error information.  That is, all LEC error keys
        # should now decode to the trivial error.  The LEC
        # can do this, because, correctness of the rectangle depends
        # only on the syndrome of its input, not on its logical state.
        # Thus, we can simply convolve the LEC output through the gate
        # gadget and then count up the decoded results of the TEC.
        
        # First, propagate the LEC through the Gadget.
        lecResult = results[self.lecName]
        ga = self[self.gaName]
        lecResult.counts, lecResult.keyMeta = ga.propagateCounts(lecResult.counts, lecResult.keyMeta)
        
        # Now convolve the Gadget and the LEC.
        tecResult = results.pop(self.tecName)
        convolved = super(ExRec, self)._convolve(results, noiseModels, pauli)
        
        # Now compute the result of propagating the LEC/CNOT input through the TEC
        # and then decoding.
        decodeCounts = convolve(convolved.counts, 
                                tecResult.counts, 
                                kMax=self.kGood[pauli], 
                                convolveFcn=self.convolveTEC,
                                splitListsInto=[2,1])
        #decodeCounts = self._convolveTEC(convolved, tecDecoder, self.kGood[pauli])
        
        
        return CountResult(decodeCounts, tecResult.keyMeta, tecResult.blocks)
    
    def _postCount(self, result, noiseModels, pauli):
        return result
        
    
#    def _convolveTEC(self, inputResult, tecDecoder, kMax):
#        
#        decodeCounts = [{}] * (kMax + 1)
#        
#        for k,countsK in enumerate(inputResult.counts):
#            for key, count in countsK.iteritems():
#                decoded = tecDecoder.getCounts(key)
#                
#                # TODO multiply by count
#            
#                for j in range(kMax - k + 1):
#                    decodeCounts[k+j] = listutils.addDicts(decodeCounts[j+k], decoded[j])
#               
#        return CountResult(decodeCounts, None, None)
    
    @staticmethod
    def convolveTEC(inputCounts, tecCounts):
        
        decodeCounts = {}
        for key, count in inputCounts.iteritems():
            decoded = tecCounts[key]
            decoded = {key: val * count for key, val in decoded.iteritems()}
            decodeCounts = listutils.addDicts(decodeCounts, decoded)
           
        return decodeCounts     

class CnotExRec(Component):
    
    lecCtrlName = 'LEC (ctrl)'
    lecTargName = 'LEC (targ)'
    tecCtrlName = 'TEC (ctrl)'
    tecTargName = 'TEC (targ)'
    cnotName = 'cnot'
    
    def __init__(self, kGood, kGoodCnot, lecCtrl, lecTarg, tecCtrl, tecTarg):
        # Construct a transversal CNOT component from the two input codes.
        ctrlCode = lecCtrl.outBlocks()[0].getCode()
        targCode = lecTarg.outBlocks()[0].getCode()
        cnot = TransCnot(kGoodCnot, ctrlCode, targCode)

        subcomponents = {
                         self.cnotName: cnot,
                         self.lecCtrlName: lecCtrl,
                         self.lecTargName: lecTarg,
                         self.tecCtrlName: tecCtrl,
                         self.tecTargName: tecTarg
                        }
        super(CnotExRec, self).__init__(kGood, subcomponents=subcomponents)
        
    def _count(self, noiseModels, pauli):
        pass
        
    def _convolve(self, results, pauli):
        
        cnot = self[self.cnotName]
        
        lecCtrlResult = results[self.lecCtrlName]
        lecTargResult = results[self.lecTargName]
        cnotResult = results[self.cnotName]
        
        # First, propagate the input results through the CNOT    
        lecCtrlResult.counts, lecCtrlResult.keyMeta = cnot.propagateCounts(lecCtrlResult.counts,
                                                             lecCtrlResult.keyMeta,
                                                             cnot.ctrlName)
        lecTargResult.counts, lecTargResult.keyMeta = cnot.propagateCounts(lecTargResult.counts,
                                                             lecTargResult.keyMeta,
                                                             cnot.targName)
        
        lecCtrlResult.blocks = cnotResult.blocks
        lecTargResult.blocks = cnotResult.blocks