'''
Created on 2011-10-25

@author: adam
'''
from counting import key, probability
from counting.component.base import Component, SequentialComponent
from counting.component.transversal import TransCnot
from counting.countParallel import convolve
from counting.result import CountResult
from util import listutils
from qec.error import Pauli
from counting.block import Block
from qec import qecc
from counting.key import IdentityManipulator, SyndromeKeyDecoder, KeyManipulator
from util.cache import fetchable

class ExRec(SequentialComponent):
    
#    lecName = 'LEC'
#    gaName = 'Ga'
#    tecName = 'TEC'
    
    def __init__(self, kGood, lec, gadget, tec):
        subs = (lec, gadget, tec)
        super(ExRec, self).__init__(kGood, subcomponents=subs)
        
    @fetchable
    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
        return super(ExRec, self).count(noiseModels, pauli, inputResult, kMax)
                
#    def prAccept(self, noiseModels, kMax=None):
#        prLEC = self[self.lecName].prAccept(noiseModels, self.kGood)
#        prGa = self[self.gaName].prAccept(noiseModels, self.kGood)
#        
#        # There are two cases. Either the TEC depends on the input,
#        # or it does not.  Assume first that it does not depend on
#        # the input and ask for its acceptance probability directly.
#        try:
#            prTEC = self[self.tecName].prAccept(noiseModels, kMax=self.kGood)
#        except TypeError:
#            # The TEC depends on the input.  We may calculate the 
#            # acceptance probability from the convolution of the
#            # TEC rejected counts and the LEC-Ga (accepted) counts which
#            # is calculated by the count() method.  Effectively we
#            # have pushed acceptance from the TEC up to the entire
#            # exRec.
#            #prTEC = super(ExRec, self).self.prAccept(noiseModels, kMax=self.kGood)
#            
#            # TODO: hardcoded for now.  need to parameterize and probably refactor
#            # to a different class.
#            pauli = Pauli.Y
#            rejected = self.count(noiseModels, pauli).rejected
#            prBad = self.prBad(noiseModels[pauli], pauli, kMax)
#            locTotals = self.locations(pauli).getTotals()
#
#            prTEC = 1 - probability.upperBoundPoly(rejected, prBad, locTotals, noiseModels[pauli])
#        
#        return prLEC * prGa * prTEC

#    def _convolve(self, results, noiseModels, pauli, inResult):
#        # The LEC should have already stripped all counts of the
#        # logical error information.  That is, all LEC error keys
#        # should now decode to the trivial error.  The LEC
#        # can do this, because, correctness of the rectangle depends
#        # only on the syndrome of its input, not on its logical state.
#        # Thus, we can simply convolve the LEC output through the gate
#        # gadget and then count up the decoded results of the TEC.
#        
#        # First, propagate the LEC through the Gadget.
#        lecResult = results[self.lecName]
#        ga = self[self.gaName]
#        lecResult.counts, lecResult.keyMeta = ga.propagateCounts(lecResult.counts, lecResult.keyMeta)
#        
#        # Now convolve the Gadget and the LEC.
#        tecResult = results.pop(self.tecName)
#        convolved = super(ExRec, self)._convolve(results, noiseModels, pauli)
#        
#        # Now compute the result of propagating the LEC/CNOT input through the TEC
#        # and then decoding.
#        decodeCounts = convolve(convolved.counts, 
#                                tecResult.counts, 
#                                kMax=self.kGood[pauli], 
#                                convolveFcn=self.convolveTEC,
#                                splitListsInto=[2,1])
#        #decodeCounts = self._convolveTEC(convolved, tecDecoder, self.kGood[pauli])
#        
#        rejectedCounts = convolve(convolved.counts,
#                                  tecResult.rejected,
#                                  kMax=self.kGood[pauli],
#                                  convolveFcn=self.convolveTECRejected,
#                                  splitListsInto=[2,1])
#        
#        
#        return CountResult(decodeCounts, tecResult.blocks, rejectedCounts=rejectedCounts)
#    
#    def _postCount(self, result, noiseModels, pauli):
#        result = super(ExRec, self)._postCount(result, noiseModels, pauli)
#        # TODO remove counts for rectangle correctness??
#        return result
        
    
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
    
#    @staticmethod
#    def convolveTEC(inputCounts, tecCounts):
#        
#        decodeCounts = {}
#        for key, count in inputCounts.iteritems():
#            decoded = tecCounts[key]
#            decoded = {key: val * count for key, val in decoded.iteritems()}
#            decodeCounts = listutils.addDicts(decodeCounts, decoded)
#           
#        return decodeCounts     
#    
#    @staticmethod
#    def convolveTECRejected(inputCounts, tecRejected):
#        # tecRejected is just a single count indexed by [input key][0]
#        rejectKey = key.rejectKey
#        rejectCounts = sum(tecRejected[key][rejectKey] * count for key, count in inputCounts.iteritems())
#        
#        return {rejectKey: rejectCounts}
        
        
#class ExRecBackward(ExRec):
#    pass
#        
#class ExRecForward(ExRec):
#    
#    def _convolve(self, results, noiseModels, pauli):
#        # The LEC should have already stripped all counts of the
#        # logical error information.  That is, all LEC error keys
#        # should now decode to the trivial error.  The LEC
#        # can do this, because, correctness of the rectangle depends
#        # only on the syndrome of its input, not on its logical state.
#        # Thus, we can simply convolve the LEC output through the gate
#        # gadget and then count up the decoded results of the TEC.
#        
#        # First, propagate the LEC through the Gadget.
#        lecResult = results[self.lecName]
#        ga = self[self.gaName]
#        lecResult.counts, lecResult.keyMeta = ga.propagateCounts(lecResult.counts, lecResult.keyMeta)
#        
#        # Now convolve the Gadget and the LEC.
#        tecResult = results.pop(self.tecName)
#        convolved = super(ExRec, self)._convolve(results, noiseModels, pauli)
#        
#        # Now filter out the nonzero syndromes.
##        zeroKey = tuple([0]*convolved.keyMeta.nblocks) # hack!
##        for count in convolved.counts:
##            count.pop(zeroKey, 0)
#        
#        # Now compute the result of propagating the LEC/CNOT input through the TEC
#        decodeCounts = convolve(convolved.counts, 
#                                tecResult.counts, 
#                                kMax=self.kGood[pauli], 
#                                convolveFcn=self.convolveTEC,
#                                splitListsInto=[2,1])
#        #decodeCounts = self._convolveTEC(convolved, tecDecoder, self.kGood[pauli])
#        
#        rejectedCounts = convolve(convolved.counts,
#                                  tecResult.rejected,
#                                  kMax=self.kGood[pauli],
#                                  convolveFcn=self.convolveTECRejected,
#                                  splitListsInto=[2,1])
#        
#        
#        return CountResult(decodeCounts, tecResult.blocks, rejectedCounts=rejectedCounts)
    
#
#class CnotExRec(Component):
#    
#    lecCtrlName = 'LEC (ctrl)'
#    lecTargName = 'LEC (targ)'
#    tecCtrlName = 'TEC (ctrl)'
#    tecTargName = 'TEC (targ)'
#    cnotName = 'cnot'
#    
#    def __init__(self, kGood, kGoodCnot, lecCtrl, lecTarg, tecCtrl, tecTarg):
#        # Construct a transversal CNOT component from the two input codes.
#        ctrlCode = lecCtrl.outBlocks()[0].getCode()
#        targCode = lecTarg.outBlocks()[0].getCode()
#        cnot = TransCnot(kGoodCnot, ctrlCode, targCode)
#
#        subcomponents = {
#                         self.cnotName: cnot,
#                         self.lecCtrlName: lecCtrl,
#                         self.lecTargName: lecTarg,
#                         self.tecCtrlName: tecCtrl,
#                         self.tecTargName: tecTarg
#                        }
#        super(CnotExRec, self).__init__(kGood, subcomponents=subcomponents)
#        
#    def _count(self, noiseModels, pauli):
#        pass
#        
#    def _convolve(self, results, pauli):
#        
#        cnot = self[self.cnotName]
#        
#        lecCtrlResult = results[self.lecCtrlName]
#        lecTargResult = results[self.lecTargName]
#        cnotResult = results[self.cnotName]
#        
#        # First, propagate the input results through the CNOT    
#        lecCtrlResult.counts, lecCtrlResult.keyMeta = cnot.propagateCounts(lecCtrlResult.counts,
#                                                             lecCtrlResult.keyMeta,
#                                                             cnot.ctrlName)
#        lecTargResult.counts, lecTargResult.keyMeta = cnot.propagateCounts(lecTargResult.counts,
#                                                             lecTargResult.keyMeta,
#                                                             cnot.targName)
#        
#        lecCtrlResult.blocks = cnotResult.blocks
#        lecTargResult.blocks = cnotResult.blocks
        
        
class Rectangle(SequentialComponent):
    
#    lecName = 'LEC'
#    gaName = 'Ga'
#    tecName = 'TEC'
    
    def __init__(self, kGood, lec, gadget):
        subs = (lec, gadget)
        super(Rectangle, self).__init__(kGood, subcomponents=subs)