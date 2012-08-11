'''
Created on 2011-10-25

@author: adam
'''
from counting.component.base import Component
import logging

logger = logging.getLogger('component')

class VerifyX(Component):
    
    cmName = 'cnotMeas'
    
    def __init__(self, kGood, kBest, zeroPrepA, zeroPrepB):
        cnotMeas = 0 # TODO
        subcomponents = [zeroPrepA, zeroPrepB, cnotMeas]
        super(VerifyX, self).__init__(kGood, kBest, subcomponents=subcomponents)
        self.prepAName = zeroPrepA.nickname()
        self.prepBName = zeroPrepB.nickname()
        
#    def _convolve(self, counts):
#        results = {}
#        
#        countsA = counts[self.prepAName]
#        countsB = counts[self.prepBName]
#        countsCM = counts[self.cmName]
#        
#        # First, X-error counts.
#        countsBCM = convolve(countsCM[Pauli.X], 
#                             countsB[Pauli.X], 
#                             convolveFcn=convolveABB, 
#                             kGood=self.kGood)
#        
#        countsVerified = convolve(countsA[Pauli.X], 
#                                  countsBCM, 
#                                  convolveFcn=convolveCountsPostselectX, 
#                                  kGood=self.kGood)
#        
#        results[Pauli.X] = countsVerified
#        
#
#        # Now Z errors.
#        
#        # TODO.  Need to reduce CM counts over two blocks to
#        # marginal counts over just block A.  This will require
#        # constructing the 2-block counts a bit differently.
#        # i.e., 2-block counts are indexed by [s1][s2] rather
#        # than a single index for the entire syndrome.
#        countsVerify = convolve(countsPrepB, countsC, kGood=kGood)
#        countsVerified = convolve(countsPrepA, countsVerify, kGood=kGood)
#    
#        
#        
#        return CountResult(str(countsA), countsA.getCode, results)