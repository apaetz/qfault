'''
Created on 2011-10-25

Rectangle and Extended Rectangle components.

@author: adam
'''
from counting.component.base import SequentialComponent
from util.cache import fetchable

class ExRec(SequentialComponent):
    
    def __init__(self, kGood, lec, gadget, tec):
        subs = (lec, gadget, tec)
        super(ExRec, self).__init__(kGood, subcomponents=subs)
        
    @fetchable
    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
        return super(ExRec, self).count(noiseModels, pauli, inputResult, kMax)                    
        
class Rectangle(SequentialComponent):
    
    def __init__(self, kGood, lec, gadget):
        subs = (lec, gadget)
        super(Rectangle, self).__init__(kGood, subcomponents=subs)
        
    @fetchable
    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
        return super(Rectangle, self).count(noiseModels, pauli, inputResult, kMax)                    
    