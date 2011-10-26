'''
Created on 2011-10-25

@author: adam
'''
from counting.component.base import Component
from counting.component.transversal import TransCnot

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