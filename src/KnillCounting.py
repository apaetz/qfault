'''
Created on 2011-09-12

@author: adam
'''

from counting import component
from qec import ed422, qecc, error
from qec.error import Pauli
from settings.noise import NoiseModelXZSympy, NoiseModelXSympy, NoiseModelZSympy

if __name__ == '__main__':
    from counting import countParallel
    countParallel.setPool(countParallel.DummyPool())
    
    prepZ = component.Prep(1, 1, ed422.prepare(Pauli.Z, Pauli.X), 
                           #ed422.ED422ZeroPlus())
                           qecc.StabilizerState(ed422.ED422Code(), [error.zType, error.xType]))
    noise = { Pauli.X: NoiseModelXSympy(),
              Pauli.Z: NoiseModelZSympy(),
              Pauli.Y: NoiseModelXZSympy() 
             }
    print prepZ.count(noise).counts()
    
    #tcnot = TransCNOT(2, 2, ed422.ED422Code(), ed422.ED422Code())
    #print tcnot.count(noise).counts()
    
    #TODO next: construct a 1-ED gadget.