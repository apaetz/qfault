'''
Created on 2011-09-12

@author: adam
'''

from counting import component
from qec import ed422
from qec.error import Pauli
from settings.noise import NoiseModelXZSympy, NoiseModelXSympy, NoiseModelZSympy

if __name__ == '__main__':
    from counting import countParallel
    countParallel.setPool(countParallel.DummyPool())
    
    prepZ = component.PrepZero(1, 1, ed422.prepare(Pauli.Z, Pauli.X), ed422.ED422ZeroPlus())
    noise = { Pauli.X: NoiseModelXSympy(),
              Pauli.Z: NoiseModelZSympy(),
              Pauli.Y: NoiseModelXZSympy() 
             }
    print prepZ.count(noise).counts()