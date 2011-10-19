'''
Created on 2011-09-12

@author: adam
'''

from counting import component
from qec import ed422, qecc, error
from qec.error import Pauli
from settings.noise import NoiseModelXZSympy, NoiseModelXSympy, NoiseModelZSympy


def run():
    prepZ = component.Prep(3, 3, ed422.prepare(Pauli.Z, Pauli.X), 
                           #ed422.ED422ZeroPlus())
                           qecc.StabilizerState(ed422.ED422Code(), [error.zType, error.xType]))
    prepX = component.Prep(3, 3, ed422.prepare(Pauli.X, Pauli.Z), 
                           qecc.StabilizerState(ed422.ED422Code(), [error.xType, error.zType]))
    noise = { Pauli.X: NoiseModelXSympy(),
              Pauli.Z: NoiseModelZSympy(),
              Pauli.Y: NoiseModelXZSympy() 
             }
    prepBlock = prepZ.count(noise)
    print prepBlock.keyGenerators()
    print prepBlock.counts()
    
    bell = component.Bell(3, 3, prepX, prepZ)
    bellBlock = bell.count(noise)
    print bellBlock.keyGenerators()
    print bellBlock.counts()
    
    #TODO next: construct a 1-ED gadget.


if __name__ == '__main__':
    from counting import countParallel
    countParallel.setPool(countParallel.DummyPool())
    
    profile = False
        
    if profile:
        import cProfile
        cProfile.run('run()', sort='cumulative')
    else:
        run()
    