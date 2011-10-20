'''
Created on 2011-09-12

@author: adam
'''

from counting import component
from qec import ed422, qecc, error
from qec.error import Pauli
from settings.noise import NoiseModelXZSympy, NoiseModelXSympy, NoiseModelZSympy


def run():
    kMax = {pauli: 2 for pauli in [Pauli.X, Pauli.Z, Pauli.Y]}
    prepZ = component.Prep(kMax, ed422.prepare(Pauli.Z, Pauli.X), 
                           #ed422.ED422ZeroPlus())
                           qecc.StabilizerState(ed422.ED422Code(), [error.zType, error.xType]))
    prepX = component.Prep(kMax, ed422.prepare(Pauli.X, Pauli.Z), 
                           qecc.StabilizerState(ed422.ED422Code(), [error.xType, error.zType]))
    noise = { Pauli.X: NoiseModelXSympy(),
              Pauli.Z: NoiseModelZSympy(),
              Pauli.Y: NoiseModelXZSympy() 
             }
    prepBlock = prepZ.count(noise)
    print prepBlock.keyGenerators()
    print prepBlock.counts()
    print prepZ.prBad(noise)
    
    bell = component.BellPair(kMax, prepX, prepZ, kMax)
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
    