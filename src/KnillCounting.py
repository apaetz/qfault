'''
Created on 2011-09-12

@author: adam
'''

from counting.component import base, bell, teleport
from qec import ed422, qecc, error
from qec.error import Pauli
from settings.noise import NoiseModelXZSympy, NoiseModelXSympy, NoiseModelZSympy
import util.cache


def run():
    kMax = {pauli: 2 for pauli in [Pauli.X, Pauli.Z, Pauli.Y]}
    
    prepZ = base.Prep(kMax, ed422.prepare(Pauli.Z, Pauli.X), 
                           #ed422.ED422ZeroPlus())
                           qecc.StabilizerState(ed422.ED422Code(), [error.zType, error.xType]))
    prepX = base.Prep(kMax, ed422.prepare(Pauli.X, Pauli.Z), 
                           qecc.StabilizerState(ed422.ED422Code(), [error.xType, error.zType]))
    noises = { Pauli.X: NoiseModelXSympy(),
              Pauli.Z: NoiseModelZSympy(),
              Pauli.Y: NoiseModelXZSympy() 
             }
    
    pauli = Pauli.X
    noise = noises[pauli]
    
    prepBlock = prepZ.count(noise, pauli)
    print prepBlock.keyMeta
    print prepBlock.counts
    print prepZ.prBad(noise, pauli)
    print prepZ.prBad(noise, pauli)(.001)
    
    bellPair = bell.BellPair(kMax, prepX, prepZ, kMax)
    bellMeas = bell.BellMeas(kMax, ed422.ED422Code(), kMax, kMax, kMax)
    data = base.Empty(ed422.ED422Code()).count(noise, pauli)
    teleportED = teleport.TeleportED(kMax, data, bellPair, bellMeas)
    tedblock = teleportED.count(noise, pauli)
    print tedblock.counts
    print teleportED.prAccept(noises)(.001)

if __name__ == '__main__':
    from counting import countParallel
    countParallel.setPool(countParallel.DummyPool())
    
    profile = False
    fetch = False
    util.cache.enableFetch(fetch)
        
    if profile:
        import cProfile
        cProfile.run('run()', sort='cumulative')
    else:
        run()
    