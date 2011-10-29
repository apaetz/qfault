'''
Created on 2011-09-12

@author: adam
'''

from counting.component import base, bell, teleport
from qec import ed422, qecc, error
from qec.error import Pauli
from settings.noise import NoiseModelXZSympy, NoiseModelXSympy, NoiseModelZSympy, \
    CountingNoiseModelX, CountingNoiseModelZ
import util.cache
from counting.component.base import InputAdapter, ConcatenatedComponent
from counting.component.transversal import TransCnot
from counting.component.exrec import ExRec
from counting.component.ec import TECDecodeAdapter, ConcatenatedTEC


def makeED(kGood):
    
    prepZ = base.Prep(kGood, 
                      ed422.prepare(Pauli.Z, Pauli.X), 
                      qecc.StabilizerState(ed422.ED422Code(), [error.zType, error.xType]))
    prepX = base.Prep(kGood, 
                      ed422.prepare(Pauli.X, Pauli.Z), 
                      qecc.StabilizerState(ed422.ED422Code(), [error.xType, error.zType]))
    
    bellPair = bell.BellPair(kGood, prepX, prepZ, kGood)
    bellMeas = bell.BellMeas(kGood, ed422.ED422Code(), kGood, kGood, kGood)
    
    teleportED = teleport.TeleportED(kGood, bellPair, bellMeas)
    
    return teleportED


def run():
    kGood = {pauli: 1 for pauli in [Pauli.X, Pauli.Z, Pauli.Y]}
    

#    noises = { Pauli.X: NoiseModelXSympy(),
#              Pauli.Z: NoiseModelZSympy(),
#              Pauli.Y: NoiseModelXZSympy() 
#             }
    noises = {Pauli.X: CountingNoiseModelX(),
              Pauli.Z: CountingNoiseModelZ(),
              Pauli.Y: None,
              }
    
    pauli = Pauli.X
    code = ed422.ED422Code()

    data = base.Empty(ed422.ED422Code()).count(noises, pauli)
    
    ed = makeED(kGood)
    led = InputAdapter(ed, (0,))
    led = ConcatenatedComponent(kGood, led, led)
    cnot = TransCnot(kGood, code, code)
    ted = TECDecodeAdapter(ed)
    ted = ConcatenatedTEC(kGood, ted, ted)
    
    exRec = ExRec(kGood, led, cnot, ted)
    
    result = exRec.count(noises, pauli)
    print result.counts

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
    