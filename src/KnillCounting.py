'''
Created on 2011-09-12

@author: adam
'''

from counting.component import base, bell, teleport
from qec import ed422, qecc, error
from qec.error import Pauli
from settings.noise import NoiseModelXZLowerSympy, NoiseModelXSympy, NoiseModelZSympy, \
    CountingNoiseModelX, CountingNoiseModelZ
import util.cache
from counting.component.base import InputAdapter, ConcatenatedComponent
from counting.component.transversal import TransCnot
from counting.component.exrec import ExRec
from counting.component.ec import TECDecodeAdapter, ConcatenatedTEC
from counting import probability


def makeED(kGood):
    
    prepZ = base.Prep(kGood, 
                      ed422.prepare(Pauli.Z, Pauli.X), 
                      qecc.StabilizerState(ed422.ED412Code(gaugeType=error.xType), [error.zType]))
    prepX = base.Prep(kGood, 
                      ed422.prepare(Pauli.X, Pauli.Z), 
                      qecc.StabilizerState(ed422.ED412Code(gaugeType=error.xType), [error.xType]))
    
    bellPair = bell.BellPair(kGood, prepX, prepZ, kGood)
    bellMeas = bell.BellMeas(kGood, ed422.ED412Code(), kGood, kGood, kGood)
    
    teleportED = teleport.TeleportED(kGood, bellPair, bellMeas)
    
    return teleportED


def run():
    kGood = {pauli: 1 for pauli in [Pauli.X, Pauli.Z, Pauli.Y]}
    

#    noises = { Pauli.X: NoiseModelXSympy(),
#              Pauli.Z: NoiseModelZSympy(),
#              Pauli.Y: NoiseModelXZLowerSympy() 
#             }
    noises = {Pauli.X: NoiseModelXSympy(),
              Pauli.Z: NoiseModelZSympy(),
              Pauli.Y: NoiseModelXZLowerSympy(),
              }
    
    pauli = Pauli.Y
    code = ed422.ED412Code()

    id = base.Empty(ed422.ED412Code())
    
    ed = makeED(kGood)
    led = InputAdapter(ed, (0,))
    led2 = ConcatenatedComponent(kGood, led, led)
    cnot = TransCnot(kGood, code, code)
    ted = TECDecodeAdapter(ed)
    ted2 = ConcatenatedTEC(kGood, ted, ted)
    
    exRec = ExRec(kGood, led, id, ted)
    
    result = exRec.count(noises, pauli)
    countPoly = probability.countsToPoly(result.counts, exRec.locations(pauli).getTotals(), noises[pauli])
    prBad = exRec.prBad(noises[pauli], pauli)
    prAccept = exRec.prAccept(noises)
    pr = countPoly / prAccept + prBad
    print pr(0.001)

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
    