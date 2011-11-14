'''
Created on 2011-09-12

@author: adam
'''

from counting.component import base, bell, teleport
from qec import ed422, qecc, error
from qec.error import Pauli
from settings.noise import NoiseModelXZLowerSympy, NoiseModelXSympy, NoiseModelZSympy, \
    CountingNoiseModelX, CountingNoiseModelZ, NoiseModelXZUpperSympy
import util.cache
from counting.component.adapter import InputAdapter
from counting.component.base import ConcatenatedComponent
from counting.component.transversal import TransCnot
from counting.component.exrec import ExRec, ExRecForward
from counting.component.ec import TECDecodeAdapter, ConcatenatedTEC,\
    LECSyndromeAdapter
from counting import probability
import logging


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
    kGood = {pauli: 3 for pauli in [Pauli.X, Pauli.Z, Pauli.Y]}
    

#    noises = { Pauli.X: NoiseModelXSympy(),
#              Pauli.Z: NoiseModelZSympy(),
#              Pauli.Y: NoiseModelXZLowerSympy() 
#             }
    noises = {Pauli.X: NoiseModelXSympy(),
              Pauli.Z: NoiseModelZSympy(),
              Pauli.Y: NoiseModelXZUpperSympy(),
              }
    
    pauli = Pauli.Y
    code = ed422.ED412Code()

    id = base.Empty(ed422.ED412Code())
    
    ed = makeED(kGood)
    led = InputAdapter(ed, (0,))
    #led = LECSyndromeAdapter(led)
    led = ConcatenatedComponent(kGood, led, led)
    cnot = TransCnot(kGood, code, code)
    ted = TECDecodeAdapter(ed)
    ted = ConcatenatedTEC(kGood, ted, ted)
    
    exRec = ExRecForward(kGood, led, cnot, ted)
    
    p = .005
    gamma = p/15

    prBad = exRec.prBad(noises[pauli], pauli)
    print 'Pr[bad]({0})'.format(p), prBad(gamma)
    
    raise Exception
    
    result = exRec.count(noises, pauli)
    prAccept = exRec.prAccept(noises)
    print 'Pr[accept]({0})'.format(p), prAccept(gamma)
    locTotals = exRec.locations(pauli).getTotals()
    noise = noises[pauli]
    
    print 'result.counts=', result.counts
    for pauli in (Pauli.X, Pauli.Z, Pauli.Y, Pauli.I):
        counts = [{pauli:c.get(pauli, 0)} for c in result.counts]
        countPoly = probability.countsToPoly(counts, locTotals, noise)
        pr = countPoly / prAccept + prBad
        print 'Pr[', pauli, ']({0})='.format(p), pr(gamma)
        print counts

if __name__ == '__main__':
    from counting import countParallel
    countParallel.setPool(countParallel.DummyPool())
    
    
    logger = logging.getLogger('counting.component')
    #logger.setLevel(logging.DEBUG)
    
    logger = logging.getLogger('counting.probability')
    #logger.setLevel(logging.DEBUG)
    
    profile = False
    fetch = True
    util.cache.enableFetch(fetch)
    util.cache.enableMemo(fetch)
        
    if profile:
        import cProfile
        cProfile.run('run()', sort='cumulative')
    else:
        run()
    