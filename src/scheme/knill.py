'''
Created on 2011-11-14

@author: adam
'''
from counting.component.adapter import InputAdapter
from counting.component.base import Prep, Concatenator, FixedOutput
from counting.component.bell import BellPair, BellMeas
from counting.component.ec import TECDecodeAdapter, ConcatenatedTEC, TECAdapter
from counting.component.exrec import ExRecForward, ExRec
from counting.component.teleport import TeleportED, UCTSyndromeOut,\
    TeleportEDFilter
from counting.component.transversal import TransCnot
from qec import ed422, error
from qec.error import Pauli, PauliError
from qec.qecc import StabilizerState
from scheme import Scheme
from counting.key import MultiBlockSyndromeKeyGenerator
import logging
from util.polynomial import SymPolyWrapper, sympoly1d
from counting import probability

logger = logging.getLogger('scheme.knill')

class KnillScheme(Scheme):
    '''
    classdocs
    '''


    def __init__(self, kPrep, kCnot, kEC, kExRec):
        '''
        Constructor
        '''
        
        prepZ = Prep(kPrep, 
                      ed422.prepare(Pauli.Z, Pauli.X), 
                      StabilizerState(ed422.ED412Code(gaugeType=error.xType), [error.zType]))
        prepX = Prep(kPrep, 
                          ed422.prepare(Pauli.X, Pauli.Z), 
                          StabilizerState(ed422.ED412Code(gaugeType=error.xType), [error.xType]))
        
        bellPair = BellPair(kEC, prepX, prepZ, kCnot)
        bellMeas = BellMeas(kEC, ed422.ED412Code(), kGoodCnot=kCnot, kGoodMeasX=kCnot, kGoodMeasZ=kCnot)
        
        ed = TeleportED(kEC, bellPair, bellMeas)
        uct = UCTSyndromeOut(kEC, bellPair, bellMeas)
        
        self.bp = bellPair
        self.bm = bellMeas
        self.ed = ed
        self.uct = uct
        self.kExRec = kExRec
        self.kCnot = kCnot
        self.kEC = kEC
        
    def count(self):
        code = ed422.ED412Code(gaugeType=error.xType)
        keyGen = MultiBlockSyndromeKeyGenerator(self.ed.inBlocks())
        blockname = self.ed.inBlocks()[0].name
        
#        inKeys = set()
#        for eX in range(1 << 4):
#            for eZ in range(1 << 4):
#                pauliError = PauliError(xbits=eX, zbits=eZ)
#                inKeys.add(keyGen.getKey({blockname: pauliError}))

        #inputs = self.getInputs(self.ed)
        inputs = {keyGen.getKey({blockname: Pauli.X}): self.prInputSyndrome(code.getSyndrome(Pauli.X))}

        cnot = TransCnot(self.kCnot, code, code)
        ted = TECDecodeAdapter(self.ed)
        ted = ConcatenatedTEC(self.kExRec, ted, ted)
        
        #inKeys = {keyGen.getKey({blockname: Pauli.X})}
        #exrecs = {'cnot': [ExRecForward(self.kExRec, led, cnot, ted)]}
        
        for inKeyA, inPrA in inputs.iteritems():
            ledA = InputAdapter(self.ed, [{inKeyA: 1}])                
            for inKeyB, inPrB in inputs.iteritems():
                ledB = InputAdapter(self.ed, [{inKeyB: 1}])
            
                led = Concatenator(self.kExRec, ledA, ledB)
                exrec = ExRecForward(self.kExRec, led, cnot, ted)
                logger.info('Counting CNOT exRec for inputs: %s, %s', inKeyA, inKeyB)
                self.countExRec(exrec, self.defaultNoiseModels, prInput=inPrA * inPrB)                
        
#        # Generate CNOT exRecs for all possible configurations of the leading EDs.
#        for haveLEDa, haveLEDb in [
#                                   (False, False),
#                                   (False, True),
#                                   (True, False),
#                                   (True, True)
#                                   ]:
#            if haveLEDa:
#                ledsA = [InputAdapter(self.ed, keyGen.getKey({blockname: Pauli.I}))]
#            else:
#                #ledsA = [FixedOutput(code, [{}, {key:inWeight}], keyGen.keyMeta()) for key in inKeys]
#                ledsA = [self.uct]
#                
#            if haveLEDb:
#                ledsB = [InputAdapter(self.ed, keyGen.getKey({blockname: Pauli.I}))]
#            else:
#                #ledsB = [FixedOutput(code, [{key:inWeight}], keyGen.keyMeta()) for key in inKeys]
#                ledsB = [self.uct]
#                
#            leds = []
#            for ledA in ledsA:
#                for ledB in ledsB:
#                    led = Concatenator(self.kExRec, ledA, ledB)
#                    leds.append(led)
#                    
#            exrecs['cnot'] = exrecs.get('cnot', []) + [ExRecForward(self.kExRec, led, cnot, ted) for led in leds]

    
    def prInputSyndrome(self, s):
        if 0 == s:
            return 1
        # TODO: this is a crude estimate of Pr[S_in!=0].  Need to get an upper bound.
        #return SymPolyWrapper(sympoly1d([(76*15*28*15), 0, 0]))
    
        code = ed422.ED412Code(gaugeType=error.xType)
        led = UCTSyndromeOut(self.kEC, self.bp, self.bm)
        led = Concatenator(self.kExRec, led, led)
        ted = TeleportEDFilter(self.kEC, self.bp, self.bm)
        ted = TECAdapter(ted)
        ted = ConcatenatedTEC(self.kExRec, ted, ted)
        cnot = TransCnot(self.kCnot, code, code)
        exrec = ExRecForward(self.kExRec, led, cnot, ted)
        
        result = exrec.count(self.defaultNoiseModels, Pauli.Y)
        sCounts = [{s:0} for _ in range(len(result.counts))]
        sKey = s # TODO: hack!
        for k,counts in enumerate(result.counts):
            for key, count in counts.iteritems():
                key1, key2 = key # TODO: hack!
                if sKey == key1 or sKey == key2: 
                    sCounts[k][s] += count
            
        
        pr = probability.countsToPoly(sCounts, exrec.locations(Pauli.Y).getTotals(), self.defaultNoiseModels[Pauli.Y])
        # TODO: compute denominator
        pr0 = probability.prMinFailures(0, self.ed.locations(Pauli.Y), self.defaultNoiseModels[Pauli.Y], kMax=0)
        prCnot0 = probability.prMinFailures(0, cnot.locations(Pauli.Y), self.defaultNoiseModels[Pauli.Y], kMax=0)
        pr = pr / (1 - 2/prCnot0 * (1/pr0 - 1))
        return pr    
    
if __name__ == '__main__':
    from counting import countParallel
    countParallel.setPool(countParallel.DummyPool())
    
    logging.getLogger('counting.threshold').setLevel(logging.DEBUG)
    
    kPrep = {Pauli.Y: 2}
    kCnot = {Pauli.Y: 2}
    kEC = {Pauli.Y: 4}
    kExRec = {Pauli.Y: 11}
    scheme = KnillScheme(kPrep, kCnot, kEC, kExRec)
    scheme.count()