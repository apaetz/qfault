'''
Created on 2011-11-14

@author: adam
'''
from counting import probability
from counting.component.adapter import InputAdapter
from counting.component.base import Prep, ParallelComponent, Empty
from counting.component.bell import BellPair, BellMeas
from counting.component.ec import DecodeAdapter, ConcatenatedTEC, TECAdapter
from counting.component.exrec import ExRec, Rectangle
from counting.component.teleport import TeleportED, \
    TeleportEDFilter, Teleport, EDInputFilter
from counting.component.transversal import TransCnot
from counting.key import MultiBlockSyndromeKeyGenerator
from counting.result import CountResult
from qec import ed422, error
from qec.error import Pauli, PauliError
from qec.qecc import StabilizerState
from scheme import Scheme
from util.polynomial import SymPolyWrapper, sympoly1d
import logging

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
#        uct = UCTSyndromeOut(kEC, bellPair, bellMeas)
        
        self.bp = bellPair
        self.bm = bellMeas
        self.ed = ed
#        self.uct = uct
        self.kExRec = kExRec
        self.kCnot = kCnot
        self.kEC = kEC
        self.code = ed422.ED412Code(gaugeType=error.xType)
        
    def count(self):
        # TODO: prEgivenAccept should take Pauli types as input, rather than syndrome keys.
        pr = self.prEgivenAccept((1,0), self.defaultNoiseModels)
        print pr(0.007/15)


    def prEgivenAccept(self, e, noiseModels):
        r'''
        Returns an upper bound on :math:`\Pr[E=e \vert \text{global accept}`, the probability that a rectangle
        induces logical error :math:`e` given acceptance of all EDs.
        '''
        
        code = ed422.ED412Code(gaugeType=error.xType)

        inputs = self.getInputs(self.ed)
        #inputs = {keyGen.getKey({blockname: Pauli.X}): self.prInputSyndrome(code.getSyndrome(Pauli.X))}

        cnot = TransCnot(self.kCnot, code, code)
        led = ParallelComponent(self.kExRec, self.ed, self.ed)
        rec = Rectangle(self.kExRec, led, cnot)
        
        #inKeys = {keyGen.getKey({blockname: Pauli.X})}
        #exrecs = {'cnot': [ExRecForward(self.kExRec, led, cnot, ted)]}
        
        prTable = {}
        for inKeyA in inputs.values():           
            for inKeyB in inputs.values():
                logger.info('Counting CNOT 1-Rec for inputs: %s, %s', inKeyA, inKeyB)
                prS = self.prEgivenS(cnot, e, inKeyA + inKeyB, noiseModels) * \
                      self.prInputSyndrome(inKeyA) * self.prInputSyndrome(inKeyB)
                
                prTable[(inKeyA, inKeyB)] = prS
                
        for key, pr in prTable.iteritems():
            print key, pr(0.007/15)
            
        pr = sum(prTable.values())
                
        # Compute Pr[E=e, ED_R accept, good | Sin=s]
        # TODO: big-time kludges here
        led = ParallelComponent(self.kExRec, self.ed, self.ed)
        
        prBad = rec.prBad(noiseModels[Pauli.Y], Pauli.Y)
        rectLocs = cnot.locations() + self.ed.locations() + self.ed.locations()
        prK0 = probability.prMinFailures(0, rectLocs, noiseModels[Pauli.Y], 0)
        omega = self.omega(noiseModels[Pauli.Y], cnot)
        
        print 'Pr[bad](0.007)=', prBad(0.007/15)
        print 'pr[K=0](0.007)=', prK0(0.007/15)
        print 'omega(0.007)=', omega(0.007/15)
        
        pr = (pr + prBad) / (prK0 * (1 - omega))
        
        return pr
    
    def omega(self, noiseModel, ga):
        Lambda = SymPolyWrapper(sympoly1d([7,0])) + probability.prMinFailures(2, self.ed.locations(), noiseModel)
        e = 2.72
        prK0Ga = probability.prMinFailures(0, ga.locations(), noiseModel, kMax=0)
        prK0ED = probability.prMinFailures(0, self.ed.locations(), noiseModel, kMax=0)
        
        # (6e) 4 * Lambda / (Pr[K_Ga = 0] Pr[K_ED = 0])
        return 6*e*4*Lambda / (prK0Ga * prK0ED)
        
        
    
    def prEgivenS(self, ga, e, s, noiseModels):
        r'''
        Returns an upper bound on :math:`\Pr[E=e, Ex_R~\text{accept}, \text{good} \vert S_\text{in}=s]` for
        input syndrome :math:`s` and logical error :math:`e`.
        '''
        
        # Compute Pr[E=e, ED_R accept, good | Sin=s]
        eds = [self.ed for _ in s]
        led = ParallelComponent(self.kExRec, *eds)
        rec = Rectangle(self.kExRec, led, ga)
        
        # We're really counting the rectangle. So we want the output of the Ga conditioned on
        # acceptance of the TED, and not the output of the TED itself.
        teds = [DecodeAdapter(EDInputFilter(ed)) for ed in eds]
        ted = ParallelComponent(self.kExRec, *teds) 
        exrec = ExRec(self.kExRec, led, ga, ted)
        
        # TODO: hack!
        inputResult = CountResult([{s: 1}], rec.inBlocks())

        result = exrec.count(noiseModels, Pauli.Y, inputResult=inputResult)
        counts = [{e: count.get(e,0)} for count in result.counts]
        prE = probability.countsToPoly(counts, exrec.locations().getTotals(), noiseModels[Pauli.Y])
        
        return prE
    
    def prInputSyndrome(self, s):
        r'''
        Returns an upper bound on :math:`\Pr[S_\text{in} = s]`, the probability that the input syndrome to 
        (a single block of) the rectangle is equivalent to :math:`s`.
        '''
        # TODO: count all possible input rectangle combinations and compute an upper bound.
        if (0,) == s:
            return 1
    
        code = ed422.ED412Code(gaugeType=error.xType)
        led = Teleport(self.kEC, self.bp, self.bm)
#        led = ParallelComponent(self.kExRec, led, led)
#        ted = TeleportEDFilter(self.kEC, self.bp, self.bm)
#        ted = TECAdapter(ted)
#        ted = ConcatenatedTEC(self.kExRec, ted, ted)
        cnot = TransCnot(self.kCnot, code, code)
#        exrec = ExRecForward(self.kExRec, led, cnot, ted)
        
        # TODO: crude upper bound of Pr[S_in!=0] <= Pr[K != 0]
        # A better bound can be obtained by counting
        n = 2 * len(led.locations()) + len(cnot.locations())
        pr = SymPolyWrapper(sympoly1d([n, 0]))
#        
#        result = exrec.count(self.defaultNoiseModels, Pauli.Y)
#        sCounts = [{s:0} for _ in range(len(result.counts))]
#        sKey = s # TODO: hack!
#        for k,counts in enumerate(result.counts):
#            for key, count in counts.iteritems():
#                key1, key2 = key # TODO: hack!
#                if sKey == key1 or sKey == key2: 
#                    sCounts[k][s] += count
#            
#        
#        pr = probability.countsToPoly(sCounts, exrec.locations(Pauli.Y).getTotals(), self.defaultNoiseModels[Pauli.Y])
#        # TODO: compute denominator
#        pr0 = probability.prMinFailures(0, self.ed.locations(Pauli.Y), self.defaultNoiseModels[Pauli.Y], kMax=0)
#        prCnot0 = probability.prMinFailures(0, cnot.locations(Pauli.Y), self.defaultNoiseModels[Pauli.Y], kMax=0)
#        pr = pr / (1 - 2/prCnot0 * (1/pr0 - 1))

        print 'Pr[s=', s, '](0.007)=', pr(0.007/15)
        return pr    
    
if __name__ == '__main__':
    from counting import countParallel
    countParallel.setPool(countParallel.DummyPool())
    
    logging.getLogger('counting.threshold').setLevel(logging.DEBUG)
    
    kPrep = {Pauli.Y: 3}
    kCnot = {Pauli.Y: 3}
    kEC = {Pauli.Y: 4}
    kExRec = {Pauli.Y: 6}
    scheme = KnillScheme(kPrep, kCnot, kEC, kExRec)
    scheme.count()