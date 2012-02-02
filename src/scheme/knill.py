'''
Created on 2011-11-14

Schemes based on Knill's C4/C6 error-detection and teleportation based proposals.

@author: adam
'''
from counting import probability
from counting.component.adapter import DecodeAdapter
from counting.component.base import Prep, ParallelComponent
from counting.component.bell import BellPair, BellMeas
from counting.component.exrec import ExRec, Rectangle
from counting.component.teleport import TeleportED, EDInputFilter
from counting.component.transversal import TransCnot
from counting.result import CountResult
from qec import ed422, error
from qec.error import Pauli
from qec.qecc import StabilizerState, TrivialStablizerCode
from scheme import Scheme
from util.polynomial import SymPolyWrapper, sympoly1d
import logging
from counting.key import SyndromeKeyGenerator

logger = logging.getLogger('scheme.knill')

class KnillScheme(Scheme):
    '''
    An error-detection scheme based on [AGP 08], which was in turn based on [Knill 04].
    '''


    def __init__(self, kPrep, kCnot, kEC, kExRec):
        '''
        Constructor
        '''

        # Use the same gauge type for each code block  Using |+> gauge qubits gives
        # slightly better Z-error protection.  In principle it is possible to even out
        # the X and Z error protection by using different gauge types, as suggested
        # by AGP and Knill.  However, it does not seem to make a significant difference
        # in the logical error probabilities.
        self.code = ed422.ED412Code(gaugeType=error.xType)
        
        prepZ = Prep(kPrep, ed422.prepare(Pauli.Z, Pauli.X), StabilizerState(self.code, [error.zType]))
        prepX = Prep(kPrep, ed422.prepare(Pauli.X, Pauli.Z), StabilizerState(self.code, [error.xType]))
        
        bellPair = BellPair(kEC, prepX, prepZ, kCnot)
        bellMeas = BellMeas(kEC, self.code, kGoodCnot=kCnot, kGoodMeasX=kCnot, kGoodMeasZ=kCnot)
        
        self.ed = TeleportED(kEC, bellPair, bellMeas)
        
        self.kExRec = kExRec
        self.kCnot = kCnot
        self.kEC = kEC
        
    def count(self):
        # TODO: prEgivenAccept should take Pauli types as input, rather than syndrome keys.
        
        for eA in [Pauli.I, Pauli.X, Pauli.Z, Pauli.Y]:
            for eB in [Pauli.I, Pauli.X, Pauli.Z, Pauli.Y]:
                pr = self.prEgivenAccept(eA, eB, self.defaultNoiseModels)
                print 'Pr[{0}{1}](0.004)={2}'.format(eA,eB,pr(0.004/15))


    def prEgivenAccept(self, eA, eB, noiseModels):
        r'''
        Returns an upper bound on :math:`\Pr[E=e \vert \text{global accept}]`, the probability that a rectangle
        induces logical error :math:`e` given acceptance of all EDs.
        '''
        
        generator = SyndromeKeyGenerator(TrivialStablizerCode(), None)
        e = (generator.getKey(eA), generator.getKey(eB))

        inputs = self.getInputs(self.ed)

        cnot = TransCnot(self.kCnot, self.code, self.code)
        led = ParallelComponent(self.kExRec, self.ed, self.ed)
        rec = Rectangle(self.kExRec, led, cnot)
        
        prTable = {}
        for inKeyA in inputs.values():
            for inKeyB in inputs.values():
                logger.info('Counting CNOT 1-Rec for inputs: %s, %s', inKeyA, inKeyB)
                prS = self.prEgivenS(cnot, e, inKeyA + inKeyB, noiseModels) * \
                      self.prSin(inKeyA) * self.prSin(inKeyB)
                
                prTable[(inKeyA, inKeyB)] = prS
                
        for key, pr in prTable.iteritems():
            print key, pr(0.004/15)
            
        pr = sum(prTable.values())
                
        # Compute Pr[E=e, ED_R accept, good | Sin=s]      
        prBad = rec.prBad(noiseModels[Pauli.Y], Pauli.Y)
        prK0 = probability.prMinFailures(0, rec.locations(), noiseModels[Pauli.Y], 0)
        omega = self.omega(noiseModels[Pauli.Y], cnot)
        
        print 'Pr[bad](0.004)=', prBad(0.004/15)
        print 'pr[K=0](0.004)=', prK0(0.004/15)
        print 'omega(0.004)=', omega(0.004/15)
        
        pr = (pr + prBad) / (prK0 * (1 - omega))
        
        return pr
    
    def omega(self, noiseModel, ga):
        prefactor = SymPolyWrapper(sympoly1d([-15,1]))**7 * SymPolyWrapper(sympoly1d([-4,1]))**15 * SymPolyWrapper(sympoly1d([-8,1]))
        Lambda = 4*prefactor*SymPolyWrapper(sympoly1d([7,0])) + probability.prMinFailures(2, self.ed.locations(), noiseModel)
        e = 2.72
        prK0Ga = probability.prMinFailures(0, ga.locations(), noiseModel, kMax=0)
        prK0ED = probability.prMinFailures(0, self.ed.locations(), noiseModel, kMax=0)
        
        # (6e) 4 * Lambda / (Pr[K_Ga = 0] Pr[K_ED = 0])
        return 6*e*4*Lambda / (prK0Ga * prK0ED)
        
        
    
    def prEgivenS(self, ga, e, s, noiseModels):
        r'''
        Returns an upper bound on :math:`\Pr[E=e, \text{Ex}_R~\text{accept}, \text{good} \vert S_\text{in}=s]` for
        input syndrome :math:`s` and logical error :math:`e`.
        '''
        
        # Compute Pr[E=e, ED_R accept, good | Sin=s]
        eds = [self.ed for _ in s]
        led = ParallelComponent(self.kExRec, *eds)
        
        # We're really counting the rectangle. So we want the output of the Ga conditioned on
        # acceptance of the TED, and not the output of the TED itself.
        teds = [DecodeAdapter(EDInputFilter(ed)) for ed in eds]
        ted = ParallelComponent(self.kExRec, *teds) 
        exrec = ExRec(self.kExRec, led, ga, ted)
        
        # TODO: hack!
        inputResult = CountResult([{s: 1}], exrec.inBlocks())

        logger.debug('Counting exRec for Pr[E={0} | Sin={1}]'.format(e, s))
        result = exrec.count(noiseModels, Pauli.Y, inputResult=inputResult)
        logger.debug('Counts before filtering for E={0}: {1}'.format(e, result.counts))
        counts = [{e: count.get(e,0)} for count in result.counts]
        logger.debug('Counts after filtering for E={0}: {1}'.format(e, counts))
        prE = probability.countsToPoly(counts, exrec.locations().getTotals(), noiseModels[Pauli.Y])
        
        return prE
    
    def prSin(self, s):
        r'''
        Returns an upper bound on :math:`\Pr[S_\text{in} = s]`, the probability that the input syndrome to 
        (a single block of) the rectangle is equivalent to :math:`s`.
        '''
        # TODO: count all possible input rectangle combinations and compute an upper bound.
        if (0,) == s:
            return 1
    
        cnot = TransCnot(self.kCnot, self.code, self.code)
    
        # TODO: this is a crude upper bound of Pr[S_in!=0] <= Pr[K != 0]
        # A better bound can be obtained by counting
        n = 2 * len(self.ed.locations()) + len(cnot.locations())
        pr = SymPolyWrapper(sympoly1d([n, 0]))

        print 'Pr[s=', s, '](0.004)=', pr(0.004/15)
        return pr
    
if __name__ == '__main__':
    from counting import countParallel
    countParallel.setPool(countParallel.DummyPool())
    
    logging.getLogger('counting.threshold').setLevel(logging.DEBUG)
    logging.getLogger('util.cache').setLevel(logging.DEBUG)
    logging.getLogger('scheme').setLevel(logging.DEBUG)
    logger.setLevel(logging.DEBUG)
    
    kPrep = {Pauli.Y: 3}
    kCnot = {Pauli.Y: 3}
    kEC = {Pauli.Y: 4}
    kExRec = {Pauli.Y: 6}
    scheme = KnillScheme(kPrep, kCnot, kEC, kExRec)
    scheme.count()