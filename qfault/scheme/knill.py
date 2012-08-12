'''
Created on 2011-11-14

Schemes based on Knill's C4/C6 error-detection and teleportation based proposals.

@author: adam
'''
from qfault.counting import probability, countParallel
from qfault.counting.component.adapter import DecodeAdapter, SyndromeAdapter
from qfault.counting.component.base import Prep, ParallelComponent
from qfault.counting.component.bell import BellPair, BellMeas
from qfault.counting.component.exrec import ExRec, Rectangle
from qfault.counting.component.teleport import TeleportED, EDInputFilter, \
    TeleportWithMeas, Teleport
from qfault.counting.component.transversal import TransCnot
from qfault.counting.convolve import convolveDict
from qfault.counting.countParallel import convolve
from qfault.counting.key import SyndromeKeyGenerator, convolveKeyCounts
from qfault.counting.probability import prMinFailures, countsAsProbability
from qfault.counting.result import CountResult
from qfault.qec import ed422, error
from qfault.qec.error import Pauli
from qfault.qec.qecc import StabilizerState, TrivialStablizerCode
from qfault.scheme import Scheme
from qfault.util.polynomial import SymPolyWrapper, sympoly1d
import logging
import operator

logger = logging.getLogger('scheme.knill')

class KnillScheme(Scheme):
    '''
    An error-detection scheme based on [AGP 08], which was in turn based on [Knill 04].
    '''


    def __init__(self, kPrep, kCnot, kEC, kExRec, enableRest=False):
        '''
        Rest locations are disabled by default because we're post-selecting on global acceptance,
        and Pauli frame corrections to teleportation EDs can be made after the gate gadget since
        we're only considering Clifford gates.
        '''

        # Use the same gauge type for each code block  Using |+> gauge qubits gives
        # slightly better Z-error protection.  In principle it is possible to even out
        # the X and Z error protection by using different gauge types, as suggested
        # by AGP and Knill.  However, it does not seem to make a significant difference
        # in the logical error probabilities.
        self.code = ed422.ED412Code(gaugeType=error.xType)
        
        prepZ = Prep(kPrep, ed422.prepare(Pauli.Z, Pauli.X), StabilizerState(self.code, [error.zType]))
        prepX = Prep(kPrep, ed422.prepare(Pauli.X, Pauli.X), StabilizerState(self.code, [error.xType]))
        
        bellPair = BellPair(kEC, prepX, prepZ, kCnot)
        bellMeas = BellMeas(kEC, self.code, kGoodCnot=kCnot, kGoodMeasX=kCnot, kGoodMeasZ=kCnot)
        
        self.bellPair = bellPair
        self.bellMeas = bellMeas
        self.ed = TeleportED(kEC, bellPair, bellMeas, enableRest)
        
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
#        inputs = {0: (0,)}

        cnot = TransCnot(self.kCnot, self.code, self.code)
        led = ParallelComponent(self.kExRec, self.ed, self.ed)
        rec = Rectangle(self.kExRec, led, cnot)
        
        prTable = {}
        for inKeyA in inputs.values():
            for inKeyB in inputs.values():
                logger.info('Counting CNOT 1-Rec for inputs: %s, %s', inKeyA, inKeyB)
                prS = self.prEgivenS(cnot, e, inKeyA + inKeyB, noiseModels) * \
                      self.prSin(inKeyA, noiseModels, Pauli.Y) * self.prSin(inKeyB, noiseModels, Pauli.Y)
                
                prTable[(inKeyA, inKeyB)] = prS
                
        for key, pr in prTable.iteritems():
            print 'Pr[E=({0},{1}), Sin={2}](0.004)={3}'.format(eA,eB, key, pr(0.004/15))
            
        pr = sum(prTable.values())
        logger.debug('Pr[E={0}](0.004)={1}'.format(e, pr(0.004/15)))
                
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
        counts = [{e: count.get(e,0)} for count in result.counts]
        logger.debug('Counts for [E={0} | Sin={1}]: {2}'.format(e, s, counts))
        prE = probability.countsToPoly(counts, exrec.locations().getTotals(), noiseModels[Pauli.Y])
#        logger.debug('Pr[E={0} | Sin={1}] = {2}'.format(e, s, prE))
        logger.debug('Pr[E={0} | Sin={1}](0.004)={2}'.format(e, s, prE(0.004/15)))
        
        return prE
    
    def prSin(self, s, noiseModels, pauli):
        r'''
        Returns an upper bound on :math:`\Pr[S_\text{in} = s]`, the probability that the input syndrome to 
        (a single block of) the rectangle is equivalent to :math:`s`.
        '''
#        # TODO: count all possible input rectangle combinations and compute an upper bound.
#        if (0,) == s:
#            return 1
    
        cnot = TransCnot(self.kCnot, self.code, self.code)
    
        # TODO: this is a crude upper bound of Pr[S_in!=0] <= Pr[K != 0]
        # A better bound can be obtained by counting
#        n = 2 * len(self.ed.locations()) + len(cnot.locations())
#        pr = SymPolyWrapper(sympoly1d([n, 0]))
        
        lec = ParallelComponent(self.kExRec, self.ed, self.ed)
        rec = Rectangle(self.kExRec, lec, cnot)
        rec = SyndromeAdapter(rec)
        
        noise = noiseModels[pauli]
        
        # Construct an upper bound on the input distribution of syndromes to the input rectangle.
        # The bound used here is, Pr[Sin=0] <= 1, Pr[Sin != 0] <= Pr[K >= 1].
        # We can bound Pr[K >= 1] <= B\gamma where B is a sum of the weights of each location in the
        # CNOT rectangle.  The input distribution is then constructed by taking the product
        # (convolution actually) of the counts for one block with itself.
    
        B = sum(noise.getWeight(loc, e) for loc in rec.locations(pauli) for e in noise.errorList(loc))
        nonZeroInputs = self.getInputs(self.ed).values()
        nonZeroInputs.remove((0,))
        inCountsA = [{(0,): 1}, {key: B for key in nonZeroInputs}]
        inCounts = countParallel.convolve(inCountsA, inCountsA, convolveFcn=convolveInCounts)
        print 'inCounts=', inCounts
        inResult = CountResult(inCounts, rec.inBlocks())
        

        # Now count the input rectangle based on the worst case distribution computed above.
        # For the CNOT, the worst case one-block result is computed by taking the maximum
        # count value over each of the two blocks.
        # TODO: need to count all possible input rectangles.        
        result = rec.count(noiseModels, pauli, inputResult=inResult)
        sCounts = []
        for count in result.counts:
            sA = sum(c for key,c in count.iteritems() if s == (key[0],))
            sB = sum(c for key,c in count.iteritems() if s == (key[1],))
            sCounts.append({s: max(sA,sB)})
        
        pr = probability.countsToPoly(sCounts, rec.locations(pauli).getTotals(), noiseModels[pauli])
        pr += rec.prBad(noiseModels[pauli], pauli)

        print 'sCounts=', sCounts
        print 'Pr[s=', s, '](0.004)=', pr(0.004/15)
        return pr
    
def convolveInCounts(counts1, counts2):
    return convolveDict(counts1, counts2, keyOp=operator.add)