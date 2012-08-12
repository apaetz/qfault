'''
Created on 2012-04-02

@author: adam
'''

from qfault.counting.block import Block
from qfault.counting.component import transversal
from qfault.counting.component.adapter import IdealDecoder, FlaggingDecoder
from qfault.counting.component.base import ParallelComponent, SequentialComponent, Prep, \
    ConcatenationFilter, Component, Filter
from qfault.counting.component.bell import BellPair, BellMeas
from qfault.counting.component.block import BlockDiscard, BlockPermutation,\
    BlockCombine
from qfault.counting.component.teleport import Teleport, TeleportED, TeleportWithMeas
from qfault.counting.countErrors import error, mapCounts, maxCount
from qfault.counting.key import IdentityManipulator, KeyMerger, SyndromeKeyGenerator,\
    KeyManipulator
from qfault.qec import ed422
from qfault.qec.error import Pauli, PauliError, xType, zType, dualType
from qfault.qec.qecc import StabilizerState, ConcatenatedCode
from qfault.scheme import Scheme
from qfault.util import bits, listutils
import logging
from qfault.util.cache import memoize, fetchable
from qfault.util.plotting import plotList
from qfault.noise import NoiseModelXSympy, NoiseModelZSympy, NoiseModelXZSympy,\
    CountingNoiseModelXZ
from qfault.counting import probability, countErrors
import qfault.counting
from qfault.counting.result import CountResult

logger = logging.getLogger('scheme.fibonacci')

class FibonacciSchemeAP09(object):
    '''
    Fibonacci scheme analyzer based on Aliferis & Preskill 2009.
    '''
    
    def __init__(self, epsilon_scale=1, disable_sbt=[2]):
        self._epsilon_scale = epsilon_scale
        self._disable_sbt = set(disable_sbt)

    @staticmethod
    def Decode(fIn, dfIn, dfBarIn):
        '''
        Decodes a block by transforming "quasi-independent" error strengths in the qubits to logical
        error strengths on the block. See [AP09] page 7.
        
        :param fIn: The probability that a qubit is flagged
        :param dfIn: The probability that a qubit is flagged and contains an error
        :param dfBarIn: The probability that a qubit is not flagged and contains an error
        '''
        
        fOut = 4 * (dfBarIn + fIn**2)
        dfOut = 2 * dfBarIn + 4 * (fIn * dfIn + fIn**2 * (2 * dfBarIn + 3 * dfIn))
        dfBarOut = 4 * (dfBarIn**2 + 2 * fIn * dfBarIn)
        
    #    logger.debug('Decode: {0}->{1}'.format((fIn, dfIn, dfBarIn), (fOut,dfOut, dfBarOut)))
        return fOut, dfOut, dfBarOut
    

    @memoize
    def Epsilon_jj(self, m, j, e):
        '''
        Returns e_m(j,j|j-1), the probability of a level-j error on the output 
        of a j-BP conditioned on acceptance of the input (j-1)-BPs.
        See AP09 Table II.
        '''
        
        f, bF, bFbar = self.BlockTeleport(m, j-1, j, e)
            
        # b^~_m(j,j^~f|j-1)
        _, _, bFbar = self.Decode(f, bF, bFbar)
        
        logger.debug('b_{0}({1},{2},f=0|{3})={4}'.format(m,j,j,j-1,bFbar))
        return bFbar
    
    @memoize
    def Epsilon_j1_j(self, m, j, e):
        '''
        Returns e_m(j-1,j|j-1), the probability of a level-(j-1) error on the output 
        of a j-BP conditioned on acceptance of the input (j-1)-BPs.
        See AP09 Table III.
        '''
        
        if 1 == j:
            # Special case (Eq. 15): e_m(0,1) <= e + C1_m e_m(0,0)
            eps = e * self._epsilon_scale + self.C1(m,0,j,e) * self.Epsilon(m, 0, 0, 0, e)
        elif j in self._disable_sbt:
            # Special case optimization (see Sec. VII A) when sub-block teleportation is disabled
            # for level-j.
            # Here level-(j-1) errors in the j-BP prep propagate just like level-0 errors
            # in the 1-BP prep.  But we can ignore the physical gate errors.
            logger.debug('No sub-block teleportation at level {0}'.format(j))
            eps = self.C1(m,1,j,e) * self.Epsilon(m,j-1,j-1,j-1,e)
        else:
            _, _, sFbar = self.SubblockTeleport(m, j-1, j, e)
            eps = sFbar + self.Epsilon(m, j-1, j-1, j-1, e)
            
        return eps

    @memoize
    def Epsilon(self, m, i, jOut, jIn, e):
        '''
        Returns e_m(i,jOut|jIn), the probability of level-i error on the output of a jOut-BP 
        conditioned on acceptance of either:
        1. The input (j-1)-BPs (jIn = jOut-1)
        2. The jOut-BP (jIn = jOut)
        '''
        
        # Sanity check
        if i > jOut or jIn > jOut:
            raise RuntimeError('Invalid parameters: m={0} i={1} jOut={2} jIn={3} e={4}'.format(m,i,jOut,jIn,e))
        
        if jIn == jOut:
            if i == jOut and 0 == i:
                #e_m(0,0) = e_0
                epsOut = e * self._epsilon_scale
    
            elif jOut in self._disable_sbt and 0 == i:
                # Special case.  There are no sub-block teleportations at this level.
                # Instead, physical qubits undergo 2 rounds of rest noise
                # (and noise from the block teleportation CNOT).
                logger.debug('No sub-block teleportation at level {0}'.format(jOut))
                epsOut = (3*e * self._epsilon_scale + self.C1(m,i,jOut,e) * self.Epsilon(m,i,jOut,jIn-1,e)) / self.PrAccept(jOut,e)
            else:    
                # e_m(i,j|j) = e_m(i,j|j-1) / p(j|j-1)
                epsOut = self.Epsilon(m, i, jOut, jOut-1, e) / self.PrAccept(jOut, e)
                
        elif jIn == jOut - 1:
            if i == jOut:
                # Compute e_m(j,j|j-1) by recursion
                epsOut = self.Epsilon_jj(m, jOut, e)
            elif i == jIn:
                # Compute e_m(j-1,j|j-1) by recursion
                epsOut = self.Epsilon_j1_j(m, jOut, e)
            else:
                # Eq. 23 e_m(i,j|j-1) = e_m(i,j-1|j-1)
                epsOut = self.Epsilon(m, i, jIn, jIn, e)
                            
        logger.debug('e_{0}({1},{2}|{3})={4}'.format(m, i, jOut, jIn, epsOut))
        
        return epsOut

    def EpsilonCSS(self, j, e):
        '''
        Returns e_css(j), the effective noise strength of a level-j CSS
        gadget.
        '''
        
        e_css = 0
        for r in ('c', 't'):
            for m in (xType, zType):
                # r_m(0,j^~f|j) = 3*e + C3_rm * e_m(0,j|j)
                f = rF = 0
                rFbar = 3*e * self._epsilon_scale + self.C3(m,r) * self.Epsilon(m, 0, j, j, e)
                
                for i in range(j):
                    f, rF, rFbar = self.Decode(f, rF, rFbar)
                    propagated = self.C3(m,r) * self.Epsilon(m, i+1, j, j, e)
                    rF += propagated
                    rFbar += propagated
                
                # r_m(j,j|j) = r_m(j,j^f|j) + r_m(j,j^~f|j)
                e_css += rF + rFbar
            
        return e_css

    @memoize
    def C1(self, m, i, j, e):
        '''
        See Eq. 37.
        '''
        
        # For even j, swap the coefficients (see section VII A)
        if not (j % 2):
            m = dualType(m)
        
        if 0 == i:    
            # For level-0, use optimization from Eq. (37)
            c1m = 16 * (self.Epsilon(m, i, i, i, e) + e * self._epsilon_scale)
            if zType == m:
                c1m += 1
        else:
            c1m = {xType: 1, zType: 2}[m]
        
        return c1m
    
    @memoize
    def C2(self, m, j): 
        c2m = {xType: 4, zType: 3}
        
        if not (j % 2):
            m = dualType(m)
        
        return c2m[m]
    
    @memoize
    def C3(self, m, r):
        c3mr = {xType: {'c': 2, 't': 3}, zType: {'c': 3, 't': 3}}
        return c3mr[m][r]
        
    @memoize
    def PrAccept(self, j, e):
        '''
        Returns p(j|j-1), the probability of accepting a j-BP given acceptance of all input (j-1)-BPs.
        '''
        
        
        if 1 == j:
            # p(1) >= (1-e)^N
            N = 72
            p = (1 - (e * self._epsilon_scale))**N
        else:
            # p(j|j-1) >= 1 - sum_m (2 * f_m^b(j,j|j-1) + (8 * f_m^s(j-1,j|j-1))
            
            p = 1
            for m in (xType, zType):
                
                fmb, bmF, bmFbar = self.BlockTeleport(m, j-1, j, e)
                    
                # f^b_m(j,j|j-1)
                fmb, _, _ = self.Decode(fmb, bmF, bmFbar)
                logger.debug('f^b_{0}({1},{2}|{3})={4}'.format(m,j,j,j-1,fmb))
                    
                if j in self._disable_sbt:
                    fms = 0
                else:
                    fms, _, _ = self.SubblockTeleport(m, j-1, j, e)
                    
                p -=  (2 * fmb + 8 * fms)
            
        logger.debug('p({0}|{1})={2}'.format(j,j-1,p))
        return p

    @memoize
    def SubblockTeleport(self, m, i, j, e):
        '''
        Returns {f^s_m(i,j|j-1), s_m(i,j^f|j-1), s_m(i,j^~f|j-1)}, the flagging
        and block error probabilities for level-j subblock teleportation.
        '''
        
        if 0 == i:
            f = sf = 0
            sf_bar = 3*e * self._epsilon_scale + (1 + self.C1(m,0,j,e)) * self.Epsilon(m, 0, j-1, j-1, e) 
        else:
            f, sf, sf_bar = self.SubblockTeleport(m, i-1, j, e)
            
            f, sf, sf_bar = self.Decode(f, sf, sf_bar)
            propagated = (1 + self.C1(m,i,j,e)) * self.Epsilon(m, i, j-1, j-1, e)
            sf += propagated
            sf_bar += propagated    
        
        logger.debug('f^s_{0}({1},{2}|{3})={4}'.format(m,i,j,j-1,f))
        logger.debug('s_{0}({1},{2},f=1|{3})={4}'.format(m,i,j,j-1,sf))    
        logger.debug('s_{0}({1},{2},f=0|{3})={4}'.format(m,i,j,j-1,sf_bar))
        return f, sf, sf_bar


    @memoize
    def BlockTeleport(self, m, i, j, e):
        '''
        Returns {f^b_m(i,j|j-1), b_m(i,j^f|j-1), b_m(i,j^~f|j-1)}, the flagging
        and block error probabilities for level-j block teleportation.
        '''
        
        if 0 == i:
            return 0, 0, 4*e * self._epsilon_scale + self.C2(m,j) * self.Epsilon(m, 0, j-1, j-1, e) 
    
        f, bf, bf_bar = self.BlockTeleport(m, i-1, j, e)
        f, bf, bf_bar = self.Decode(f, bf, bf_bar)
        eps = self.Epsilon(m, i, j-1, j-1, e)
        bf += self.C2(m,j) * eps
        bf_bar += self.C2(m,j) * eps
        
        logger.debug('f^b_{0}({1},{2}|{3})={4}'.format(m,i,j,j-1,f))
        logger.debug('b_{0}({1},{2},f=1|{3})={4}'.format(m,i,j,j-1,bf))    
        logger.debug('b_{0}({1},{2},f=0|{3})={4}'.format(m,i,j,j-1,bf_bar))
        return f, bf, bf_bar


class FibonacciSchemeSyndrome(FibonacciSchemeAP09):
    '''
    Fibonacci scheme analyzer based on Aliferis & Preskill, but enhanced to use
    syndrome counting for level-1 errors.
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self._noise_models = {Pauli.X: NoiseModelXSympy(),
                              Pauli.Z: NoiseModelZSympy(),
                              Pauli.Y: NoiseModelXZSympy(),
                             }
        
        self._code = ed422.ED412Code()
        
        self._kGood = {Pauli.Y: 6}
        
        super(FibonacciSchemeSyndrome, self).__init__(epsilon_scale=8/15.)

        
    def DecodeLevel1Result(self, result, locations):
        '''
        Returns polynomials corresponding to the 
        (flagging, flagged error, unflagged error) strengths resulting
        from decoding the given level-1 syndrome counts.
        '''
        decoder = FlaggingDecoder(self._code)
        result = decoder.count(inputResult=result)
        
        # Count keys are of the form ((flag, error),)
        flag_counts = [0] * len(result.counts)
        flagged_error_counts = [0] * len(result.counts)
        unflagged_error_counts = [0] * len(result.counts)
        
        for k, counts in enumerate(result.counts):
            for key, count in counts.iteritems():
                flag, error = key[0]
                if flag:
                    flag_counts[k] += count
                    if error:
                        flagged_error_counts[k] += count
                elif error:
                    unflagged_error_counts[k] += count
            
        pr_flag = probability.ConstructPolynomial(flag_counts, locations.getTotals(), self._noise_models[Pauli.Y])
        pr_flagged_error = probability.ConstructPolynomial(flagged_error_counts, locations.getTotals(), self._noise_models[Pauli.Y])
        pr_unflagged_error = probability.ConstructPolynomial(unflagged_error_counts, locations.getTotals(), self._noise_models[Pauli.Y])
        
        return pr_flag, pr_flagged_error, pr_unflagged_error


    def BP1(self):
        '''
        Returns a 1-BP component.
        '''     
        bp1 = BP1Max(self._kGood)
        return bp1
    
#    def Epsilon1(self, syndrome_error, j_out, j_in, e):
#        '''
#        Returns epsilon_e(1,j_out,good|j_in) the probability of level-1 syndrome
#        'e' at the output of a (j_out)-BP in the 'good' case, and conditioned on 
#        acceptance of (j_in)-BPs.
#        '''
#        if 1 == j_out:
#            # Base case.  Count the 1-BP.
#            bp_1 = self.BP1()
#            counts = bp_1.count(self._noise_models, Pauli.Y).counts
#            
#            # Take the maximum count for each fault order k.
#            counts = [{syndrome_error: count.get((syndrome_error,), 0)} for count in counts]
#            
#            logger.debug('Counts for [E={0}]: {1}'.format(syndrome_error, counts))
#            epsilon_e = probability.countsToPoly(counts, bp_1.locations().getTotals(), self._noise_models[Pauli.Y])
#            epsilon_e /= self.PrAccept(j_out, e)
#            
#        return epsilon_e(e)
    
    @memoize
    def Epsilon_j1_j(self, m, j, e):
        '''
        Returns e_m(j-1,j|j-1), the probability of a level-(j-1) error on the output 
        of a j-BP conditioned on acceptance of the input (j-1)-BPs.
        See AP09 Table III.
        '''
        
        if 1 == j:
            raise RuntimeError('level-0 not supported')
        elif 2 == j:
            # Special case optimization (see Sec. VII A)
            # Here level-1 errors in the 2-BP prep propagate just like level-0 errors
            # in the 1-BP prep.  But we can ignore the physical gate errors.
            return self.C1(m,1,j,e) * self.Epsilon(m,1,1,1,e)
        
        _, _, sFbar = self.SubblockTeleport(m, j-1, j, e)
        
        return sFbar + self.Epsilon(m, j-1, j-1, j-1, e)
    
    @memoize
    def Epsilon(self, m, i, jOut, jIn, e):
        '''
        Returns e_m(i,jOut|jIn), the probability of level-i error on the output of a jOut-BP 
        conditioned on acceptance of either:
        1. The input (j-1)-BPs (jIn = jOut-1)
        2. The jOut-BP (jIn = jOut)
        '''
        
        if 0 == i and 0 != jOut:
            raise RuntimeError('Level-0 errors are tracked as level-1 syndromes')       
        
        if 1 == i:
            raise RuntimeError('Level-1 errors strengths are not directly available')
                
        return super(FibonacciSchemeSyndrome, self).Epsilon(m, i, jOut, jIn, e)
    
    def EpsilonCSS(self, j, e):
        e_css = 0
        for r in ('c', 't'):
            for m in (xType, zType):
                # Compute the level-1 flagging and error probabilities
                f, rf, rf_bar = [x(e/15.) for x in self._CnotExRecPoly(m, j, r)]
                
                for i in range(1,j):
                    f, rf, rf_bar = self.Decode(f, rf, rf_bar)
                    propagated = self.C3(m,r) * self.Epsilon(m, i+1, j, j, e)
                    rf += propagated
                    rf_bar += propagated
                
                # r_m(j,j|j) = r_m(j,j^f|j) + r_m(j,j^~f|j)
                e_css += rf + rf_bar
        
        return e_css
        
    @memoize
    def PrAccept(self, j, e):
        '''
        Returns p(j|j-1), the probability of accepting a j-BP given acceptance of all input (j-1)-BPs.
        '''
        if 1 == j:
            p = self.BP1().prAccept(self._noise_models)(e/15.)
            logger.debug('p({0}|{1})={2}'.format(j,j-1,p))
        elif 2 == j:
            # This gives us p(j) *without* conditioning on acceptance of (j-1)-BPs.
            # i.e. we have p(2|1) p(1)**3
            p = self._SubblockTeleportComponent(xType, j).prAccept(self._noise_models)(e/15.)
            
            # There are three 1-BPs in the sub-block teleportation.  Divide to
            # get the correct conditional probability.  
            p /= self.PrAccept(1, e)**3
            
            logger.debug('p({0}|{1})={2}'.format(j,j-1,p))
        else:
            p = super(FibonacciSchemeSyndrome, self).PrAccept(j, e)
            
        return p
    
    def _SubblockTeleportComponent(self, m, j):
        if j == 1:
            return self.BP1()
        
#        bp1j = BP1LevelJ(k_good, self.BP1(), j-1)
        
        if xType == m:
            # Select the Z-basis measurements
            tp_output_block = 1
        else:
            # Select X-basis measurements
            tp_output_block = 0             
            
        sbt_j_1 = self._SubblockTeleportComponent(m, j-1)

        k_good = {key: (3*val)*(j) for key, val in self._kGood.iteritems()}        
        sbt = SubblockTeleport(k_good, 
                               ParallelComponent(k_good, sbt_j_1, sbt_j_1), 
                               j % 2,
                               teleport_output_block=tp_output_block)
        return sbt
    
    @memoize
    def _SubblockTeleportCountDecode(self, m, j):
        sbt = self._SubblockTeleportComponent(m, j)
        return self.DecodeLevel1Result(sbt.count(self._noise_models, Pauli.Y), sbt.locations())
    
    @memoize
    def _SubblockTeleportPrBad(self, m, j):
        sbt = self._SubblockTeleportComponent(m, j)
        return sbt.prBad(self._noise_models[Pauli.Y], Pauli.Y)
    
    def _BlockTeleportComponent(self, m, j):
        k_good = {key: val + (j-1) for key, val in self._kGood.iteritems()}
        bp1j = BP1LevelJ(k_good, self.BP1(), j-1)           
        
        if xType == m:
            block_to_meas = (j % 2) + 1
        else:
            block_to_meas = ((j+1) % 2) + 1
        
        bt = BlockTeleport(k_good, ParallelComponent(self._kGood, bp1j, bp1j), m, block_to_meas)
        return bt

    @memoize
    def _BlockTeleportCountDecode(self, m, j):
        bt = self._BlockTeleportComponent(m, j)
        return self.DecodeLevel1Result(bt.count(self._noise_models, Pauli.Y), bt.locations())
    
    @memoize
    def _BlockTeleportPrBad(self, m, j):
        bt = self._BlockTeleportComponent(m, j)
        return bt.prBad(self._noise_models[Pauli.Y], Pauli.Y)

    @memoize
    def _CnotExRecPoly(self, m, j, r):
            output_block = {'c': 0, 't': 1}
            exrec = CnotExRec(self._kGood, m, j, output_block[r])
            return self.DecodeLevel1Result(exrec.count(self._noise_models, Pauli.Y), exrec.locations())
        
    
    @memoize
    def SubblockTeleport(self, m, i, j, e):
        '''
        Returns {f^s_m(i,j|j-1), s_m(i,j^f|j-1), s_m(i,j^~f|j-1)}, the flagging
        and block error probabilities for level-j subblock teleportation.
        '''
        if 0 == i:
            # shouldn't get here
            raise RuntimeError('level-0 strengths are not computed')
        if 1 == j:
            raise RuntimeError('No sub-block teleportation at level-1')
            
        if 1 == i:
            # Level-1 is special because we keep track of syndromes counts
            # on the entire block rather than error probabilities on each
            # qubit.
            # 1. Construct the level-j sub-block teleportation sequence
            # 2. decode the level-1 counting result
            # 3. normalize by the acceptance probabilities
            # 4. Add Pr[bad] (TODO)           
            
            # The sub-block teleportation subcircuit requires three input (j-1)-BPs so
            # normalization is the third power of p(j-1).
            norm = self.PrAccept(j-1, e)**3
            logger.debug('{0}-BP level-1 sub-block teleport normalization = {1}'.format(j, 1/norm))
                    
            f, sf, sf_bar = self._SubblockTeleportCountDecode(m, j)
            pbad = self._SubblockTeleportPrBad(m, j)(e/15.)
            
            f, sf, sf_bar = [x(e/15.) / norm + pbad for x in (f, sf, sf_bar)]
        
            logger.debug('f^s_{0}({1},{2}|{3})={4}'.format(m,i,j,j-1,f))
            logger.debug('s_{0}({1},{2},f=1|{3})={4}'.format(m,i,j,j-1,sf))    
            logger.debug('s_{0}({1},{2},f=0|{3})={4}'.format(m,i,j,j-1,sf_bar))
            
        else:
            f, sf, sf_bar = super(FibonacciSchemeSyndrome, self).SubblockTeleport(m, i, j, e)
            
        if f > 1 or sf > 1 or sf_bar > 1:
            logger.warn('Probability > 1: f^s_{0}({1},{2}|{3}={4}, s_{0}({1},{2},f=1|{3})={5}, s_{0}({1},{2},f=0|{3})={6}'.format(m,i,j,j-1,f,sf,sf_bar))
            
        return f, sf, sf_bar

    @memoize
    def BlockTeleport(self, m, i, j, e):
        '''
        Returns {f^s_m(i,j|j-1), b_m(i,j^f|j-1), b_m(i,j^~f|j-1)}, the flagging
        and block error probabilities for level-j block teleportation.
        '''
        if 0 == i:
            # shouldn't get here
            raise RuntimeError('level-0 strengths are not computed')
            
        if 1 == i:
            # TODO
            # 1. Construct the level-(j-1) sub-block teleportation sequence
            # 1b. Append the level-j block teleportation sequence
            # 2. decode the level-1 counting result
            # 3. normalize by the acceptance probabilities
            # 4. Add Pr[bad]           
            
            # We need to normalize by dividing by the acceptance probability p(j-1)
            # The block teleportation subcircuit requires four input (j-1)-BPs so
            # normalization is the fourth power of p(j-1).
            norm = self.PrAccept(j-1, e)**4
            logger.debug('{0}-BP level-1 block teleport normalization = {1}'.format(j, 1/norm))

            f, bf, bf_bar = self._BlockTeleportCountDecode(m, j)
            pbad = self._BlockTeleportPrBad(m, j)(e/15.)
            
            f, bf, bf_bar = [x(e/15.) / norm + pbad for x in (f, bf, bf_bar)]
 
            logger.debug('f^b_{0}({1},{2}|{3})={4}'.format(m,i,j,j-1,f))
            logger.debug('b_{0}({1},{2},f=1|{3})={4}'.format(m,i,j,j-1,bf))    
            logger.debug('b_{0}({1},{2},f=0|{3})={4}'.format(m,i,j,j-1,bf_bar))
        
        else:
            f, bf, bf_bar = super(FibonacciSchemeSyndrome, self).BlockTeleport(m, i, j, e)
            
        if f > 1 or bf > 1 or bf_bar > 1:
            logger.warn('Probability > 1: f^b_{0}({1},{2}|{3}={4}, b_{0}({1},{2},f=1|{3})={5}, b_{0}({1},{2},f=0|{3})={6}'.format(m,i,j,j-1,f,bf,bf_bar))
        
        return f, bf, bf_bar
    
class SubblockTeleport(SequentialComponent):
    
    def __init__(self, kGood, bell_pair, block_to_teleport=0, teleport_output_block=2, postselect=False):
        code = bell_pair.outBlocks()[0].getCode()
        
        cnot = transversal.TransCnot(kGood, code, code)
        discard = BlockDiscard(cnot.outBlocks(), [block_to_teleport ^ 1])
        bm = BellMeas(kGood, code, kGood, kGood, kGood)
        if postselect:
            if 2 != teleport_output_block:
                raise RuntimeError('teleportation with postselection supports output block 2 only')
            teleport = TeleportED(kGood, bell_pair, bm)
        else:
            teleport = Teleport(kGood, bell_pair, bm, output_block=teleport_output_block)
        
        super(SubblockTeleport, self).__init__(kGood, subcomponents=[bell_pair, cnot, discard, teleport])
        
    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
        r = super(SubblockTeleport, self).count(noiseModels, pauli, inputResult, kMax)
        return r
        
class BlockTeleport(SequentialComponent):
    
    def __init__(self, kGood, bell_pair, meas_basis, block_to_meas):
        code = bell_pair.outBlocks()[0].getCode()
        cnot = transversal.TransCnot({p: 1 for p in kGood.keys()}, code, code)
        
        bpc = SequentialComponent(kGood, subcomponents=[bell_pair, cnot])
        
        discard1 = BlockDiscard(cnot.outBlocks(), [0])
        if 1 == block_to_meas:
            discard2 = BlockDiscard(cnot.outBlocks(), [1])
        else:
            discard2 = BlockDiscard(cnot.outBlocks(), [0])

        bp_cnot1 = SequentialComponent(kGood, [bpc, discard1])
        bp_cnot2 = SequentialComponent(kGood, [bpc, discard2])
        bp_cnot_parallel = ParallelComponent(kGood, bp_cnot1, bp_cnot2)
            
        discard3 = BlockDiscard(cnot.outBlocks(), [block_to_meas - 1])
        
        if xType == meas_basis:
            m = Pauli.X
        else:
            m = Pauli.Z
        meas = transversal.TransMeas(kGood, code, m)
        
        super(BlockTeleport, self).__init__(kGood, subcomponents=[bp_cnot_parallel, cnot, discard3, meas])
        
    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
        r = super(BlockTeleport, self).count(noiseModels, pauli, inputResult, kMax)
        return r
        
        
class CnotGadget(SequentialComponent):
    
    def __init__(self, kGood, j):
        bp1 = BP1Max(kGood)
        bpj = BP1LevelJ(kGood, bp1, j)
        code = bpj.outBlocks()[0].getCode()
        
        bp_in = ParallelComponent(kGood, bpj, bpj)
        cnot = transversal.TransCnot(kGood, code, code)

        super(CnotGadget, self).__init__(kGood, subcomponents=[bp_in, cnot])
        
class CnotExRec(SequentialComponent):
    
    def __init__(self, kGood, m, j, output_block):
        gadget = CnotGadget(kGood, j)
        discard = BlockDiscard(gadget.outBlocks(), [output_block ^ 1])
        
        code = gadget.outBlocks()[output_block].getCode()
        
        bp1 = BP1Max(kGood)
        bpj = BP1LevelJ(kGood, bp1, j)
        
        bell_pair = ParallelComponent(kGood, bpj, bpj)
        bm = BellMeas(kGood, code)
        if xType == m:
            output_block = 1
        else:
            output_block = 0
        teleport = Teleport(kGood, bell_pair, bm, output_block=output_block)
        
        super(CnotExRec, self).__init__(kGood, subcomponents=[gadget, discard, teleport])

        
        
class BP1Max(SequentialComponent):
    '''
    Component that outputs just a single block of the 1-BP.  This block
    contains syndrome counts representing the maximum over both 1-BP blocks. 
    '''
    
    def __init__(self, kGood):
        
        code = ed422.ED412Code(gaugeType=None)
        
        prepZ = Prep(kGood, ed422.prepare(Pauli.Z, Pauli.X), code)
        prepX = Prep(kGood, ed422.prepare(Pauli.X, Pauli.Z), code)
        
        bp = BellPair(kGood, prepX, prepZ, kGood)
        
        bell_meas = BellMeas(kGood, code, kGoodCnot=kGood, kGoodMeasX=kGood, kGoodMeasZ=kGood)
        teleport = TeleportED(kGood, bp, bell_meas)
        parallel_teleport = ParallelComponent(kGood, teleport, teleport)
        combine = BlockCombine(parallel_teleport.outBlocks())
        
        # Temp: remove this!
#        discard = BlockDiscard(bp.outBlocks(), [1])
#        super(BP1Max, self).__init__(kGood, subcomponents=[bp, discard])
#        super(BP1Max, self).__init__(kGood, subcomponents=[prepZ])
        super(BP1Max, self).__init__(kGood, subcomponents=[bp, parallel_teleport, combine])
        
#    def outBlocks(self):
#        return (self[-1].outBlocks()[0],)
#        
#    @fetchable
#    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
#        result = super(BP1Max, self).count(noiseModels, pauli, inputResult, kMax)
#        
#        discard1 = BlockDiscard(self[-1].outBlocks(), [0])
#        discard2 = BlockDiscard(self[-1].outBlocks(), [1])
#        
#        # Separate into two single blocks
#        counts1 = discard2.count(noiseModels, Pauli.Y, result).counts
#        counts2 = discard1.count(noiseModels, Pauli.Y, result).counts
#        
#        # Take the maximum count for each fault order k and each syndrome.
#        counts = countErrors.maxCount(counts1, counts2)
#        
#        return CountResult(counts, self.outBlocks() + result.blocks[2:])

class GaugeFilter(Filter):
    
    def __init__(self, old_code, new_code):
        self._new_code = new_code
        self._old_code = old_code
        self._parity_checks = SyndromeKeyGenerator(self._old_code, '').parityChecks()
        self._gauge_operators = []
        for ops in new_code.gaugeOperators():
            self._gauge_operators += list(ops.values())
        
        super(GaugeFilter, self).__init__()
    
    def inBlocks(self):
        return (Block('', self._old_code),)
    
    def outBlocks(self):
        return (Block('', self._new_code),)
    
    def propagateCounts(self, inputResult):
        result = super(GaugeFilter, self).propagateCounts(inputResult)
        result.blocks[0].code = self._new_code
        return result
    
    def keyPropagator(self, subPropagator=IdentityManipulator()):
        return self.GaugelessPropagator(subPropagator, self._parity_checks, self._gauge_operators)
        
    class GaugelessPropagator(KeyManipulator):
        
        def __init__(self, manipulator, parity_checks, gauge_operators):
            super(GaugeFilter.GaugelessPropagator, self).__init__()
            gauge_operators = set(gauge_operators)
            self._gauge_bits = [i for i in range(len(parity_checks)) if parity_checks[i] in gauge_operators]
            self._n = len(parity_checks)
            
        def _manipulate(self, key):
            key_bit_list = listutils.remove_subsequence(bits.bitsToList(key[0], self._n), self._gauge_bits)
            return (bits.listToBits(key_bit_list),) + key[1:]
        
class BP1LevelJ(SequentialComponent):
    
    def __init__(self, kGood, bp1, j):
        if 2 == j:
            # Level-2 sub-block teleportation is just level-1 teleportation with postselection.
            # So we can count it precisely.
            bp = ParallelComponent(kGood, bp1, bp1)    
            block_to_teleport = j % 2
            sub = SubblockTeleport(kGood, bp, block_to_teleport=block_to_teleport, postselect=True)
        else:
            # For j=1, this is obviously correct.  For j > 2, the sub-block teleportation outputs
            # the bottom half of a level-(j-1) BP.  But logical corrections are based on level-(j-1)
            # (which is greater than one), so there is nothing else to do here.
            sub = bp1                
            
        super(BP1LevelJ, self).__init__(kGood, subcomponents=[sub])
        
        
def PlotPaccept(fibonacci, epsilons, j_max):
    paccept = []
    for j in range(1,j_max+1):
        paccept.append([fibonacci.PrAccept(j, e) for e in epsilons])
        
#    print epsilons
#    print 'p(j|j-1):'
#    for j, pj in enumerate(paccept):
#        print j+1, pj

    plotList(epsilons, paccept, labelList=('1','2','3','4','5'), xLabel=r'$\epsilon$', yLabel=r'$p(j|j-1)$', legendLoc='lower left', filename='fibonacci-paccept',
             xscale='log')

def PlotEpsilonCSS(fibonacci, epsilons, j_max):
    e_css = []
    for j in range(1,j_max+1):
        e_css.append([fibonacci.EpsilonCSS(j, e) for e in epsilons])
        
    plotList(epsilons, e_css, labelList=('1','2','3','4','5'), xLabel=r'$\epsilon$', yLabel=r'$\epsilon_{css}$', legendLoc='lower left', filename='fibonacci-e_css',
             xscale='log', yscale='log')
    
def PlotPrBadSBT(epsilons, j_max):
    fibonacci = FibonacciSchemeSyndrome()
    pbad = []
    for j in range(1,j_max+1):
        print j, len(fibonacci._SubblockTeleportComponent(xType, j).locations(Pauli.Y))
        pbad.append([fibonacci._SubblockTeleportPrBad(xType, j)(e/15.) for e in epsilons])
        
    plotList(epsilons, pbad, labelList=('1','2','3','4','5'), xLabel=r'$\epsilon$', yLabel=r'Pr[bad] (sub-block)', legendLoc='lower left', filename='fibonacci-pbad-sbt',
             xscale='log', yscale='log')   
        
        
def PlotPrBadBP1(epsilons):
    fibonacci = FibonacciSchemeSyndrome()
    pbad = []
    pr_bad = fibonacci.BP1().prBad(fibonacci._noise_models[Pauli.Y], Pauli.Y)
    for j in range(1):
        pbad.append([pr_bad(e/15.) for e in epsilons])
        
    plotList(epsilons, pbad, labelList=('1','2','3','4','5'), xLabel=r'$\epsilon$', yLabel=r'Pr[bad] (1-BP)', legendLoc='lower left', filename='fibonacci-pbad-bp1',
             xscale='log', yscale='log')
        


#############
# Old stuff
#############

class FibonacciScheme(Scheme):
    '''
    classdocs
    '''


    def __init__(self, kBPT):
        '''
        Constructor
        '''
        self.code = ed422.ED412Code(gaugeType=error.xType)
        
        prepZ = Prep(kBPT, ed422.prepare(Pauli.Z, Pauli.X), StabilizerState(self.code, [error.zType]))
        prepX = Prep(kBPT, ed422.prepare(Pauli.X, Pauli.Z), StabilizerState(self.code, [error.xType]))
        
        bellPair = BellPair(kBPT, prepX, prepZ, kBPT)
        bellMeas = BellMeas(kBPT, self.code, kGoodCnot=kBPT, kGoodMeasX=kBPT, kGoodMeasZ=kBPT)
        teleport = TeleportED(kBPT, bellPair, bellMeas)
        
        self.bp = bellPair
        self.bpt = BellPairTeleport(kBPT, bellPair, teleport)
        self.bp2 = BellPair(kBPT, PrepLevel2(kBPT, self.bpt, Pauli.X), PrepLevel2(kBPT, self.bpt, Pauli.Z), kBPT)
        
        print PrepLevel2(kBPT, self.bpt, Pauli.X).inBlocks()
        
        bm2 = BellMeas(kBPT, ConcatenatedCode(self.code, self.code), kGoodCnot=kBPT, kGoodMeasX=kBPT, kGoodMeasZ=kBPT)
        
        teleport2 = TeleportED(kBPT, self.bp2, bm2, acceptFunction=concatenated422Acceptor(ConcatenatedCode(self.code, self.code)))
        self.tp2 = teleport2
        self.bpt2 = BellPairTeleport(kBPT, self.bp2, teleport2)
        
        permutations = ([0,2,1,3], [0,3,1,2], [1,2,0,3], [1,3,0,2])
    
        self.sbtList = [SubblockTeleport(kBPT, self.bpt, perm) for perm in permutations]
        
    def count(self):
        result = self.bpt2.count(self.defaultNoiseModels, Pauli.Y)
        print self.bpt2.outBlocks()
        Component.ValidateResult(result)
        return result

#        results = [sbt.count(self.defaultNoiseModels, Pauli.Y) for sbt in self.sbtList]
#        countss = [result.counts for result in results]
#        print countss
#        countsMax = maxCount(*countss)
#        print countsMax
#        
#        return countsMax
#    
#        bptDecode = BPTLevelOneSingleBlock(self.bpt.kGood, self.bpt)
#        result = bptDecode.count(self.defaultNoiseModels, Pauli.Y)
#        print result.counts

class BellPairTeleport(SequentialComponent):
    
    def __init__(self, kGood, bellPair, teleport):
        parallelTeleport = ParallelComponent(kGood, teleport, teleport)
        super(SequentialComponent, self).__init__(kGood, subcomponents=[bellPair, parallelTeleport])
        

class PrepLevel2(SequentialComponent):
    
    def __init__(self, kGood, bellPair, eigenstate=Pauli.Z):
        twinBP = ParallelComponent(kGood, bellPair, bellPair)
        subs = [twinBP]
        
        if (Pauli.Z == eigenstate):
            # |0> is prepared by permuting the second and third logical blocks.
            perm = BlockPermutation(twinBP.outBlocks(), [0,2,1,3])
            subs.append(perm)
            
        code = subs[-1].outBlocks()[0].getCode()
        subs.append(ConcatenationFilter(code, code))
            
        super(SequentialComponent, self).__init__(kGood, subcomponents=subs)
                
    def count(self, noiseModels, pauli, inputResult=None, kMax=None):
        result = SequentialComponent.count(self, noiseModels, pauli, inputResult=inputResult, kMax=kMax)
        return result
        
class SubblockTeleportOld(SequentialComponent):
    
    def __init__(self, kGood, bellPair, permutation):
        code = bellPair.outBlocks()[0].getCode()
        
        twinBP = ParallelComponent(kGood, bellPair, bellPair)
        perm = BlockPermutation(twinBP.outBlocks(), permutation)
        discard2 = BlockDiscard(perm.outBlocks(), 2)
        cnot = transversal.TransCnot(kGood, code, code)
        discard1 = BlockDiscard(cnot.outBlocks(), 1)
        bm = BellMeas(kGood, code, kGood, kGood, kGood)
        teleport = Teleport(kGood, bellPair, bm)
        
        super(SubblockTeleport, self).__init__(kGood, subcomponents=[twinBP, perm, discard2, cnot, discard1, teleport])    
        
        
        
def concatenated422Acceptor(code):
    botStabilizers = [set(code.stabilizerGeneratorsBottom(q)) for q in range(code.n)] 
    # TODO: more robust way to get parity checks?
    parityChecks = SyndromeKeyGenerator(code, None).parityChecks()
    botSyndromeBits = [[(check in botStabilizers[q]) for check in parityChecks] for q in range(code.n)]
    
    def acceptFcn(syndrome):
        bottom = code.bottom()
        bottomSyndromes = code.bottomSyndromes(syndrome)
        bottomCorrections = [bottom.syndromeCorrection(s) for s in bottomSyndromes]
        bottomFlags = [bottom.detectSyndrome(s) for s in bottomSyndromes]
        
        bottomCorrection = sum(bottomCorrections, PauliError(0))
        topSyndromeCorrection = bits.listToBits(not gen.commutesWith(bottomCorrection) for gen in code.stabilizerGeneratorsTop())
        
        topSyndrome = code.topSyndrome(syndrome) ^ topSyndromeCorrection
                
        accept = True
        if not any(bottomFlags):
            accept = (0 == topSyndrome)
        elif 1 == sum(bottomFlags):
            accept = True
        else:
            types = bottomCorrection.types()
            if (xType in types):
                flagError = PauliError(len(bottomFlags), xbits=bits.listToBits(bottomFlags))
            else:
                flagError = PauliError(len(bottomFlags), zbits=bits.listToBits(bottomFlags))
        
                if any(l == flagError for l in bottom.logicalOperators()[0].values()):
                    accept = False
                    
#        logger.info("syndrome={0} bottomCorrection={1} bottomFlags={5} topSyndromeCorrection={2} topSyndrome={3} accept={4}".format(syndrome, 
#                                                                                                                        bottomCorrection, 
#                                                                                                                        topSyndromeCorrection,
#                                                                                                                        topSyndrome,
#                                                                                                                        accept,
#                                                                                                                        bottomFlags))

        return accept
                    
    return acceptFcn
        
class BPTLevelOneSingleBlock(SequentialComponent):
    
    def __init__(self, kGood, bpt):
        discard = BlockDiscard(bpt.outBlocks(), 1)
        decode = IdealDecoder(discard.outBlocks()[0].getCode())
        
        super(BPTLevelOneSingleBlock, self).__init__(kGood, subcomponents=[bpt, discard, decode])
