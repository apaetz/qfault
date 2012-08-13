'''
Created on 2012-04-02

@author: adam
'''

from fibonacci_components import *
from qfault.counting import probability
from qfault.counting.component.adapter import FlaggingDecoder
from qfault.noise.noise import NoiseModelXSympy, NoiseModelZSympy, \
    NoiseModelXZSympy
from qfault.qec.error import zType, dualType
from qfault.util.cache import memoize

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