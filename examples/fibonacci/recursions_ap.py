'''
Created on 2012-04-16

Recursion relations for the Fibonacci scheme as given by Aliferis and Preskill
arXiv:0809.5063.

@author: adam
'''
from qfault.util.cache import memoize
import logging
from qfault.util.plotting import plotList
import math

logger = logging.getLogger('fibonacci')

xType = 'X'
zType = 'Z'
_dualTable = {xType: zType,
             zType: xType}

def dualType(etype):
    '''
    Returns the dual of the given error type.
    >>> dualType(xType)
    'Z'
    >>> dualType(zType)
    'X'
    '''
    return _dualTable[etype]

def decode(fIn, dfIn, dfBarIn):
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
    
#    logger.debug('decode: {0}->{1}'.format((fIn, dfIn, dfBarIn), (fOut,dfOut, dfBarOut)))
    return fOut, dfOut, dfBarOut

@memoize
def recursionAP_jj(m, j, e):
    '''
    Returns e_m(j,j|j-1), the probability of a level-j error on the output 
    of a j-BP conditioned on acceptance of the input (j-1)-BPs.
    See AP09 Table II.
    '''
    
    f, bF, bFbar = blockTeleport(m, j-1, j, e)
        
    # b^~_m(j,j^~f|j-1)
    _, _, bFbar = decode(f, bF, bFbar)
    
    logger.debug('b_{0}({1},{2},f=0|{3})={4}'.format(m,j,j,j-1,bFbar))
    return bFbar
    
@memoize
def recursionAP_j1j(m, j, e):
    '''
    Returns e_m(j-1,j|j-1), the probability of a level-(j-1) error on the output 
    of a j-BP conditioned on acceptance of the input (j-1)-BPs.
    See AP09 Table III.
    '''
    
    if 1 == j:
        # Special case (Eq. 15): e_m(0,1) <= e + c1_m e_m(0,0)
        return e + c1(m,0,j,e) * epsilon(m, 0, 0, 0, e)
    elif 2 == j:
        # Special case optimization (see Sec. VII A)
        # Here level-1 errors in the 2-BP prep propagate just like level-0 errors
        # in the 1-BP prep.  But we can ignore the physical gate errors.
        return c1(m,1,j,e) * epsilon(m,1,1,1,e)
    
    _, _, sFbar = subblockTeleport(m, j-1, j, e)
    
    return sFbar + epsilon(m, j-1, j-1, j-1, e)

@memoize
def epsilon(m, i, jOut, jIn, e):
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
            epsOut = e

        elif 2 == jOut and 0 == i:
            # Special case.  There are no sub-block teleportations
            # at level j=2.  Instead, physical qubits undergo 2 rounds of rest noise
            # (and noise from the block teleportation CNOT).
            epsOut = (3*e + c1(m,i,jOut,e) * epsilon(m,i,jOut,jIn-1,e)) / p_accept(jOut,e)
        else:    
            # e_m(i,j|j) = e_m(i,j|j-1) / p(j|j-1)
            epsOut = epsilon(m, i, jOut, jOut-1, e) / p_accept(jOut, e)
            
    elif jIn == jOut - 1:
        if i == jOut:
            # Compute e_m(j,j|j-1) by recursion
            epsOut = recursionAP_jj(m, jOut, e)
        elif i == jIn:
            # Compute e_m(j-1,j|j-1) by recursion
            epsOut = recursionAP_j1j(m, jOut, e)
        else:
            # Eq. 23 e_m(i,j|j-1) = e_m(i,j-1|j-1)
            epsOut = epsilon(m, i, jIn, jIn, e)
                        
    logger.debug('e_{0}({1},{2}|{3})={4}'.format(m, i, jOut, jIn, epsOut))
    
    return epsOut

def epsilon_css(j, e):
    '''
    Returns e_css(j), the effective noise strength of a level-j CSS
    gadget.
    '''
    
    e_css = 0
    for r in ('c', 't'):
        for m in (xType, zType):
            # r_m(0,j^~f|j) = 3*e + c3_rm * e_m(0,j|j)
            f = rF = 0
            rFbar = 3*e + c3(m,r) * epsilon(m, 0, j, j, e)
            
            for i in range(j):
                f, rF, rFbar = decode(f, rF, rFbar)
                propagated = c3(m,r) * epsilon(m, i+1, j, j, e)
                rF += propagated
                rFbar += propagated
            
            # r_m(j,j|j) = r_m(j,j^f|j) + r_m(j,j^~f|j)
            e_css += rF + rFbar
        
    return e_css

@memoize
def c1(m, i, j, e):
    '''
    See Eq. 37.
    '''
    
    # For even j, swap the coefficients (see section VII A)
    if not (j % 2):
        m = dualType(m)
    
    if 0 == i:    
        # For level-0, use optimization from Eq. (37)
        c1m = 16 * (epsilon(m, i, i, i, e) + e)
        if zType == m:
            c1m += 1
    else:
        c1m = {xType: 1, zType: 2}[m]
    
    return c1m
    
@memoize
def c2(m, j): 
    c2m = {xType: 4, zType: 3}
    
    if not (j % 2):
        m = dualType(m)
    
    return c2m[m]

@memoize
def c3(m, r):
    c3mr = {xType: {'c': 2, 't': 3}, zType: {'c': 3, 't': 3}}
    return c3mr[m][r]
        
@memoize
def p_accept(j, e):
    '''
    Returns p(j|j-1), the probability of accepting a j-BP given acceptance of all input (j-1)-BPs.
    '''
    
    
    if 1 == j:
        # p(1) >= (1-e)^N
        N = 72
        p = (1-e)**N
    else:
        # p(j|j-1) >= 1 - sum_m (2 * f_m^b(j,j|j-1) + (8 * f_m^s(j-1,j|j-1))
        
        p = 1
        for m in (xType, zType):
            
            fmb, bmF, bmFbar = blockTeleport(m, j-1, j, e)
                
            # f^b_m(j,j|j-1)
            fmb, _, _ = decode(fmb, bmF, bmFbar)
            logger.debug('f^b_{0}({1},{2}|{3})={4}'.format(m,j,j,j-1,fmb))
                
            if 2 == j:
                # No sub-block teleportation at level 2.
                fms = 0
                
            else:
                fms, _, _ = subblockTeleport(m, j-1, j, e)
                
            p -=  (2 * fmb + 8 * fms)
        
    logger.debug('p({0}|{1})={2}'.format(j,j-1,p))
    return p 

@memoize
def subblockTeleport(m, i, j, e):
    '''
    Returns {f^s_m(i,j|j-1), s_m(i,j^f|j-1), s_m(i,j^~f|j-1)}, the flagging
    and block error probabilities for level-j subblock teleportation.
    '''
    
    if 2 == j:
        raise RuntimeError('No subblock teleportation at j=2')
    
    f = sF = 0
    sFbar = 3*e + (1 + c1(m,0,j,e)) * epsilon(m, 0, j-1, j-1, e) 

    for k in range(i):
        f, sF, sFbar = decode(f, sF, sFbar)
        propagated = (1 + c1(m,k+1,j,e)) * epsilon(m, k+1, j-1, j-1, e)
        sF += propagated
        sFbar += propagated
    
    logger.debug('f^s_{0}({1},{2}|{3})={4}'.format(m,i,j,j-1,f))
    logger.debug('s_{0}({1},{2},f=1|{3})={4}'.format(m,i,j,j-1,sF))    
    logger.debug('s_{0}({1},{2},f=0|{3})={4}'.format(m,i,j,j-1,sFbar))
    return f, sF, sFbar

@memoize
def blockTeleport(m, i, j, e):
    '''
    Returns {f^b_m(i,j|j-1), b_m(i,j^f|j-1), b_m(i,j^~f|j-1)}, the flagging
    and block error probabilities for level-j block teleportation.
    '''
    f = bF = 0
    bFbar = 4*e + c2(m,j) * epsilon(m, 0, j-1, j-1, e) 

    for k in range(i):
        f, bF, bFbar = decode(f, bF, bFbar)
        eps = epsilon(m, k+1, j-1, j-1, e)
        bF += c2(m,j) * eps
        bFbar += c2(m,j) * eps
    
    logger.debug('f^b_{0}({1},{2}|{3})={4}'.format(m,i,j,j-1,f))
    logger.debug('b_{0}({1},{2},f=1|{3})={4}'.format(m,i,j,j-1,bF))    
    logger.debug('b_{0}({1},{2},f=0|{3})={4}'.format(m,i,j,j-1,bFbar))
    return f, bF, bFbar

def ConcatenationLevel(e, e_targ):
    '''
    Returns the minimum concatenation level required in order to acheive an effective
    error rate of at most 'e_targ' given a physical error rate of 'e'.
    '''
    
    e_eff = e
    j = 0
    while e_eff > e_targ:
        j += 1
        e_eff = epsilon_css(j, e)
    
    return j, e_eff

def BPMultiplicity(j, e, e_targ):
    if 1 == p_accept(j,e):
        return 1
    
#    return math.ceil(math.log(e_targ / j) / math.log(1 - p_accept(j, e)))
    return math.log(e_targ / j) / math.log(1 - p_accept(j, e))

def BPMultiplicities(j, e, e_targ):
    
    multiplicities = [0] * j
    for i in range(j):
        m = BPMultiplicity(j-i, e, e_targ/(j-i))
        multiplicities[j-i-1] = m
        e_targ = (e_targ - e_targ/(j-i)) / (m + BP_j(j-i-1))
        
    return multiplicities

def BPGateOverhead(j, multiplicities):
    
    if 0 == j:
        return A_BP(j)
    
    # First, compute the number of transversal gates
    # There are 5 transversal CNOTs and 4 transversal measurements
    A_trans = A_BP(j)
    
    A = A_trans + BP_j(j) * BPGateOverhead(j-1, multiplicities)

    return multiplicities[j-1] * A

def BPBPOverhead(j, multiplicities):
    
    if 0 == j:
        return 1
    
    if j <= 2:
        # No sub-block teleportations
        A = 12 * BPBPOverhead(j-1, multiplicities)
    else:
        A = (12+8) * BPBPOverhead(j-1, multiplicities)
        
    return multiplicities[j-1] * A

#def LevelJOverhead(j, J, e, e_targ):
#    '''
#    Returns M_j, the number of j-BPs required.
#    '''
#    
#    # TODO: not sure if this is really the correct overhead
#    # calculation.    
#    # Compute the total number of j-BPs required in a 
#    # J-BP.
#    B = 1
#    for i in range(j,J):
#        B *= BP_j(i)
#        
#    if j < J:
#        M_j_plusone = LevelJOverhead(j+1, J, e, e_targ)
#    else:
#        M_j_plusone = 1
#
#    print 'B({0}, {1})={2}'.format(j, J, B)
#    
#    num = math.log(e_targ / (2 * J * B * M_j_plusone))
#    den = math.log(1 - p_accept(j, e))
#    
#    return math.ceil(num / den)

def BP_j(j):
    '''
    Returns the number of j-BPs required to implement a (j+1)-BP.
    '''
    if j < 2:
        # No sub-block teleportations for j <= 2.
        return 12
    else:
        return 12 + 2*4
    
def A_BP(j):
    '''
    Returns the number of gates required to implement a j-BP, not
    including input (and subblock teleportation) (j-1)-BPs.
    '''
    if 0 == j:
        return 3
    
    trans_count = 9
    if 2 >= j:
        # Levels 1 and 2 have rests instead of subblock teleportations
        trans_count = 9 + 2
    
    return trans_count * 4**j

def ComputeOverhead(e, e_targ):
    J, e_eff = ConcatenationLevel(e, e_targ)
#    Mj_list = [LevelJOverhead(j, J, e, e_targ - e_eff) for j in range(1,J+1)]
#    
#    A_J = sum(A_BP(j) * m for j, m in enumerate(Mj_list))


    multiplicities = BPMultiplicities(J, e, (e_targ - e_eff)/2)
    
    # The J-Rec has two J-BPs, four transversal measurements and three transversal CNOTs
    A_J = 2 * BPGateOverhead(J, multiplicities) + 7 * 4**J
        
    return A_J

def ComputeQubitOverhead(e, e_targ):
    J, e_eff = ConcatenationLevel(e, e_targ)

    multiplicities = BPMultiplicities(J, e, (e_targ - e_eff)/2)
    
    # The J-Rec has two J-BPs plus the data blocks
    A_J = 2*BPBPOverhead(J, multiplicities) + 2 * 4**J 
        
    return A_J


def PlotEcssPaccept(epsilons, j_max):
    paccept = []
    ecss = []
    for j in range(1,j_max+1):
        paccept.append([p_accept(j,e) for e in epsilons])
        ecss.append([epsilon_css(j, e) for e in epsilons])
        
    print epsilons
    print 'p(j|j-1):'
    for j, pj in enumerate(paccept):
        print j+1, pj
        
    print 'e_css(j):'
    for j, ecss_j in enumerate(ecss):
        print j+1, ecss_j

    plotList(epsilons, paccept, labelList=('1','2','3','4','5'), xLabel=r'$\epsilon$', yLabel=r'$p(j|j-1)$', legendLoc='lower left', filename='fibonacci-paccept',
             xscale='log')
    
    plotList(epsilons, ecss, labelList=('1','2','3','4','5'), xLabel=r'$\epsilon$', yLabel=r'$e_{css}(j)$', legendLoc='lower right', filename='fibonacci-e_css',
             xscale='log', yscale='log')
    
def PlotOverhead(epsilons, e_targs):
    
    
    overhead = []
    for e_targ in e_targs:
        overhead.append([ComputeOverhead(8/15. * e, e_targ) for e in epsilons])
        
    # Add text to indicate concatenation level    
    def concatentation_level_text(plt):
        plt.text(7e-6, 1e3, "2")
        plt.text(2.8e-5, 9e4, "3")
        plt.text(1.4e-4, 1e7, "4")
        plt.text(4e-4, 3e9, "5")
        plt.text(6e-4, 5e11, "6")
        
    custom = concatentation_level_text
        
    plotList(epsilons, overhead, filename='fibonacci-gate-overhead', labelList=[str(et) for et in e_targs], xLabel=r'$p$', 
#             yLabel='physical locations per logical location', 
             legendLoc='upper left', xscale='log', yscale='log', xlim=(epsilons[0], epsilons[-1]), ylim=(10e1, 10e22),
             custom_command=custom)

def PlotQubitOverhead(epsilons, e_targs):
    
    
    overhead = []
    for e_targ in e_targs:
        overhead.append([ComputeQubitOverhead(8/15. * e, e_targ) for e in epsilons])
        
    # Add text to indicate concatenation level    
    def concatentation_level_text(plt):
        plt.text(1.2e-5, 5e2, "2")
        plt.text(2.8e-5, 7e3, "3")
        plt.text(1.4e-4, 2e6, "4")
        plt.text(4e-4, 3e8, "5")
        plt.text(6e-4, 5e10, "6")
        
    custom = concatentation_level_text
        
    plotList(epsilons, overhead, filename='fibonacci-qubit-overhead', labelList=[str(et) for et in e_targs], xLabel=r'$p$', 
#             yLabel='physical locations per logical location', 
             legendLoc='upper left', xscale='log', yscale='log', xlim=(epsilons[0], epsilons[-1]), ylim=(10e1, 10e22),
             custom_command=custom)


if __name__ == '__main__':
   
    logging.basicConfig()
    logger.setLevel(logging.INFO)
        
    epsilons = []
    e = 1e-4
    while e < 8e-4:
        epsilons.append(e)
        e += 5e-5
        
#    epsilons = (1e-4,)
    j_max = 5
    
    epsilons = [1e-6, 5e-6, 1e-5, 2e-5, 4e-5, 1e-4, 2e-4, 4e-4, 6e-4, 8e-4, 9e-4, 1e-3, 1.1e-3, 1.2e-3]
    e_targs = [1e-12, 1e-10, 1e-9, 1e-6]
    PlotQubitOverhead(epsilons, e_targs)
    PlotOverhead(epsilons, e_targs)
    
    print [A_BP(j) for j in range(1,5)]
    
    print [BP_j(j) for j in range(1,5)]
    
    print [A_BP(j) + BP_j(j)*A_BP(j-1) for j in range(1,5)]
    print [BPGateOverhead(j, [1]*5) for j in range(1,5)]
    
    print [(e,ConcatenationLevel(8/15. * e, 1e-6)) for e in epsilons]