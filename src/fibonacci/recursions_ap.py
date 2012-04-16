'''
Created on 2012-04-16

@author: adam
'''
from util.cache import memoize
import logging
from util.plotting import plotList

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
def recursionAP_jj(m, j, e0):
    '''
    Returns e_m(j,j|j-1), the probability of a level-j error on the output 
    of a j-BP conditioned on acceptance of the input (j-1)-BPs.
    See AP09 Table II.
    '''
    
    f, bF, bFbar = blockTeleport(m, j-1, j, e0)
        
    # b^~_m(j,j^~f|j-1)
    _, _, bFbar = decode(f, bF, bFbar)
    
    logger.debug('b_{0}({1},{2},f=0|{3})={4}'.format(m,j,j,j-1,bFbar))
    return bFbar
    
@memoize
def recursionAP_j1j(m, j, e0):
    '''
    Returns e_m(j-1,j|j-1), the probability of a level-(j-1) error on the output 
    of a j-BP conditioned on acceptance of the input (j-1)-BPs.
    See AP09 Table III.
    '''
    
    if 1 == j:
        # Special case (Eq. 15): e_m(0,1) <= e0 + c1_m e_m(0,0)
        return e0 + c1(m,0,j,e0) * epsilon(m, 0, 0, 0, e0)
    elif 2 == j:
        # Special case optimization (see Sec. VII A)
        # Here level-1 errors in the 2-BP prep propagate just like level-0 errors
        # in the 1-BP prep.  But we can ignore the physical gate errors.
        return c1(m,1,j,e0) * epsilon(m,1,1,1,e0)
    
    _, _, sFbar = subblockTeleport(m, j-1, j, e0)
    
    return sFbar + epsilon(m, j-1, j-1, j-1, e0)

@memoize
def epsilon(m, i, jOut, jIn, e0):
    '''
    Returns e_m(i,jOut|jIn), the probability of level-i error on the output of a jOut-BP 
    conditioned on acceptance of either:
    1. The input (j-1)-BPs (jIn = jOut-1)
    2. The jOut-BP (jIn = jOut)
    '''
    
    # Sanity check
    if i > jOut or jIn > jOut:
        raise RuntimeError('Invalid parameters: m={0} i={1} jOut={2} jIn={3} e0={4}'.format(m,i,jOut,jIn,e0))
    
    if jIn == jOut:
        if i == jOut and 0 == i:
            #e_m(0,0) = e_0
            epsOut = e0

        elif 2 == jOut and 0 == i:
            # Special case.  There are no sub-block teleportations
            # at level j=2.  Instead, physical qubits undergo 2 rounds of rest noise
            # (and noise from the block teleportation CNOT).
            epsOut = (3*e0 + c1(m,i,jOut,e0) * epsilon(m,i,jOut,jIn-1,e0)) / pj_j1(jOut,e0)
        else:    
            # e_m(i,j|j) = e_m(i,j|j-1) / p(j|j-1)
            epsOut = epsilon(m, i, jOut, jOut-1, e0) / pj_j1(jOut, e0)
            
    elif jIn == jOut - 1:
        if i == jOut:
            # Compute e_m(j,j|j-1) by recursion
            epsOut = recursionAP_jj(m, jOut, e0)
        elif i == jIn:
            # Compute e_m(j-1,j|j-1) by recursion
            epsOut = recursionAP_j1j(m, jOut, e0)
        else:
            # Eq. 23 e_m(i,j|j-1) = e_m(i,j-1|j-1)
            epsOut = epsilon(m, i, jIn, jIn, e0)
                        
    logger.debug('e_{0}({1},{2}|{3})={4}'.format(m, i, jOut, jIn, epsOut))
    
    return epsOut

@memoize
def c1(m, i, j, e0):
    '''
    See Eq. 37.
    '''
    
    # For even j, swap the coefficients (see section VII A)
    if not (j % 2):
        m = dualType(m)
    
    if 0 == i:    
        # For level-0, use optimization from Eq. (37)
        c1m = 16 * (epsilon(m, i, i, i, e0) + e0)
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
def pj_j1(j, e0):
    '''
    Returns p(j|j-1), the probability of accepting a j-BP given acceptance of all input (j-1)-BPs.
    '''
    
    
    if 1 == j:
        # p(1) >= (1-e0)^N
        N = 72
        p = (1-e0)**N
    else:
        # p(j|j-1) >= 1 - sum_m (2 * f_m^b(j,j|j-1) + (8 * f_m^s(j-1,j|j-1))
        
        p = 1
        for m in (xType, zType):
            
            fmb, bmF, bmFbar = blockTeleport(m, j-1, j, e0)
                
            # f^b_m(j,j|j-1)
            fmb, _, _ = decode(fmb, bmF, bmFbar)
            logger.debug('f^b_{0}({1},{2}|{3})={4}'.format(m,j,j,j-1,fmb))
                
            if 2 == j:
                # No sub-block teleportation at level 2.
                fms = 0
                
            else:
                fms, _, _ = subblockTeleport(m, j-1, j, e0)
                
            p -=  (2 * fmb + 8 * fms)
        
    logger.debug('p({0}|{1})={2}'.format(j,j-1,p))
    return p 

@memoize
def subblockTeleport(m, i, j, e0):
    '''
    Returns {f^s_m(i,j|j-1), s_m(i,j^f|j-1), s_m(i,j^~f|j-1)}, the flagging
    and block error probabilities for level-j subblock teleportation.
    '''
    
    if 2 == j:
        raise RuntimeError('No subblock teleportation at j=2')
    
    f = sF = 0
    sFbar = 3*e0 + (1 + c1(m,0,j,e0)) * epsilon(m, 0, j-1, j-1, e0) 

    for k in range(i):
        f, sF, sFbar = decode(f, sF, sFbar)
        propagated = (1 + c1(m,k+1,j,e0)) * epsilon(m, k+1, j-1, j-1, e0)
        sF += propagated
        sFbar += propagated
    
    logger.debug('f^s_{0}({1},{2}|{3})={4}'.format(m,i,j,j-1,f))
    logger.debug('s_{0}({1},{2},f=1|{3})={4}'.format(m,i,j,j-1,sF))    
    logger.debug('s_{0}({1},{2},f=0|{3})={4}'.format(m,i,j,j-1,sFbar))
    return f, sF, sFbar

@memoize
def blockTeleport(m, i, j, e0):
    '''
    Returns {f^b_m(i,j|j-1), b_m(i,j^f|j-1), b_m(i,j^~f|j-1)}, the flagging
    and block error probabilities for level-j block teleportation.
    '''
    f = bF = 0
    bFbar = 4*e0 + c2(m,j) * epsilon(m, 0, j-1, j-1, e0) 
#    logger.debug('b_{0}({1},{2},f=0|{3})={4}'.format(m, 0, j, j-1, bFbar))

    for k in range(i):
        f, bF, bFbar = decode(f, bF, bFbar)
        eps = epsilon(m, k+1, j-1, j-1, e0)
        bF += c2(m,j) * eps
        bFbar += c2(m,j) * eps
    
    logger.debug('f^b_{0}({1},{2}|{3})={4}'.format(m,i,j,j-1,f))
    logger.debug('b_{0}({1},{2},f=1|{3})={4}'.format(m,i,j,j-1,bF))    
    logger.debug('b_{0}({1},{2},f=0|{3})={4}'.format(m,i,j,j-1,bFbar))
    return f, bF, bFbar
        
if __name__ == '__main__':
   
    logging.basicConfig()
    logger.setLevel(logging.DEBUG)
    
    epsilons = (1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4)
#    epsilons = (1e-3,)
    paccept = []
    for j in range(1,6):
        paccept.append([pj_j1(j,e0) for e0 in epsilons])
        
#    print epsilons
#    for pj in paccept:
#        print pj

    plotList(epsilons, paccept, labelList=('1','2','3','4','5'), xLabel=r'$\epsilon$', yLabel=r'$p(j|j-1)$', legendLoc='lower left', filename='fibonacci-paccept')

#    e0 = 0.001
#    print epsilon(xType, 0, 0, 0, e0)
#    print epsilon(xType, 0, 1, 0, e0)
#    print epsilon(xType, 1, 1, 0, e0)
#    print epsilon(xType, 1, 1, 1, e0)

#    print f_s(xType, 2, 2, e0)