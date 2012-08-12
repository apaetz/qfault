'''
Created on 2010-04-16

@author: adam
'''
from error import Pauli
from qecc import CssCode
from qfault.util import bits, iteration
from qfault.util.bits import weight
from qfault.util.cache import fetchable
import logging

logger = logging.getLogger('Golay')

class GolayCode(CssCode):
    
    def __init__(self):
        super(GolayCode, self).__init__('Golay', 23, 1, 7)
        
        # TODO: really should make this static        
        self.corrections = GolayCode.generateCorrectionsTable()
    
    def hashError(self, e, eType, logical=False):
    
        """Reduces the input X or Z error modulo the stabilizers and possibly the logical operator.  
        Returns a minimum-Hamming-weight representative.
        """
        syndrome = self.getSyndrome(e)
        # TODO is there an existing package that does bit counting/manipulation?
        parity = bits.parity(e, 23)
        eMin = minErrors[syndrome][parity]
        if not logical: 
            return eMin
        
        # The logical X/Z operator effectively measures the parity of the error.
        # If it is included in the stabilizer, then even/odd parity errors with
        # the same syndrome are equivalent.
        eMinL = minErrors[syndrome][parity^1]
        
        w1 = bits.weight(eMin, 23)
        w2 = bits.weight(eMinL, 23)
        return eMin if w1 < w2 else eMinL
        
    def getCorrection(self, e, eType):
        reduced = self.hashError(e, eType)        
        return self.corrections[self.getSyndrome(reduced)]
    
    def decodeError(self, e, eType):
        # We have a logical error if the result of the correction anti-commutes
        # with the dual logical operator.  For example the error ZI... anti-commutes
        # with the logical X-operator XX...
        # This effectively boils down to calculating the parity of the result of the
        # correction.
        return bool(bits.parity(e ^ self.getCorrection(e, eType), self.n))
    
    @staticmethod
    def getSyndrome(pattern): 
        """Compute the syndrome corresponding to the given pattern.
        
        Returns which parity checks are violated, an 11 bit number from 0 to (1<<10)-1.
        That is, the remainder after dividing the pattern (when considering it as 
        the vector representation of a polynomial) by the generator polynomial, GENPOL.
        In the program this pattern has several meanings: (1) pattern = information
        bits, when constructing the encoding table; (2) pattern = error pattern,
        when constructing the decoding table; and (3) pattern = received vector, to
        obtain its syndrome in decoding.
        """
        
        #TODO Generalize for any stabilizer code with a generator polynomial
        
        X22    = 1<<22                # vector representation of X^22
        X11    = 1<<11                # vector representation of X^11
        MASK12 = (1<<23)-(1<<11)    # auxiliary vector for testing
        GENPOL = 0xc75                # generator polynomial, g(x) = x^11+x^10+x^6+x^5+x^4+x^2+1
        
        aux = X22
        if pattern >= X11:
            while pattern & MASK12: 
                while not (aux & pattern): 
                    aux = aux >> 1
                pattern ^= (aux/X11) * GENPOL
        return pattern
    
    @staticmethod
    def generateCorrectionsTable(): 
        """
        It would take a while to precompute the corrections for all 2^23 different X errors.  
        Instead, we precompute the corrections for all 2^11 syndromes only.  This can be done 
        by simply computing the syndromes for all 0-, 1-, 2- or 3-bit errors, and using that 
        the code is perfect.
        """
        corrections = [-1] * (1<<11)
        
        corrections[0] = 0
        for weight in range(1, 4): 
            for t in iteration.SubsetIterator(range(23), weight): 
                tint = reduce(lambda x, y: x+y, map(lambda x: 1<<x, t))
                corrections[GolayCode.getSyndrome(tint)] = tint
        return corrections
    


        
class GolayZero(GolayCode):
    
    def __init__(self):
        super(self.__class__, self).__init__()
        self.name = 'GolayZero'
        
    def hashError(self, e, eType):
        logical = False
        
        if Pauli.Z == eType:
            logical = True
            
        return super(GolayZero, self).hashError(e, eType, logical)
    
    def getCorrection(self, e, eType):
        # TODO
        if Pauli.Z == eType:
            return super(self.__class__, self).getCorrection(e, eType)
        
        raise Exception
        
class GolayPlus(GolayCode):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.name = 'GolayPlus'
        
    def hashError(self, e, eType):
        logical = False
        
        if Pauli.X == eType:
            logical = True
            
        return super(GolayPlus, self).hashError(e, eType, logical)
    
    def getCorrection(self, e, eType):
        # TODO
        if Pauli.X == eType:
            return super(self.__class__, self,).getCorrection(e, eType)
        
        raise Exception
    
    
    

@fetchable
def calcMinErrors():
    '''
    Returns a list of minimum weight errors indexed by syndrome, and then by parity.
    
    >>> mins = calcMinErrors()
    >>> syndrome = 0
    >>> parity = 1
    >>> mins[syndrome][parity]
    3189
    '''
    # A list of error pairs (even-parity error, odd-parity error), indexed by syndrome.
    # Each pair contains a tuple that records the pattern and weight of the 
    # minimum weight error for that syndrome.
    #minWeightForSyndrome = [[(None,24),(None,24)]] * (1<<11)
    minWeightForSyndrome = [[24, 24] for _ in range (1<<11)]
    minErrorForSyndrome = [[None, None] for _ in range (1<<11)] 
    
    logger.info('Calculating minimum weight error for each syndrome')
    onePercent = (1<<23)/100
    percentComplete = 0
    
    for error in range(1<<23):
        
        if not(error % onePercent):            
            logger.info(str(percentComplete) + ' percent done.')
            percentComplete += 1
            
        rawWeight = weight(error, 23)
        parity = rawWeight % 2
        syndrome = GolayCode.getSyndrome(error)
        
        #print 'error={0:b} weight={1} parity={2} syndrome={3:b}\n'.format(error, rawWeight, parity, syndrome)
        
        currentMin = minWeightForSyndrome[syndrome][parity] 
        
        if(currentMin > rawWeight):
            minWeightForSyndrome[syndrome][parity] = rawWeight
            minErrorForSyndrome[syndrome][parity] = error                
    
    return minErrorForSyndrome

#
# Module Initialization
#
# Load/Pre-Compute the minimum weight error for each syndrome.
#minErrors = calcMinErrors()