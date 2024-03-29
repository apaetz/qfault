'''
Created on Mar 3, 2010

@author: adam
'''
from error import Pauli, PauliError
import error as error
from qfault.util import bits, listutils

class Qecc(object):
    '''
    Abstract base class for quantum error correcting codes.
    '''

    def __init__(self, name, n, k, d):
        '''
        Constructor
        '''
        self.name = name
        self.n = n
        self.k = k
        self.d = d
    
    def getCorrection(self, e):
        '''
        :param e: The error, given in descending qubit order. 
        :type e: :class:`qec.error.PauliError`
        '''
        return 0
    
    def decodeError(self, e):
        '''
        Returns the Pauli error resulting from perfectly decoding a block
        with error e.  That is, returns the logical error corresponding
        to the physical error e.
        :param e: The error, given in descending qubit order. 
        :type e: :class:`qec.error.PauliError`
        '''
        return e
    
    def detectError(self, e):
        '''
        Returns True if the Pauli error e is detectable and False otherwise.
        '''
        return False
           
    def blockLength(self):
        return self.n
            
    def __str__(self):
        return self.name
    
class QeccNone(Qecc):
    
    def __init__(self, n):
        super(self.__class__, self).__init__('QeccNone', n, n, 0)
             
    
class StabilizerCode(Qecc):
    
    @staticmethod
    def Syndrome(e, generators):
        '''
        Returns the syndrome, in descending generator order, corresponding to error e.
        :param e: The error, given in descending qubit order. 
        :type e: :class:`qec.error.PauliError`
        :param sequence generators: The list of stabilizer generators.
        :rtype: int
        
        >>> from qec.error import Pauli
        >>> e = Pauli.Y + Pauli.X
        >>> generators = [Pauli.X+Pauli.X, Pauli.Z+Pauli.Z, Pauli.Y+Pauli.Z]
        >>> str('{0:b}'.format(StabilizerCode.Syndrome(e, generators)))
        '101'
        '''
        return bits.listToBits((not e.commutesWith(s)) for s in generators)

    def __init__(self, name, n, k, d):
        super(StabilizerCode, self).__init__(name, n, k, d)

    
    def getSyndrome(self, e):
        return self.Syndrome(e, self.stabilizerGenerators())
    
#    def hashError(self, e):
#        return self.Syndrome(e, self.stabilizerGenerators() + self.normalizerGenerators())
#    
#    def unHashError(self, key):
#        normalizers = self.normalizerGenerators()
#        nNorms = len(normalizers)
#        normalizerChecks = key & ((1 << nNorms) - 1)
#        syndrome = key >> nNorms
#        
#        e = self.syndromeCorrection(syndrome)
#        
#        normalizers = self.normalizerGenerators()
#        for i,norm in enumerate(normalizers):
#            if (normalizerChecks & 1) == e.commutesWith(norm):
#                # The normalizer list is formatted so that X and Z
#                # operators are adjacent.  Obtain the dual operator
#                # by flipping the LSB of the index.
#                dual = normalizers[i ^ 1]
#                e *= dual
#            normalizerChecks >>= 1
#            
#        return e
        
    def syndromeLength(self):
        '''
        Returns the number of bits in the syndrome.
        '''
        return len(self.stabilizerGenerators())
       
    def stabilizerGenerators(self):
        '''
        Returns an iterable containing the stabilizer generators.
        Each generator operator is given in descending qubit order.
        '''
        raise NotImplementedError
    
    def normalizerGenerators(self):
        '''
        A sequence of normalizer generators.
        Default ordering is [X1,Z1,X2,Z2,...,Xk,Zk], where Xi (Zi) is the logical X (Z) 
        operator on logical qubit i.  This ordering, however, is not required or enforced.
        Any ordering is permitted as long as it is consistent for each call to this method.
        Each normalizer operator is given in descending qubit order.
        '''
        
        # Assume that all of the logical operators are also in the normalizer.
        # (This is not true for stabilizer states, and so this method must be overridden
        # in that case.)
        logicals = self.logicalOperators()
        normalizers = []
        for k in range(len(logicals)):
            normalizers += tuple(logicals[k][eType] for eType in sorted(logicals[k].keys()))
            
        return tuple(normalizers)
    
    def logicalOperators(self):
        '''
        A sequence of dictionaries of logical operators indexed by [logical qubit][type],
        where the logical qubit can range from 0-k, and type is error.xType or error.zType.
        Each operator is given in descending (physical) qubit order.
        '''
        raise NotImplementedError
    
    def syndromeCorrection(self, s):
        '''
        Returns the Pauli recovery operation corresponding to syndrome s.
        :param int s: The syndrome. Syndrome bits must be given in descending generator order. 
                      Bit i of s corresponds to self.stabilizerGenerators[self.n - i].
        '''
        raise NotImplementedError
    
    def decodeError(self, e):
        '''
        Returns the Pauli error resulting from perfectly decoding a block
        with error e.  That is, returns the logical error corresponding
        to the physical error e.
        :param e: The error, given in descending qubit order. 
        :type e: :class:`qec.error.PauliError`
        '''
       
        syndrome = self.getSyndrome(e)
        corrected = e * self.syndromeCorrection(syndrome)
        
        decoded = Pauli.I
        for logicalOp in self.logicalOperators():
            logicalErr = Pauli.I
            Xl = logicalOp[error.xType]
            Zl = logicalOp[error.zType]
            
            if not corrected.commutesWith(Xl):
                logicalErr *= Pauli.Z
            
            if not corrected.commutesWith(Zl):
                logicalErr *= Pauli.X
                
            decoded += logicalErr
                
        return decoded
    
    def detectSyndrome(self, s):
        '''
        Returns True if the syndrome s represents a detectable error, and returns False otherwise.
        '''
        return 0 != s
    
    def detectError(self, e):
        return self.detectSyndrome(self.getSyndrome(e))
    
#    def logicalOperator(self, qubit, eType):
#        '''
#        Returns the logical operator corresponding to an error of eType,
#        on the given qubit.
#        '''
#        raise NotImplementedError

    def __eq__(self, other):
        if self.stabilizerGenerators() == other.stabilizerGenerators() and\
           self.normalizerGenerators() == other.normalizerGenerators():
            return True
        
        return False
    
class TrivialStablizerCode(StabilizerCode):
    
    def __init__(self):
        super(TrivialStablizerCode, self).__init__('Trivial', 1, 1, 0)
        
    def stabilizerGenerators(self):
        return tuple([])
    
    def logicalOperators(self):
        return ({error.xType: Pauli.X, error.zType: Pauli.Z},)
    
    def syndromeCorrection(self, s):
        return Pauli.I
    
class Codeword(Qecc):
    '''
    This class is intended as an extra interface that is implemented
    by codeword states (e.g. stabilizer states).
    '''
    def get_code(self):
        raise NotImplementedError
    
class StabilizerState(StabilizerCode, Codeword):
    '''
    Abstract class for stabilizer states.
    
    Initialize by specifying the logical operators that are
    added to the stabilizer generators.
    '''
    
    def __init__(self, code, logicalOpTypes):
        if code.k != len(logicalOpTypes):
            raise Exception('Number of logical operators ({0}) does not match k={1}'.format(len(logicalOpTypes), code.k))
        
        name = str(code) + ''.join(logicalOpTypes)
        super(StabilizerState, self).__init__(name, code.n, code.k, code.d)
        
        self.code = code
        
        logicals = code.logicalOperators()
        self.lStabs = tuple(logicals[k][etype] for k,etype in enumerate(logicalOpTypes))

    def stabilizerGenerators(self):
        return self.code.stabilizerGenerators() + self.lStabs
    
    def normalizerGenerators(self):
        return tuple([])
    
    def logicalOperators(self):
        return self.code.logicalOperators()
        
    def get_code(self):
        return self.code
            
class CssCode(StabilizerCode):
    def __init__(self, name, n, k, d):
        super(CssCode, self).__init__(name, n, k, d)
 
    def stabilizerGenerators(self, types=(error.xType, error.zType)):
        raise NotImplementedError
        
#    def getSyndrome(self, e, types=(error.xType, error.zType)):
#        '''
#        Since X and Z stabilizers are distinct for CSS codes, it
#        is possible to return just the 'X-part' or just the 'Z-part'
#        of the syndrome.  To return just the 'X-part', use paulis=[Pauli.X],
#        and similiarly for just the 'Z-part'.
#        '''
#        dualTypes = tuple(error.dualType(t) for t in types)
#        return StabilizerCode.Syndrome(e, self.stabilizerGenerators(dualTypes))
    
    def syndromeLength(self, types=(error.xType, error.zType)):
        return len(self.stabilizerGenerators(types))
    
#    def decodeError(self, e):
#        bits = self._isLogicalError(e[error.xType], Pauli.X) + 2*self._isLogicalError(e[error.zType], Pauli.Z)
#        return Pauli.FromBits(bits)
    
#    def _isLogicalError(self, e):
#        '''
#        Subclass hook.
#        Returns True if the string e represents a logical error
#        of type 'eType', and False otherwise.
#        '''
#        raise NotImplementedError
    
# TODO: create a dummy CSS code for testing.

#class CssState(CssCode, Codeword):
#    '''
#    Abstract class for stabilizer states.
#    
#    Initialize by specifying the logical operators that are
#    added to the stabilizer generators.
#    '''
#    
#    def __init__(self, code, logicalOpTypes):
#        if code.k != len(logicalOpTypes):
#            raise Exception('Number of logical operators ({0}) does not match k={1}'.format(len(logicalOpTypes), code.k))
#        
#        name = str(code) + ''.join(logicalOpTypes)
#        super(CssState, self).__init__(name, code.n, code.k, code.d)
#        
#        self.lStabs = {error.xType: [], error.zType: []} 
#        for i,ltype in enumerate(logicalOpTypes):
#            self.lStabs[ltype].append(code.logicalOperator(i, ltype)) 
#        self.code = code
#        
#    def stabilizerGenerators(self, types=(error.xType, error.zType)):
#        return self.code.stabilizerGenerators(types) + self.logicalStabilizers(types)
#
#    def logicalStabilizers(self, types=(error.xType, error.zType)):
#        stabs = []
#        for t in types:
#            stabs += self.lStabs[t]
#        return tuple(stabs)
#               
#    def logicalOperator(self, qubit, eType):
#        return self.code.logicalOperator(qubit,eType)
#
#    def get_code(self):
#        return self.code


class ConcatenatedCode(StabilizerCode):
    
    def __init__(self, top, bottom):
        assert(1 == bottom.k)  # Only support one-qubit concatenation for now.
        
        n = top.n * bottom.n
        k = top.k
        d = top.d * bottom.d
        name = str(top) + "^" + str(bottom)
        super(ConcatenatedCode, self).__init__(name, n, k, d)
        
        self._bottom = bottom
        self._top = top
        
    def top(self):
        return self._top
    
    def bottom(self):
        return self._bottom
    
    def stabilizerGeneratorsBottom(self, qubit):
        Id = Pauli.I << self.n-1
        generators = [Id * (gen << self._bottom.n * qubit) for gen in self._bottom.stabilizerGenerators()]
        return tuple(generators)
    
    def stabilizerGeneratorsTop(self):
        return tuple(self._expandTopOperator(gen) for gen in self._top.stabilizerGenerators())

    def stabilizerGenerators(self):
        '''
        Returns the stabilizer generators.  The bottom code generators are given first, ordered by
        ascending qubit number.  The top code stabilizers are given last.
        '''
        generators = ()
        for q in range(self._top.n):
            generators += self.stabilizerGeneratorsBottom(q)
                            
        generators += self.stabilizerGeneratorsTop()
        return generators
    
    def bottomSyndromes(self, s):
        '''
        Returns a tuple of bottom code syndromes in ascending qubit order.
        '''
        # Get rid of the top code syndrome
        s >>= (self._top.n - self._top.k)
        
        l = self._bottom.n - self._bottom.k
        mask = bits.lsbMask(l)
        return tuple((s >> (l*q)) & mask for q in range(self._top.n))
    
    def topSyndrome(self, s):
        '''
        Returns the top code syndrome.
        '''
        return s & bits.lsbMask(self._top.n - self._top.k)

    
    def _expandTopOperator(self, operator):
        nBot = self._bottom.n
        logBottom = self._bottom.logicalOperators()
        pauliList = operator.asList()
        
        op = PauliError(0)
        for pauli in pauliList:
            subOp = PauliError(nBot)
            for etype in pauli.types():
                subOp ^= logBottom[0][etype]
            op += subOp

        return op

    def normalizerGenerators(self):
        return tuple(self._expandTopOperator(norm) for norm in self._top.normalizerGenerators())


    def logicalOperators(self):
        logicals = []
        for qubitOps in self._top.logicalOperators():
            logicals.append({etype: self._expandTopOperator(op) for etype, op in qubitOps.iteritems()})

        return logicals

    def detectSyndrome(self, s):
        bLen = (self._bottom.n - self._bottom.k)
        tLen = (self._top.n - self._top.k)
        l = bLen * self._top.n + tLen
        sbits = bits.bitsToList(s, l)
        
        syndromes = listutils.chop(sbits, ([bLen] * self._top.n) + [tLen])
        
        
    
class BellState(CssCode):
    
    def __init__(self, code, plusState, zeroState):
        super(BellState, self).__init__('Bell-' + str(plusState) + str(zeroState),
                                        plusState.n + zeroState.n,
                                        plusState.k + zeroState.k,
                                        min(plusState.d, zeroState.d))
        self.code = code
        self.lstabs = {
                       error.xType: [(x,x) for x in plusState.logicalStabilizers(error.xType)] + \
                                    [(Pauli.I,x) for x in plusState.logicalStabilizers(error.xType)],
                       error.zType: [(z,Pauli.I) for z in plusState.logicalStabilizers(error.zType)] + \
                                    [(z,z) for z in plusState.logicalStabilizers(error.zType)]
                      }


    def logicalStabilizers(self, types=(error.xType, error.zType)):
        stabs = []
        for t in types:
            stabs += self.lstabs[t]
        return tuple(stabs)
                   
    def logicalOperator(self, qubit, eType):
        return self.code.logicalOperator(qubit % self.code.k,eType)

    def stabilizerGenerators(self, types=(error.xType, error.zType)):
        raise NotImplementedError
        
    def getSyndrome(self, e, types=(error.xType, error.zType)):
        dualTypes = tuple(error.dualType(t) for t in types)
        return tuple(StabilizerCode.Syndrome(e[0], self.stabilizerGenerators(dualTypes)))
    
    def syndromeLength(self, types=(error.xType, error.zType)):
        return len(self.stabilizerGenerators(types))
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()