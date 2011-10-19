'''
Created on Mar 3, 2010

@author: adam
'''
from qec.error import Pauli
import qec.error as error
from util import bits

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
        return 0
    
    def decodeError(self, e):
        return bool(e)
           
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
        return bits.listToBits((not e.commutesWith(s)) for s in generators)

    def __init__(self, name, n, k, d):
        super(StabilizerCode, self).__init__(name, n, k, d)

    
    def getSyndrome(self, e):
        return self.Syndrome(e, self.stabilizerGenerators())
    
    def syndromeCorrection(self, s):
        raise NotImplementedError
    
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
        '''
        raise NotImplementedError
    
    def normalizerGenerators(self):
        '''
        A sequence of normalizer generators ordered by [X1,Z1,X2,Z2,...,Xk,Zk],
        where Xi (Zi) is the logical X (Z) operator on logical qubit i.
        '''
        raise NotImplementedError
    
#    def logicalOperator(self, qubit, eType):
#        '''
#        Returns the logical operator corresponding to an error of eType,
#        on the given qubit.
#        '''
#        raise NotImplementedError
    
class Codeword(Qecc):
    '''
    This class is intended as an extra interface that is implemented
    by codeword states (e.g. stabilizer states).
    '''
    def getCode(self):
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
        
        normalizers = code.normalizerGenerators()
        typeOffset = {error.xType: 0, error.zType: 1}
        self.lStabs = tuple(normalizers[2*i + typeOffset[etype]] for i,etype in enumerate(logicalOpTypes))

    def stabilizerGenerators(self):
        return self.code.stabilizerGenerators() + self.lStabs
    
    def normalizerGenerators(self):
        return []
    
    def logicalStabilizers(self):
        return self.lStabs
        
    def getCode(self):
        return self.code
            
class CssCode(StabilizerCode):
    def __init__(self, name, n, k, d):
        super(CssCode, self).__init__(name, n, k, d)
 
    def stabilizerGenerators(self, types=(error.xType, error.zType)):
        raise NotImplementedError
        
    def getSyndrome(self, e, types=(error.xType, error.zType)):
        '''
        Since X and Z stabilizers are distinct for CSS codes, it
        is possible to return just the 'X-part' or just the 'Z-part'
        of the syndrome.  To return just the 'X-part', use paulis=[Pauli.X],
        and similiarly for just the 'Z-part'.
        '''
        dualTypes = tuple(error.dualType(t) for t in types)
        return StabilizerCode.Syndrome(e, self.stabilizerGenerators(dualTypes))
    
    def syndromeLength(self, types=(error.xType, error.zType)):
        return len(self.stabilizerGenerators(types))
    
    def decodeError(self, e):
        bits = self._isLogicalError(e[Pauli.X], Pauli.X) + 2*self._isLogicalError(e[Pauli.Z], Pauli.Z)
        return Pauli.FromBits(bits)
    
    def _isLogicalError(self, e):
        '''
        Subclass hook.
        Returns True if the string e represents a logical error
        of type 'eType', and False otherwise.
        '''
        raise NotImplementedError
    
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
#    def getCode(self):
#        return self.code

        
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