'''
Created on Mar 3, 2010

@author: adam
'''
#from qec.Error import CompoundError
from qec.error import Pauli, PauliError
import qec.error as error
import itertools
import operator
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
        
    def reduceError(self, e):
        return e;
    
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
        
    @staticmethod
    def Reduce(e, generators):
        # TODO: more efficient way to compute this?
        w = bits.weight(e)
        genRange = range(len(generators))
        for n in genRange:
            for indices in itertools.combinations(genRange, n):
                e1 = reduce(operator.xor, [generators[i] for i in indices], e)
            w1 = bits.weight(e1)
            if w1 < w:
                e = e1
                w = w1
            
        return e

    def __init__(self, name, n, k, d):
        super(StabilizerCode, self).__init__(name, n, k, d)

    def getSyndrome(self, e):
        return self.Syndrome(e, self.stabilizerGenerators())
    
    def reduceError(self, e):
        return self.Reduce(e, self.stabilizerGenerators())
        
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
    
    def logicalOperator(self, qubit, eType):
        '''
        Returns the logical operator corresponding to an error of eType,
        on the given qubit.
        '''
        raise NotImplementedError
    
class StabilizerState(StabilizerCode):
    '''
    Abstract class for stabilizer states.
    
    Initialize by specifying the logical operators that are
    added to the stabilizer generators.
    '''
    
    def __init__(self, code, logicalOpTypes):
        if code.k != len(logicalOpTypes):
            raise Exception('Number of logical operators ({0}) does not match k={1}'.format(len(logicalOpTypes), code.k))
        
        name = str(code).join(logicalOpTypes)
        super(StabilizerState, self).__init__(name, code.n, code.k, code.d)
        self.logicalOps = [code.logicalOperator(i, ltype) for i,ltype in enumerate(logicalOpTypes)]
        self.code = code
        
    def stabilizerGenerators(self):
        return self.code.stabilizerGenerators() + self.logicalOps
        
    def reduceError(self, e):
        return StabilizerCode.Reduce(self.code.reduceError(e), self.logicalOps)
    
    def getSyndrome(self, e):
        s = self.code.getSyndrome(e)
        return (s << len(self.logicalOps)) + StabilizerCode.Syndrome(e, self.logicalOps)
    
    def logicalOperator(self, qubit, eType):
        return self.code.logicalOperator(qubit,eType)
            
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

class CssState(CssCode):
    '''
    Abstract class for stabilizer states.
    
    Initialize by specifying the logical operators that are
    added to the stabilizer generators.
    '''
    
    def __init__(self, code, logicalOpTypes):
        if code.k != len(logicalOpTypes):
            raise Exception('Number of logical operators ({0}) does not match k={1}'.format(len(logicalOpTypes), code.k))
        
        name = str(code) + ''.join(logicalOpTypes)
        super(CssState, self).__init__(name, code.n, code.k, code.d)
        
        self.logicalOps = {error.xType: [], error.zType: []} 
        for i,ltype in enumerate(logicalOpTypes):
            self.logicalOps[ltype].append(code.logicalOperator(i, ltype)) 
        self.code = code
        
    def stabilizerGenerators(self, types=(error.xType, error.zType)):
        gens = []
        for t in types:
            gens += self.logicalOps[t]
        return self.code.stabilizerGenerators(types) + tuple(gens)
               
    def logicalOperator(self, qubit, eType):
        return self.code.logicalOperator(qubit,eType)


        
#class CompoundCode(Qecc):
#    
#    def __init__(self, code1, code2):
#        super(CompoundCode, self).__init__('{0} : {1}'.format(code1.name, code2.name),
#                                           code1.n + code2.n,
#                                           code1.k + code2.k,
#                                           min(code1.d, code2.d))
#        self.codes = (code1, code2)
#
#    def reduceError(self, e, eType):
#        return CompoundError(self.codes[0].reduceError(e[0], eType), 
#                             self.codes[1].reduceError(e[1], eType))
#    
#    def getCorrection(self, e, eType):
#        return self.codes[0].getCorrection(e[0], eType), self.codes[1].getCorrection(e[1], eType)
#    
#    def decodeError(self, e, eType):
#        return self.codes[0].decodeError(e[0], eType), self.codes[1].decodeError(e[1], eType)    