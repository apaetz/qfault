'''
Created on Mar 3, 2010

@author: adam
'''
#from qec.Error import CompoundError
from qec.error import Pauli

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
    def __init__(self, name, n, k, d):
        super(StabilizerCode, self).__init__(name, n, k, d)
    
    def getSyndrome(self, e):
        raise NotImplementedError
    
class CssCode(StabilizerCode):
    def __init__(self, name, n, k, d):
        super(CssCode, self).__init__(name, n, k, d)
        
    def getSyndrome(self, e, paulis=(Pauli.X, Pauli.Z)):
        '''
        Since X and Z stabilizers are distinct for CSS codes, it
        is possible to return just the 'X-part' or just the 'Z-part'
        of the syndrome.  To return just the 'X-part', use paulis=[Pauli.X],
        and similiarly for just the 'Z-part'.
        '''
        raise NotImplementedError
    
    def decodeError(self, e):
        bits = self._isLogicalError(e[Pauli.X], Pauli.X) + 2*self._isLogicalError(e[Pauli.Z], Pauli.Z)
        return Pauli.FromBits(bits)
    
    def _decodeError(self, e, Type):
        '''
        Subclass hook.
        Returns True if the string e represents a logical error
        of type 'eType', and False otherwise.
        '''
        raise NotImplementedError
        
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