'''
Created on Mar 3, 2010

@author: adam
'''
#from qec.Error import CompoundError

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
        
    def reduceError(self, e, eType):
        return e;
    
    def getCorrection(self, e, eType):
        return 0
    
    def decodeError(self, e, eType):
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
    
class CssCode(StabilizerCode):
    def __init__(self, name, n, k, d):
        super(CssCode, self).__init__(name, n, k, d)
        
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