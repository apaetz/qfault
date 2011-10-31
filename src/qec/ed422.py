'''
Created on 2011-09-09

@author: adam
'''
from counting.location import Locations
from qec.encode.ancilla import ancillaZPrep
from qec.error import Pauli, PauliError
from qec.qecc import CssCode
import qec.error as error
from util import bits
from util.bits import weight
from util.cache import memoize

# Cnots for stabilizer state preparation circuits.
# Cnots are given as a list of rounds.
# Qubits are zero-based.
cnotPreps = {
             (Pauli.Z, Pauli.Z): [[(2,0)], [(0,3), (2,1)]],
             (Pauli.Z, Pauli.X): [[(0,2), (1,3)]],
             (Pauli.X, Pauli.Z): [[(0,1), (2,3)]],
             (Pauli.X, Pauli.X): [[(0,1)], [(3,0), (2,1)]]
            }

def prepare(eigen1, eigen2):
    '''
    Returns the set of locations for preparing an ancilla for which the first logical qubit is the
    +1 eigenstate of eigen1 and the second logical qubit is the +1 eigenstate of eigen2.
    For example, prepare(Pauli.Z, Pauli.X) prepares encoded |0>|+>.
    '''
    return Locations(ancillaZPrep(cnotPreps[(eigen1, eigen2)]), '[[4,2,2]].'+str(eigen1+eigen2))

class ED412Code(CssCode):
    '''
    The [[4,2,2]] code is the smallest quantum error *detecting* code.
    It encodes two qubits, though one of the two is often defined as a "gauge"
    qubit and is ignored.
    There are two stabilizer generators: XXXX, and ZZZZ.
    One set of logical operators is: IIXX, IZIZ (given in descending qubit order).
    The other set of logical operators is: IXIX, IIZZ.
    '''

    _generators = {error.xType: Pauli.X+Pauli.X+Pauli.X+Pauli.X,
                   error.zType: Pauli.Z+Pauli.Z+Pauli.Z+Pauli.Z}
    
    _normalizers = {error.xType: Pauli.I+Pauli.I+Pauli.X+Pauli.X, 
                    error.zType: Pauli.I+Pauli.Z+Pauli.I+Pauli.Z}
    
    _gaugeOperators = {error.xType: Pauli.I+Pauli.X+Pauli.I+Pauli.X, 
                       error.zType: Pauli.I+Pauli.I+Pauli.Z+Pauli.Z}


    
    def __init__(self, gaugeType=error.xType):
        '''
        Constructor
        '''
        super(ED412Code, self).__init__('[[4,2,2]]' + str(gaugeType), 4, 1, 2)
        self._gaugeType = gaugeType
     
    def stabilizerGenerators(self, types=(error.xType, error.zType)):
        return tuple(self._generators[t] for t in types) + (self._gaugeOperators[self._gaugeType],)
    
    def normalizerGenerators(self):
        return tuple(self._normalizers[t] for t in (error.xType, error.zType))
    
    def syndromeCorrection(self, s):
        corr = self._syndromeCorrectionTable(self._gaugeType)[s]
        #print 's=', s, 'corr=', corr
        return corr
    
    @staticmethod
    @memoize
    def _syndromeCorrectionTable(gaugeType):
        # The correction for the non-trivial syndrome is ambiguous, since
        # it can result from any single-qubit error.  Choice of a correction
        # is mostly arbitrary, though some choices are better than others
        # depending on the encoding circuit.  Here we apply a single Pauli
        # to the first qubit.
        # TODO: determine which corrections are optimal.
        zCorr = Pauli.X + (Pauli.I ** 3)
        xCorr = Pauli.Z + (Pauli.I ** 3)
        
        # A non-zero gauge syndrome is corrected by applying the dual
        # gauge operator.
        gaugeCorr = ED412Code._gaugeOperators[error.dualType(gaugeType)]
        
        # Syndrome bits are ordered from MSB to LSB as:
        # XXXX, ZZZZ, gauge
        return [Pauli.I ** 4,      # Commutes with all three stabilizers 
                gaugeCorr,         # anti-commutes with gauge stabilizer
                zCorr,             # anti-commutes with ZZZZ
                gaugeCorr * zCorr, # anti-commutes with {ZZZZ, gaugeStabilizer}
                xCorr,             # anti-commutes with XXXX
                xCorr * gaugeCorr,
                xCorr * zCorr,
                xCorr * zCorr * gaugeCorr
                ]
        
    
#    def reduceError(self, e):
#        r = dict()
#        for pauli in [Pauli.X, Pauli.Z]:
#            if bits.weight(e[pauli], 4) > 2:
#                r[pauli] = e[pauli] ^ 0xF
#
#        return PauliError(r[Pauli.X], r[Pauli.Z])
    
#    def getSyndrome(self, e, types=(error.xType, error.zType)):
#        '''
#        >>> ED412Code().getSyndrome(Pauli.X)
#        2
#        >>> ED412Code().getSyndrome(Pauli.Z)
#        1
#        >>> ED412Code().getSyndrome(Pauli.Y)
#        3
#        >>> ED412Code().getSyndrome(Pauli.Y, error.xType)
#        1
#        >>> ED412Code().getSyndrome(Pauli.Y, error.zType)
#        1
#        >>> ED412Code().getSyndrome(Pauli.X, error.zType)
#        0
#        >>> ED412Code().getSyndrome(Pauli.Z, error.xType)
#        0
#        '''
#        s = 0
#        for etype in types:
#            s <<= 1
#            s += bits.parity(e[etype], 4)
#             
#        return s 
    
#    def syndromeLength(self, types=(error.xType, error.zType)):
#        return 2 * len(types)
    
    def _isLogicalError(self, e, eType):
        # This code detects weight-one errors.  Weight-three errors are
        # equivalent to weight-one errors.  Weight-four errors are 
        # equivalent to weight-zero errors.  Thus the only undetectable
        # logical errors are weight-two.
        return 2 == bits.weight(e)
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()