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

# Cnots for stabilizer state preparation circuits.
# Cnots are given as a list of rounds.
# Qubits are zero-based.
cnotPreps = {
             (Pauli.Z, Pauli.Z): [[(3,0)], [(0,4), (3,2)]],
             (Pauli.Z, Pauli.X): [[(0,2), (1,3)]]  
            }

def prepare(eigen1, eigen2):
    '''
    Returns the set of locations for preparing an ancilla for which the first logical qubit is the
    +1 eigenstate of eigen1 and the second logical qubit is the +1 eigenstate of eigen2.
    For example, prepare(Pauli.Z, Pauli.X) prepares encoded |0>|+>.
    '''
    return Locations(ancillaZPrep(cnotPreps[(eigen1, eigen2)]), '[[4,2,2]].'+str(eigen1+eigen2))

class ED422Code(CssCode):
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
    
    _logicalOps = (
                   # Qubit 1
                   {error.xType: Pauli.I+Pauli.I+Pauli.X+Pauli.X,
                    error.zType: Pauli.I+Pauli.Z+Pauli.I+Pauli.Z},
                   # Qubit 2
                   {error.xType: Pauli.I+Pauli.X+Pauli.I+Pauli.X,
                    error.zType: Pauli.I+Pauli.I+Pauli.Z+Pauli.Z},
                   )

    def __init__(self):
        '''
        Constructor
        '''
        super(ED422Code, self).__init__('[[4,2,2]]', 4, 2, 2)
     
    def stabilizerGenerators(self, types=(error.xType, error.zType)):
        return tuple(self._generators[t] for t in types)
    
    def logicalOperator(self, qubit, eType):
        return self._logicalOps[qubit][eType]
    
#    def reduceError(self, e):
#        r = dict()
#        for pauli in [Pauli.X, Pauli.Z]:
#            if bits.weight(e[pauli], 4) > 2:
#                r[pauli] = e[pauli] ^ 0xF
#
#        return PauliError(r[Pauli.X], r[Pauli.Z])
    
#    def getSyndrome(self, e, types=(error.xType, error.zType)):
#        '''
#        >>> ED422Code().getSyndrome(Pauli.X)
#        2
#        >>> ED422Code().getSyndrome(Pauli.Z)
#        1
#        >>> ED422Code().getSyndrome(Pauli.Y)
#        3
#        >>> ED422Code().getSyndrome(Pauli.Y, error.xType)
#        1
#        >>> ED422Code().getSyndrome(Pauli.Y, error.zType)
#        1
#        >>> ED422Code().getSyndrome(Pauli.X, error.zType)
#        0
#        >>> ED422Code().getSyndrome(Pauli.Z, error.xType)
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