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


    def __init__(self):
        '''
        Constructor
        '''
        super(ED422Code, self).__init__('[[4,2,2]]', 4, 2, 2)
        
    def reduceError(self, e):
        r = dict()
        for pauli in [Pauli.X, Pauli.Z]:
            if bits.weight(e[pauli], 4) > 2:
                r[pauli] = e[pauli] ^ 0xF

        return PauliError(r[Pauli.X], r[Pauli.Z])
    
    def getSyndrome(self, e, types=(error.xType, error.zType)):
        s = 0
        for etype in types:
            s << 1
            s += bits.parity(e[etype], 4)
             
        return s 
    
    def _isLogicalError(self, e, eType):
        # This code detects weight-one errors.  Weight-three errors are
        # equivalent to weight-one errors.  Weight-four errors are 
        # equivalent to weight-zero errors.  Thus the only undetectable
        # logical errors are weight-two.
        return 2 == bits.weight(e, 4)
    
class ED422State(ED422Code):
    '''
    Abstract class for [[4,2,2]] stabilizer states.
    
    Initialize by specifying the logical operators that are
    added to the stabilizer generators.
    '''
    
    def __init__(self, logicalStabX=0, logicalStabZ=0):
        super(ED422State, self).__init__()
        self.logical = dict()
        self.logical[error.xType] = logicalStabX
        self.logical[error.zType] = logicalStabZ
        
    def reduceError(self, e):
        e1 = super(ED422State, self).reduceError(e)
        e2 = super(ED422State, self).reduceError(PauliError(e[Pauli.X] ^ self.logical[Pauli.X],
                                                            e[Pauli.Z] ^ self.logical[Pauli.Z]))
        
        if weight(e2) < weight(e1):
            return e2
        return e1
       
    def getSyndrome(self, e, types=(error.xType, error.zType)):
        s = super(ED422State, self).getSyndrome(e, types)
        for etype in types:
            s << 1
            logical = self.logical[error.dualType(etype)]
            s += bits.parity(e[etype] & logical, 4)
            
        print 'e=', str(e), 's=', s
        return s
    
    def _isLogicalError(self, e, eType):
        # A logical error occurs if the error anticommutes with
        # any of the logical operators (which are *not* stabilizer
        # generators).
        
        # We want the logical operator corresponding to the dual of eType.
        # For example, if eType == Pauli.Z and XXII is the logical-X operator
        # in the stabilizer, then the logical-X operator that we want is
        # XIXI.
        logical = self.logical[error.dualType(eType)] ^ 0b0110
        
        return bits.parity(e & logical, 4)
        
        
    
class ED422ZeroPlus(ED422State):
    '''
    Stabilizer state |0>|+> for the [[4,2,2]] code.
    '''
    
    def __init__(self):
        super(ED422ZeroPlus, self).__init__(logicalStabX=0b0101, logicalStabZ=0b0101)
        self.name += ' |0>|+>'
        
if __name__ == '__main__':
    print prepare(Pauli.Z, Pauli.X)