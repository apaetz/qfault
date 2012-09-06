'''
Created on 2012-08-11

@author: adam
'''
from qfault.qec.error import Pauli
from knill_scheme import KnillScheme
import logging

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('counting.threshold').setLevel(logging.DEBUG)
    logging.getLogger('counting.component').setLevel(logging.DEBUG)
    logging.getLogger('counting.result').setLevel(logging.DEBUG)
#    logging.getLogger('scheme.knill').setLevel(logging.DEBUG)
    
    kPrep = {Pauli.Y: 3}
    kCnot = {Pauli.Y: 3}
    kEC = {Pauli.Y: 4}
    kExRec = {Pauli.Y: 6}

    kPrep = {Pauli.Y: 1}
    kCnot = {Pauli.Y: 1}
    kEC = {Pauli.Y: 1}
    kExRec = {Pauli.Y: 1}
    
    scheme = KnillScheme(kPrep, kCnot, kEC, kExRec)
    scheme.count()