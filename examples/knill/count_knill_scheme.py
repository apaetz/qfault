'''
Created on 2012-08-11

@author: adam
'''
from qfault.counting import countParallel
from qfault.qec.error import Pauli
from qfault.scheme.knill import KnillScheme
import logging

if __name__ == '__main__':
    countParallel.setPool(countParallel.DummyPool())
#    countParallel.configureMultiProcess(16)
    
    logging.getLogger('counting.threshold').setLevel(logging.DEBUG)
    logging.getLogger('scheme.knill').setLevel(logging.DEBUG)
    
    kPrep = {Pauli.Y: 3}
    kCnot = {Pauli.Y: 3}
    kEC = {Pauli.Y: 4}
    kExRec = {Pauli.Y: 6}
    
    scheme = KnillScheme(kPrep, kCnot, kEC, kExRec)
    scheme.count()