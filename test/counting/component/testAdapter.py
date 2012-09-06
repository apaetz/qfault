'''
Created on 2012-01-30

@author: adam
'''
import unittest
import testComponent
from counting.component.adapter import IdealDecoder
from qec import ed422
from qec.error import PauliError
from counting.key import SyndromeKeyGenerator, SyndromeKeyDecoder


class TestAdapter(testComponent.ComponentTestCase):

    def testKeyPropagator(self):
        code = ed422.ED412Code()
        generator = SyndromeKeyGenerator(code, None)
        propagator = self._getComponent(None, code).keyPropagator()
        keyDecoder = SyndromeKeyDecoder(self.trivialCode)
        
        for eX in xrange(1 << code.n):
            for eZ in xrange(1 << code.n):
                pauli = PauliError(code.n, eX, eZ)
                key = (generator.get_key(pauli),)
                decoded = keyDecoder.asPauli(propagator(key)[0])
                if (decoded != code.decodeError(pauli)):
                    raise Exception('pauli={0}, key={3}, decoded={1}, code.decode={2}'.format(pauli, decoded, code.decodeError(pauli), key))
    
    def _getComponent(self, kGood, code):
        return IdealDecoder(code)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()