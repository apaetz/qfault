'''
Created on 2012-04-03

@author: adam
'''
import unittest
import qec.ed422
from qec.qecc import ConcatenatedCode


class TestConcatenatedCode(unittest.TestCase):


    def test412(self):
        code = qec.ed422.ED412Code()
        cat = ConcatenatedCode(code, code)
        
        print cat.stabilizerGenerators()
        print cat.normalizerGenerators()
        print cat.logicalOperators()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()