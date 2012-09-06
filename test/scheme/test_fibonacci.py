'''
Created on 2012-05-01

@author: adam
'''
from counting.key import SyndromeKeyGenerator
from qec import ed422
from qec.error import Pauli, PauliError, xType
from scheme.fibonacci import BP1Max
from settings.noise import CountingNoiseModelXZ
import unittest


class TestFibonacci(unittest.TestCase):
    
    def setUp(self):
        import util.cache
        util.cache.enableFetch(False)


    def testBP1Max(self):
        counting_nm = {Pauli.Y: CountingNoiseModelXZ()}
        bp1 = BP1Max({Pauli.Y: 1})
        result = bp1.count(counting_nm, Pauli.Y)
        print result.counts
        
        code = ed422.ED412Code(gaugeType=None)
        generator = SyndromeKeyGenerator(code, '')
        emap = {}
        for ix in range(5):
            if ix:
                x = 1 << (ix-1)
            else:
                x = 0
            for iz in range(5):
                if iz:
                    z = 1 << (iz-1)
                else:
                    z = 0
                e = PauliError(4, x, z)
                
                key = generator.get_key(e)
                if not emap.has_key(key) or (e.weight() < emap[key].weight()): 
                    emap[key] = e
                        
                
        print emap
        
        for key, count in result.counts[1].iteritems():
            print emap[key[0]], count


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testBP1Max']
    unittest.main()