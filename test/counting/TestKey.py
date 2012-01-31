'''
Created on 2012-01-10

@author: adam
'''
import unittest
from counting.key import convolveKeyCounts

class Test(unittest.TestCase):


    def testConvolveKeyCounts(self):
        counts1 = {(1,0): 2}
        counts2 = {(2,): 3, (4,): 4}
        convolved = convolveKeyCounts(counts1, counts2, [3,1], [3])
        print convolved


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testKey']
    unittest.main()