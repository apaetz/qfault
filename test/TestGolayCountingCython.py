'''
Created on 2010-11-04

@author: adam
'''
import unittest

from component import cython


class TestCython(unittest.TestCase):


	def followTEC(self):
		'''
		Very simple input data for which a logical error should occur on a block iff the TEC for that
		block was present.
		'''
		# Index as [k][i]=syndrome
		sLECLists = [[0], [0]]
		sCLists = [[0], [0]]
		
		# Index as [k][syndrome]=count
		countsLEC = [[1], [0,1]]
		countsC = [[1], [0,1]]
		
		lCountsTEC = [
						[ [1, 1], [0, 1] ], # k=0
						[ [0, 1], [0, 1] ], # k=1
					 ]
		
		configs = [
					[0, 0, 0, 0, 0]
				  ]
		
		errorStr = ['IX', 'XI', 'XX']
		ecStr = ['--', '-B', 'A-', 'AB']
		
		counts, _ = cython.exrec.countExRecForConfigs_c(sLECLists, sCLists, countsLEC, countsC, lCountsTEC, configs)
		print 'error ec count'
		for error in range(3):
			eA = (error+1) >> 1
			eB = (error+1) & 1
			for ec in range(4):
				ecA = ec >> 1
				ecB = ec & 1
				count = counts[ec][error]
				print errorStr[error], ecStr[ec], count
				assert len(count) == 1
				
				# A logical error on a block should only occur if the TEC on that block exists.
				assert count[0] == ((ecA and eA) and not eB*(eB ^ ecB)) or ((ecB and eB) and not eA*(eA ^ ecA))
				


if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testName']
	unittest.main()