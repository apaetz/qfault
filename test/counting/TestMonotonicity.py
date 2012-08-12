'''
Created on 2011-04-04

@author: adam
'''
import unittest
import GolayCounting
import component
from counting.count_errors import CountResult
from counting.probability import countResultAsPoly
from util.polynomial import chebyshevT


class TestMonotonicityPrBad(unittest.TestCase):

	def __init__(self, result, kGood, gMin, gMax, testMethod='testPrBadD1'):
		super(TestMonotonicityPrBad, self).__init__(testMethod)
		self.kGood = kGood
		self.gMin = gMin
		self.gMax = gMax
		self.result = result

	def testPrBadD1(self):
		'''
		This test checks the non-negativity condition for the first deriviative of the upper bound for Pr[bad].
		'''
		B = calcB(self.result.locTotals)
		assert (self.kGood+1)/self.gMax >= B


class TestMonotonicityPrAccept(unittest.TestCase):

	def __init__(self, result, kGood, gMin, gMax, testMethod='testPrAcceptD1'):
		super(TestMonotonicityPrAccept, self).__init__(testMethod)
		self.kGood = kGood
		self.gMin = gMin
		self.gMax = gMax
		self.result = result

	def testPrAcceptD1(self):
		gMin = self.gMin
		gMax = self.gMax
		counts = self.result.countsRejected
		B = calcB(self.result.locTotals)
		term = lambda k: k*(calcCPrime(min, gMax, counts,k) + calcCPrime(max, gMin,counts,k)) - \
		                 B * (calcC(max, gMax, counts,k) + calcC(min, gMin, counts, k))
		terms = [term(k) for k in range(1,len(counts))]
		assert sum(terms) >= 0
		
class TestMonotonicityPrMalig(unittest.TestCase):

	def __init__(self, result, noise, kGood, gMin, gMax, testMethod='testPrMaligD4'):
		super(TestMonotonicityPrMalig, self).__init__(testMethod)
		self.kGood = kGood
		self.gMin = gMin
		self.gMax = gMax
		self.result = result
		self.noise = noise

		# Don't care about Pr[bad] for this test.		
		self.result = CountResult(result.counts, result.countsRejected, result.locTotals, result.prAccept, 0)
		
	def testMinOrder(self):
		pass

	def testPrMaligD4(self):
		gMin = self.gMin
		gMax = self.gMax
		result = self.result
		prMalig = countResultAsPoly(self.result, self.noise)
		epsilon = max(gMax/1000, gMin)
		
		pMax = lambda xMin, xMax: TestMonotonicityPrMalig.prMaligUB(result.counts, result.locTotals, xMin, xMax, result.prAccept)
		assert isMonotonic(prMalig, pMax, sum(result.locTotals), epsilon, gMax)
	
		prUBEpsilon = TestMonotonicityPrMalig.prMaligUB(result.counts, result.locTotals, gMin, epsilon, result.prAccept)
		Tn = chebyshevT(sum(result.locTotals))
		Tn5 = Tn.diff(order=5)
		Tn5_1 = Tn5(1)
		d5Max = Tn5_1 * prUBEpsilon
		d4 = prMalig.diff(order=4)
		
		assert d4(0) > abs(d5Max * (epsilon-gMin))
				
	@staticmethod
	def prMaligUB(counts, n, gMin, gMax, prAccept):
		term = lambda k: calcC(max, gMax, counts, k) + calcC(min, gMin, counts, k)
		return calcA(n, gMin) * sum(term(k) for k in range(len(counts))) / prAccept(gMax)
		
#		
#	def testPrBadD2(self):
#		TestMonotonicity.CheckD2(self.kGood, self.gMin, self.gMax, self.result.locTotals)
#
#	def testPrAcceptD1(self):
#		rejectCounts = self.result.countsRejected
#		
#		# Not all components involve postselection.
#		if (None != rejectCounts) and (0 != rejectCounts):
#			print rejectCounts
#			print self.result.prAccept
#			TestMonotonicity.CheckNonNegativity(rejectCounts)
#			TestMonotonicity.CheckD1Tight(self.kGood, self.gMin, self.gMax, self.result.locTotals, rejectCounts)
#
#	def testPrAcceptD2(self):
#		rejectCounts = self.result.countsRejected
#		
#		# Not all components involve postselection.
#		if (None != rejectCounts) and (0 != rejectCounts):
#			#TestMonotonicity.CheckD2Tight(self.kGood, self.gMin, self.gMax, self.result.locTotals, rejectCounts)
#			TestMonotonicity.CheckD2(0, self.gMin, self.gMax, self.result.locTotals)
#
#	@staticmethod
#	def CheckD1(kGood, gMin, gMax, n):
#		B = TestMonotonicity.CalcB(n)
#		delta = TestMonotonicity.CalcDelta(gMin, gMax)
#		assert kGood + 1 >= B * gMax / delta
#		
#	@staticmethod
#	def CheckD1Tight(kGood, gMin, gMax, n, counts):
#		num = sum(counts[k] * (gMin/(1-12*gMin))**k * (k/gMax)/(1-12*gMin) for k in range(len(counts)))
#		den = sum(counts[k] * (gMin/(1-12*gMin))**k for k in range(len(counts)))
#		
#		kPrime = num/den - 1
#		TestMonotonicity.CheckD1(kPrime, gMin, gMax, n)
#		
#	@staticmethod
#	def CheckD2(kGood, gMin, gMax, n):
#		B = TestMonotonicity.CalcB(n)
#		C = TestMonotonicity.CalcC(n)
#		delta = TestMonotonicity.CalcDelta(gMin, gMax)
#		#assert kGood * (kGood + 1) >= (gMax / delta)**2 * ((B*(1-12*gMax)/(1-4*gMin))**2 - C)
#		assert (B*(1-12*gMax))**2 + ((kGood+1)/gMax)**2 >= C + 2*B*(kGood + 1)/gMax
#
#	@staticmethod
#	def CheckD2Tight(kGood, gMin, gMax, n, counts):
#		num = sum(counts[k] * (gMin/(1-12*gMin))**k * (k-1)/(gMax * (1-12*gMin))**2 for k in range(len(counts)))
#		den = sum(counts[k] * (gMin/(1-12*gMin))**k for k in range(len(counts)))
#		
#		kPrime = num/den - 1
#		TestMonotonicity.CheckD2(kPrime, gMin, gMax, n)
#		
#
#		
#
#	
#	@staticmethod
#	def CalcC(n):	
#		nMeas = n.measX + n.measZ
#		nPrep = n.prepX + n.prepZ
#		return 144*n.cnot + 64*n.rest + 16*(nPrep + nMeas)
#	
#	@staticmethod
#	def CalcDelta(gMin, gMax):
#		#return (1-12*gMax)/(1-12*gMin)
#		return (1-12*gMax)
	
#class TestMonotonicityWithEvents(TestMonotonicity):
#	
#	def testPrEventD1(self):
#		TestMonotonicity.CheckD1Tight(self.kGood, self.gMin, self.gMax, self.result.locTotals, self.result.counts)
#		
#	def testPrEventD2(self):
#		print self.result.counts
#		TestMonotonicity.CheckNonNegativity(self.result.counts)
#		#TestMonotonicity.CheckD2Tight(self.kGood, self.gMin, self.gMax, self.result.locTotals, self.counts)
#		TestMonotonicity.CheckD2(0, self.gMin, self.gMax, self.result.locTotals)
		
		
def calcA(n, gamma):
	nMeas = n.measX + n.measZ
	nPrep = n.prepX + n.prepZ
	nPM = nMeas + nPrep
	return (1-12*gamma)**n.cnot * (1-8*gamma)**n.rest * (1-4*gamma)**(nPM)
			
def calcB(n):
	nMeas = n.measX + n.measZ
	nPrep = n.prepX + n.prepZ
	return 12*n.cnot + 8*n.rest + 4*(nPrep + nMeas)
		
def calcC(minOrMax, gamma, counts,k):
	Gamma = gamma/(1-12*gamma)
	return minOrMax(counts[k],0) * (Gamma**k)

def calcCPrime(minOrMax, gamma, counts,k):
	return minOrMax(counts[k],0) * gamma**(k-1) / (1-12*gamma)**k


def isMonotonic(p, pMax, maxDegree, xMin, xMax):
	'''
	Tests whether univariate polynomial p with maximum degree maxDegree 
	is monotone non-decreasing over the interval [xMin, xMax].
	pMax must be a function that upper bounds the absolute value of p over
	the interval and takes two arguments (xMin, xMax).
	'''
	dP = p.diff()
	Tn = chebyshevT(maxDegree)
	d2Tn = Tn.diff(2)
	d2Tn1 = d2Tn(1)
	
	return checkMonotonicity(pMax, dP, xMin, xMax, d2Tn1)

def getMinDegree(p, maxDegree):
	dP = p
	deg = 0
	while (p(0) == 0) and (deg <= maxDegree):
		dP = dP.diff()
		deg += 1
		
	return deg

def checkMonotonicity(pMax, dP, xMin, xMax, d2Tn1):
	'''
	Recursively checks for monotonicity of the polynomial with
	derivative dP over the interval [xMin, xMax] by breaking the
	interval into smaller sub-intervals.
	'''
	#print 'checking [{0},{1}]'.format(xMin, xMax)	
	dP0 = dP(xMin)
	pMaxInterval = pMax(xMin,xMax)
	if abs(d2Tn1 * pMaxInterval) <= dP0:
		return True
	
	# Base case.  Prevent possible infinite recursion.
	# Value of 1e-10 is arbitrary.
	if abs(xMax-xMin) < 1e-10:
		return False
	
	# Monotonicity condition failed over this interval.  Break it into
	# two sub-intervals and check again.
	return checkMonotonicity(pMax, dP, xMin, xMax/2, d2Tn1) and checkMonotonicity(pMax, dP, xMax/2, xMax, d2Tn1)
	

def golayRandomSuite():
	preps = GolayCounting.getSteaneRandomPreps()
	settings = GolayCounting.globalSettings
	return monotonicitySuite(preps, settings)
	
def monotonicitySuite(preps,settings):
	settingsEC = settings.child('ec')
	settingsVZ = settingsEC.child('vz')
	settingsVX = settingsVZ.child('vx')
	settingsPrep = settingsVX.child('zero')
	
	noise = settings['noise']
	gMin = 0
	#gMax = settings['pMax']/15
	gMax = .0018/15
	
	pair1 = (preps[0], preps[1])
	pair2 = (preps[2], preps[3])
	
	verifyX1_X = component.xverify.countXVerifyXOnly(preps[0], preps[1], settingsVX, noise)
	verifyX2_X = component.xverify.countXVerifyXOnly(preps[2], preps[3], settingsVX, noise)
	#verifyX1_Z = component.xverify.countXVerifyZOnly(preps[0], preps[1], settingsVX, noise)
	#verifyX2_Z = component.xverify.countXVerifyZOnly(preps[2], preps[3], settingsVX, noise)
	#verifyZ_X = component.zverify.countZVerifyXOnly(pair1, pair2, settingsVZ, noise)
	verifyZ_Z = component.zverify.countZVerifyZOnly(pair1, pair2, settingsVZ, noise)
	#ec_X = component.ec.countEC_xOnly((pair1,pair2), settingsEC, noise)
	#ec_Z = component.ec.countEC_zOnly((pair1,pair2), settingsEC, noise)
	exRec_X = component.exrec.exRecXOnly((pair1,pair2), settings)
	exRec_Z = component.exrec.exRecZOnly((pair1,pair2), settings)
	
	suite = unittest.TestSuite()
	loader = unittest.TestLoader()
	testNames = loader.getTestCaseNames(TestMonotonicityPrBad)
	for prep in preps:
		#prep = prep.filterAgainst('prepX').filterAgainst('measX')
		#print prep.getTotals()
		result = CountResult(None,None,prep.getTotals(),0,0)
		tests = [TestMonotonicityPrBad(result, settingsPrep['kGood'], gMin, gMax, test) \
				 for test in testNames]
		suite.addTests(tests)
	
	testNames = loader.getTestCaseNames(TestMonotonicityPrAccept)	
	tests = [TestMonotonicityPrAccept(verifyX1_X, settingsVX['kGood'], gMin, gMax, test)
			 for test in testNames]
	suite.addTests(tests)
	
	tests = [TestMonotonicityPrAccept(verifyX2_X, settingsVX['kGood'], gMin, gMax, test)
			 for test in testNames]
	suite.addTests(tests)
	
	tests = [TestMonotonicityPrAccept(verifyZ_Z, settingsVZ['kGood'], gMin, gMax, test)
			 for test in testNames]
	suite.addTests(tests)

#	tests = [TestMonotonicity(ec_X, settingsEC['kGood'], gMin, gMax, test)
#			 for test in testNames]
#	suite.addTests(tests)
#	
#	resultLEC = CountResult(None,None,ec_X.locTotals,0,0)
#	tests = [TestMonotonicity(resultLEC, settings['kGood-LEC-CNOT'][0], gMin, gMax, test)
#			 for test in testNames]
#	suite.addTests(tests)
#	
#	resultLEC = CountResult(None,None,ec_X.locTotals,0,0)
#	tests = [TestMonotonicity(resultLEC, settings['kGood-LEC-CNOT'][1], gMin, gMax, test)
#			 for test in testNames]
#	suite.addTests(tests)

	
#	tests = [TestMonotonicity(verifyX1_Z, settingsVX['kGood'], gMin, gMax, test)
#			 for test in testNames]
#	suite.addTests(tests)
#	
#	tests = [TestMonotonicity(verifyX2_Z, settingsVX['kGood'], gMin, gMax, test)
#			 for test in testNames]
#	suite.addTests(tests)
#	
#	tests = [TestMonotonicity(verifyZ_Z, settingsVZ['kGood'], gMin, gMax, test)
#			 for test in testNames]
#	suite.addTests(tests)

#	tests = [TestMonotonicity(ec_Z, settingsEC['kGood'], gMin, gMax, test)
#			 for test in testNames]
#	suite.addTests(tests)
#
#	resultLEC = CountResult(None,None,ec_Z.locTotals,0,0)
#	tests = [TestMonotonicity(resultLEC, settings['kGood-LEC-CNOT'][0], gMin, gMax, test)
#			 for test in testNames]
#	suite.addTests(tests)
#	
#	resultLEC = CountResult(None,None,ec_Z.locTotals,0,0)
#	tests = [TestMonotonicity(resultLEC, settings['kGood-LEC-CNOT'][1], gMin, gMax, test)
#			 for test in testNames]
#	suite.addTests(tests)
#
#
	testNames = loader.getTestCaseNames(TestMonotonicityPrMalig)
	for E in exRec_X:
		for result in [E[0]]:
			tests = [TestMonotonicityPrMalig(result, noise['X'], settings['kGood'], gMin, gMax, test)
					 for test in testNames]
			suite.addTests(tests)
	for E in exRec_Z:
		for result in [E[0]]:
			tests = [TestMonotonicityPrMalig(result, noise['Z'], settings['kGood'], gMin, gMax, test)
					 for test in testNames]
			suite.addTests(tests)
	
	

	return suite


if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testName']
	
	# Evaluation of Chebyshev polynomials was exceeding default recursion limit.
	#import sys
	#sys.setrecursionlimit(10000)
	
	suite = golayRandomSuite()
	unittest.TextTestRunner(verbosity=2).run(suite)