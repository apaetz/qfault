from util.counterUtils import SubsetIterator
from golayCode import permuteList, permuteListByInverse
#from GolayOverlapSuperXpairs import superXpairs
from counting.countErrors import reduceAllErrorSyndromesZero, countErrors1Block
from settings.noise import NoiseModelXZSympy
from counting.location import Locations
from counting.countParallel import configureMultiProcess

prep = [[[9, 22], [4, 21], [7, 20], [1, 19], [3, 17], [5, 16], [0, 14], [6, 13]], [[4, 22], [8, 21], [2, 20], [6, 19], [9, 16], [1, 15], [7, 14], [10, 13], [3, 12]], [[20, 17], [19, 16], [14, 15], [5, 22], [3, 21], [8, 18], [0, 13], [10, 12], [2, 11]], [[10, 22], [9, 21], [8, 20], [6, 17], [4, 15], [2, 12], [13, 11]], [[22, 20], [21, 19], [17, 18], [16, 14], [15, 12], [5, 11], [8, 13]], [[1, 22], [0, 21], [7, 19], [5, 18], [9, 17], [2, 16], [8, 15], [10, 14], [3, 13], [4, 11]], [[3, 22], [2, 21], [0, 18], [10, 17], [8, 16], [5, 15], [1, 13]]]  # schedule from Appendix A of Adam's project
permA1 = [18, 15, 16, 20, 0, 17, 22, 1, 13, 11, 12, 14, 4, 2, 6, 19, 21, 9, 7, 8, 3, 10, 5]

permsA2A3 = [
([20, 21, 22, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], [2, 8, 22, 0, 7, 6, 18, 4, 12, 9, 11, 16, 19, 1, 17, 13, 3, 15, 21, 14, 10, 5, 20]),
([20, 21, 22, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], [16, 6, 14, 18, 2, 5, 13, 8, 1, 19, 12, 11, 15, 21, 3, 9, 4, 10, 0, 22, 7, 20, 17]),
([20, 21, 22, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], [0, 19, 14, 22, 15, 9, 21, 3, 5, 10, 17, 7, 11, 4, 1, 6, 2, 8, 18, 16, 20, 13, 12]),
([2, 8, 22, 0, 7, 6, 18, 4, 12, 9, 11, 16, 19, 1, 17, 13, 3, 15, 21, 14, 10, 5, 20], [20, 21, 22, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]),
([2, 8, 22, 0, 7, 6, 18, 4, 12, 9, 11, 16, 19, 1, 17, 13, 3, 15, 21, 14, 10, 5, 20], [0, 19, 14, 22, 15, 9, 21, 3, 5, 10, 17, 7, 11, 4, 1, 6, 2, 8, 18, 16, 20, 13, 12]),
([16, 6, 14, 18, 2, 5, 13, 8, 1, 19, 12, 11, 15, 21, 3, 9, 4, 10, 0, 22, 7, 20, 17], [20, 21, 22, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]),
([0, 19, 14, 22, 15, 9, 21, 3, 5, 10, 17, 7, 11, 4, 1, 6, 2, 8, 18, 16, 20, 13, 12], [20, 21, 22, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]),
([0, 19, 14, 22, 15, 9, 21, 3, 5, 10, 17, 7, 11, 4, 1, 6, 2, 8, 18, 16, 20, 13, 12], [2, 8, 22, 0, 7, 6, 18, 4, 12, 9, 11, 16, 19, 1, 17, 13, 3, 15, 21, 14, 10, 5, 20]),
]


bestXZset = (
	([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22], 
	 [20, 0, 11, 19, 8, 4, 15, 12, 22, 18, 16, 2, 14, 7, 1, 10, 9, 3, 21, 5, 13, 17, 6]), 
	([14, 10, 11, 13, 15, 7, 12, 3, 19, 20, 8, 22, 16, 18, 6, 1, 2, 21, 4, 0, 5, 9, 17], 
	 [12, 0, 10, 11, 17, 3, 1, 19, 8, 6, 18, 20, 4, 2, 14, 7, 13, 9, 22, 5, 15, 16, 21])			
)




#def checkSuperX(e1, Xweights):
#	noPermSyndromes = [golayCode.getSyndrome(golayCode.permuteBinaryString(syndrome, range(23)))
#														   for syndrome in range(2**11)]
#	for perm in superXfaultTolerantPerms:
#		permSyndromes = [golayCode.getSyndrome(golayCode.permuteBinaryString(syndrome, perm))
#														   for syndrome in range(2**11)]
#		
#		ancillaVerification0checkXerrorsForGivenErrorOrdersAndSyndromePermutations(e1, noPermSyndromes, permSyndromes, Xweights, True)


mask11 = (2**11) - 1


def getOverlapPrep(perm):
	permPrep = []
	for round in prep:
		permRound = []
		for src, tgt in round:
			permCnot = [perm[src], perm[tgt]]
			permRound.append(permCnot)
		permPrep.append(permRound)

	return permPrep	



def ancillaVerification0checkXerrorsForGivenErrorOrdersAndSyndromePermutations(e1, 
																			   permutedSyndromes1, 
																			   permutedSyndromes2,
																			   Xweights, shortcircuit=10):

									
	# index by [weight][order]									 
	stats = [[0] * 9 for _ in range(9)]
	failures = 0
	foo4 = 0
	for s in range(2**11):
		sp1 = permutedSyndromes1[s]
		sp2 = permutedSyndromes2[s]
		for w in range(2):
			errorOrder = e1[sp1][w] + e1[sp2][w]
			stats[Xweights[s][w]][errorOrder] +=1
			if (s !=0 and w!=0 and errorOrder <= 3 and errorOrder <= Xweights[s][w]):
			#if (errorOrder <= 3 and errorOrder < Xweights[s][w]):
				failures += 1
				if shortcircuit: return ()	  # return infinity as soon as any fault-tolerance violation is found
				if failures >= shortcircuit: return failures

			elif(errorOrder == 4 and errorOrder <= Xweights[s][w]):
				foo4 += 1
	
	#print 'Stats, by error weight:'
	#for i,stat in enumerate(stats):
		#print i, stat
	#print [sum(stats[w][o] for w in range(len(stats))) for o in range(len(stats))]
	#print sum(stats[w][4] for w in range(len(stats)))
	return failures

def checkXFaultTolerance(errorOrders1, errorOrders2, Xweights):

	for s in range(2**11):
		for w in range(2):
			errorOrder = errorOrders1[s][w] + errorOrders2[s][w]
			if (s !=0 and w!=0 and errorOrder <= 3 and errorOrder <= Xweights[s][w]):
			#if (errorOrder <= 3 and errorOrder < Xweights[s][w]):
				return False
	
	return True

def ancillaVerification0checkZerrors(errorOrdersZ1,
									 errorOrdersZ2,
									 Zweights):
	stats = [[0] * 8 for _ in range(4)]
	for s in range(2**11):
		errorOrder = errorOrdersZ1[s] + errorOrdersZ2[s]
		stats[Zweights[s]][errorOrder] +=1
		if errorOrder <= 3 and errorOrder < Zweights[s]:
			return
		elif s != 0 and errorOrder <= 2 and errorOrder == Zweights[s]:
			return
		
	for i, stat in enumerate(stats):
		print i, stat
	return 0

def isPairFaultTolerant(errorOrders1, errorOrders2, weights, maxCorrectableWeight, superOrder=0):
	
	stats = [[0] * 10 for _ in range(max(weights)+1)]
	for s in xrange(len(weights)):
		errorOrder = errorOrders1[s] + errorOrders2[s]
		stats[weights[s]][errorOrder] += 1
		if (s !=0 and errorOrder <= maxCorrectableWeight and errorOrder < weights[s]):
			return False
		elif (s != 0 and errorOrder <= superOrder and errorOrder <= weights[s]):
			return False
		
#		elif s != 0 and errorOrder <= 2 and errorOrder == weights[s]:
#			return
		
	for w, stat in enumerate(stats):
		print w, stat
	return True


def getPermutedSyndromes(perm):
	return [golayCode.getSyndrome(golayCode.permuteBinaryString(syndrome, perm))
													   for syndrome in range(2**11)]

def findXTolerantForPerm(perm, findMax=10, permMax=50000):
	permutedErrorOrders = permuteErrorOrdersX(e1, perm)
	
	found = []
	for permutationIndex in xrange(permMax):
		if permutationIndex%100 == 0: print "permutationIndex:", permutationIndex
		permutation = permuter.pseudorandomBitPermutation()
		if permutationIndex < 23:	   # start by checking each of the 23 cyclic shifts (which are Golay-code symmetries)
				permutation[23-permutationIndex:] = range(permutationIndex)
				permutation[:23-permutationIndex] = range(permutationIndex, 23)
		permutedErrorOrders2 = permuteErrorOrdersX(e1, permutation)
		if isPairFaultTolerant(permutedErrorOrders, permutedErrorOrders2, Xweights, 3, 2):
			found.append(permutation)
			if len(found) == findMax:
				break
	
	return found


def permuteSyndromeIndexedList(list, qubitPerm):
	permutedSyndromes = [golayCode.getSyndrome(golayCode.permuteBinaryString(syndrome, qubitPerm))
													   for syndrome in range(2**11)]
	return [list[permutedSyndromes[s]] for s in range(len(list))]

def combineZ(errorOrdersZ1, errorOrdersZ2):
	eOrdersNew = [8 for _ in range(2**11)]
	for s1, order1 in enumerate(errorOrdersZ1):
		for s2, order2 in enumerate(errorOrdersZ2):
			sNew = s1 ^ s2
			eOrdersNew[sNew] = min(eOrdersNew[sNew], order1 + order2)
			
	return eOrdersNew


def findZTolerantSets(xTolerantPairs, errorOrdersZ):

	combinedErrorOrdersZ = []
	for pair in xTolerantPairs:
		xPermSyndromes = [permuteZSyndromes(perm) for perm in pair]
		permErrorOrders = [permuteListByInverse(errorOrdersZ, syndromes) for syndromes in xPermSyndromes]
		combinedErrorOrdersZ.append(combineZ(permErrorOrders[0], permErrorOrders[1]))
												
	for i1, i2 in SubsetIterator(range(len(combinedErrorOrdersZ)), 2):
		eOrders1 = combinedErrorOrdersZ[i1]
		eOrders2 = combinedErrorOrdersZ[i2]
		if isPairFaultTolerant(eOrders1, eOrders2, Zweights, 3, 2):
			print [xTolerantPairs[i1], xTolerantPairs[i2]]

def sumErrorOrdersByWeight(errorOrders, weights):
	maxOrder = max(errorOrders)
	maxWeight = max(weights)
	weightTable = [0] * (maxWeight+1)
	for w in range(maxWeight+1):
		weightList = [0] * (maxOrder+1)
		for s in [s for s in xrange(len(weights)) if weights[s] == w]:
			weightList[errorOrders[s]] += 1
			
		weightTable[w] = weightList
		
	return weightTable


def permuteXSyndromes(perm):
	parityBit = (1<<11)
	permSyndromes = [golayCode.getSyndrome(golayCode.permuteBinaryString(syndrome, perm))
													   for syndrome in range(2**11)]
	permSyndromes1 = [parityBit + s for s in permSyndromes]
	return permSyndromes + permSyndromes1

def permuteErrorOrdersX(errOrders, perm):
	permSyndromes = permuteXSyndromes(perm)
	permOrders = permuteList(errOrders, permSyndromes)
	return permOrders

def permuteErrorOrdersZ(errOrders, perm):
	permSyndromes = permuteZSyndromes(perm)
	permOrders = permuteList(errOrders, permSyndromes)
	return permOrders



def permuteZSyndromes(perm):
	return [golayCode.getSyndrome(golayCode.permuteBinaryString(syndrome, perm))
													   for syndrome in range(2**11)]

def checkAncillaSet(perms):
	print perms
	permXsyndromes = [permuteXSyndromes(perm) for perm in perms]	
	#permZsyndromes = [permuteZSyndromes(perm) for perm in perms]
	#xErrorOrders = [permuteListByInverse(errorOrdersX, syndromes) for syndromes in permXsyndromes]
	#zErrorOrders = [permuteListByInverse(errorOrdersZ, syndromes) for syndromes in permZsyndromes]
	xErrorOrders = [permuteList(errorOrdersX, syndromes) for syndromes in permXsyndromes]
#	zErrorOrders = [permuteList(errorOrdersZ, syndromes) for syndromes in permZsyndromes]
	
#	zErrorOrders12 = combineZ(zErrorOrders[0], zErrorOrders[1])
#	zErrorOrders34 = combineZ(zErrorOrders[2], zErrorOrders[3])
#	
	#print errorOrdersX
	#print xErrorOrders[0]
	#print xErrorOrders[1]
	#print permXsyndromes[0][498]
	#print permXsyndromes[1][498]
	


	badSyndrome = 3004
	invSyndromes2 = permuteListByInverse(permXsyndromes[2], permXsyndromes[2])
	print permXsyndromes[2]
	for i, s in enumerate(permXsyndromes[2]):
		if s == badSyndrome:
			print '3004 comes from {0}'.format(i)
		if xErrorOrders[2][s] != errorOrdersX[i]:
			raise Exception("bad error order! sPerm={0}, permOrder={1}, s={2}, order={3}".format(s, xErrorOrders[2][s], i, errorOrdersX[i]))
			
	perm2Syndrome = invSyndromes2[badSyndrome]
	print 'perm2syndrome=', perm2Syndrome, 'order=', errorOrdersX[perm2Syndrome]
	print xErrorOrders[2][badSyndrome]
	print xErrorOrders[3][badSyndrome]

	if not isPairFaultTolerant(xErrorOrders[0], xErrorOrders[1], Xweights, 3):
		return False
	if not isPairFaultTolerant(xErrorOrders[2], xErrorOrders[3], Xweights, 3):
		return False
	print 'set is X fault-tolerant.'
#	if not isPairFaultTolerant(zErrorOrders12, zErrorOrders34, Zweights, 3):
#		return False
#	
#	print 'set is Z fault-tolerant'
#	return True

if __name__ == '__main__':
	from util.listutils import nonZeroIndices
	from golay import golayCode
	
	configureMultiProcess(1)
	
	from golay.ancillaPrep import randomAncillaZPrep#, ancilla0checkXerrors
	errorOrdersX = [4 for _ in range(2**12)]
	errorOrdersZ = [4 for _ in range(2**11)]	
	weights = golayCode.errorWeights()
	Xweights = [weights[s & mask11][s>>11] for s in xrange(1<<12)]
	Zweights = [min(w) for w in weights]
	
	cnots = randomAncillaZPrep(prep, range(7))
	reduceAllErrorSyndromesZero(cnots)
	cnots = Locations(cnots, 'GolayPrepOverlap')

#	orderWeightCounts = [[0] * 8 for i in range(8)]
#	for s in range(2**11):
#		for w in range(2):
#			order, weight = errorOrdersX[s][w], Xweights[s][w]
#			orderWeightCounts[order][weight] += 1
#	print orderWeightCounts

#	
	noise = NoiseModelXZSympy()
	for k in range(4):
		countsX = countErrors1Block(k, cnots, '0', 12, 'X', noise)
		for s in nonZeroIndices(countsX):
			errorOrdersX[s] = min(errorOrdersX[s], k)
			
		countsZ = countErrors1Block(k, cnots, '0', 11, 'Z', noise)
		for s in nonZeroIndices(countsZ):
			errorOrdersZ[s] = min(errorOrdersZ[s], k)
		
	weightTable = sumErrorOrdersByWeight(errorOrdersX, Xweights)
	for w, orders in enumerate(weightTable):
		print w, orders

	permuter = golayCode.BitPermutations()
	numPermutationsToTry = 10000
	e1 = errorOrdersX
#	
#	checkAncillaSet([bestXZset[0][0],bestXZset[0][1], bestXZset[1][0], bestXZset[1][1]])
#	raise Exception

#	checkSuperX(e1, Xweights)
	
	#unpermutedSyndromes = range(2**11)
	found = findXTolerantForPerm(range(23), findMax=30)
	print found
	found += [range(23)]
	#found = findXTolerantForPerm(superXfaultTolerantPerms[2])
	
#	findZTolerantSets(superXpairs, errorOrdersZ)
#	raise Exception
#	
#	for _ in xrange(10):
#		perm1 = permuter.pseudorandomBitPermutation()
#		perm2 = findXTolerantForPerm(perm1, 1, 100000)
#		print [perm1, perm2]
#		

#   a1Perm = found.pop()
#	print a1Perm
#	eOrdersZA0VX = combineZ(errorOrdersZ, permuteSyndromeIndexedList(errorOrdersZ, a1Perm))
	
#	for k in range(4):
#		orderK = filter(lambda order: order == k, eOrdersZA0VX)
#		print 'Order', k, ':', len(orderK) 
			

#	found = superXfaultTolerantPerms[1:20]
#	
	foundX = []
	for perm1 in found:
		perm1Orders = permuteErrorOrdersX(errorOrdersX, perm1)
		for perm2 in found:
			perm2Orders = permuteErrorOrdersX(errorOrdersX, perm2)
			if isPairFaultTolerant(perm1Orders, perm2Orders, Xweights, 3, 2):
				foundX.append((perm1, perm2))
	
	print 'Found {0} X-tolerant pairs to try'.format(len(foundX))
	print 'Computing Z error order convolution...'
	eOrdersZForPairs = []
	for pair in foundX:
		eOrdersZForPairs.append(combineZ(permuteErrorOrdersZ(errorOrdersZ, pair[0]),
								permuteErrorOrdersZ(errorOrdersZ, pair[1])))
	pairsAndOrders = zip(foundX, eOrdersZForPairs)
	print 'Finding Z tolerant sets'
	for pair1, eOrdersZ1 in pairsAndOrders:
		for pair2, eOrdersZ2 in pairsAndOrders:
			if isPairFaultTolerant(eOrdersZ1, eOrdersZ2, Zweights, 3, 2):
				print (pair1, pair2)