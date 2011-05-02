#!/usr/bin/python
# by Ben Reichardt, created 12/28/2009
# 
# This code simulates a CNOT extended rectangle (exRec).  Given error rates for the different 
# locations, it estimates the probability that the decoders can't be pulled across the rectangle.  
#
# Currently, the code is mostly set up for using with the Golay code, for which counting is computationally 
# expensive.  Because the ancillas used for error correction need to pass verification tests before they can 
# be used, it really does the full simulation, without any precomputation.  
#
# Extensions: 
#	* Some code needs to be refactored, for example moving stuff into Utils.py.
#	* It would be useful to keep track of resources used. 
#	* Look into the quality of the random number generator.
from counting.countParallel import configureMultiProcess, DummyPool

print "Location counter running"

import util.counterUtils
import sim.simulateUtils
from sim.simulateUtils import cnot, measX, measZ, prepX, prepZ, rest


#from golay import golayCode
from golaySim import Corrector, prepareAncilla, prepareAncillaNew, prepareAncillaOverlap,\
	prepareAncillaSteaneRandom
n = 23

prepareAncilla = prepareAncillaOverlap

#import BaconShor9
#from BaconShor9 import Corrector, prepareAncilla
#n = 9

#import BaconShor25
#from BaconShor25 import Corrector, prepareAncilla
#n = 25


import random

def checkMalignantUsingDecoder(errors, corrector):
	"""Checks whether the decoders can be pulled past, as in the AGP paper.
	
	It applies the decoders to the first and last checkpoints, and sees if they agree. 
	"""
	errorsX, errorsZ = decode1RecErrors(errors, corrector)
	isXMalignant = bool(errorsX)
	isZMalignant = bool(errorsZ)
	return isXMalignant, isZMalignant

def decode1RecErrors(errors, corrector):
	"""Determines and decodes logical errors introduced by the 1-Rec.
	It applies the decoders to the first and last checkpoints, and sees if they agree. 
	"""
	decode = lambda block: [corrector.decodeError(errors[etype][block], etype) for etype in ['X', 'Z']]
	Da1, Db1, Da3, Db3 = map(decode, ['Da1', 'Db1', 'Da3', 'Db3'])
	# having decoded the blocks at each checkpoint, we see whether the decoders commute
	# for example, if Da1 has a logical X error, and Db1 doesn't, then Da3 and Db3 both should
	errorX = (Da3[0] != Da1[0]) + ((Db3[0] != (Da1[0] ^ Db1[0])) << 1)
	errorZ = (Da3[1] != (Da1[1] ^ Db1[1])) + ((Db3[1] != Db1[1]) << 1)
		# it is also interesting to see what location sets are malignant just for X errors, or 
		# are simultaneously malignant for X and Z errors; see the following commented-out lines
	#isMalignant = bool(errorX)
	#isMalignant = errorX and errorZ
	#if isMalignant: 
	#	print Da1, Db1, Da3, Db3
	#	print util.counterUtils.errorsToStr(combinedErrors, 9)
	return errorX, errorZ


def sampleExRec(errorRates):
	'''
	Samples from the exRec error distribution by simulating one execution of the exRec.
	Returns the X and Z errors on the two blocks.
	'''
	
	# should move this code into golayCode.py
	# it would also be useful to keep track of expended resources
	ecXa1 = prepareAncilla(errorRates, 'ecXa1', 'X')
	ecZa1 = prepareAncilla(errorRates, 'ecZa1', 'Z')
	ecXb1 = prepareAncilla(errorRates, 'ecXb1', 'X')
	ecZb1 = prepareAncilla(errorRates, 'ecZb1', 'Z')
	ecXa2 = prepareAncilla(errorRates, 'ecXa2', 'X')
	ecZa2 = prepareAncilla(errorRates, 'ecZa2', 'Z')
	ecXb2 = prepareAncilla(errorRates, 'ecXb2', 'X')
	ecZb2 = prepareAncilla(errorRates, 'ecZb2', 'Z')
		# join together all the error blocks
	errors = reduce(util.counterUtils.errorsMerge, [ecXa1, ecZa1, ecXb1, ecZb1, ecXa2, ecZa2, ecXb2, ecZb2])
	
	errors['X']['Da1'] = errors['Z']['Da1'] = 0
	errors['X']['Db1'] = errors['Z']['Db1'] = 0
	
	# now apply the various CNOTs, apply and propagate the error corrections, 
	for i in range(n):
		cnot(errorRates, errors, 'Da1', i, 'ecXa1', i)
		measZ(errorRates, errors, 'ecXa1', i)
		cnot(errorRates, errors, 'ecZa1', i, 'Da1', i)
		measX(errorRates, errors, 'ecZa1', i)
		
		cnot(errorRates, errors, 'Db1', i, 'ecXb1', i)
		measZ(errorRates, errors, 'ecXb1', i)
		cnot(errorRates, errors, 'ecZb1', i, 'Db1', i)
		measX(errorRates, errors, 'ecZb1', i)
	
	#print util.counterUtils.errorsToStr(errors)
	cXa1 = corrector.correctXError(errors['X']['ecXa1'])
	cZa1 = corrector.correctZError(errors['Z']['ecZa1'])
	cXb1 = corrector.correctXError(errors['X']['ecXb1'])
	cZb1 = corrector.correctZError(errors['Z']['ecZb1'])
	errors['X']['Da1'] ^= cXa1
	errors['Z']['Da1'] ^= cZa1
	errors['X']['Db1'] ^= cXb1
	errors['Z']['Db1'] ^= cZb1
	
	# keep a snapshot of the blocks after the first ECs
	errors['X']['Da3'] = errors['X']['Da1']
	errors['Z']['Da3'] = errors['Z']['Da1']
	errors['X']['Db3'] = errors['X']['Db1']
	errors['Z']['Db3'] = errors['Z']['Db1']
	
	for i in range(n):
		cnot(errorRates, errors, 'Da3', i, 'Db3', i)
	
	for i in range(n):
		cnot(errorRates, errors, 'Da3', i, 'ecXa2', i)
		measZ(errorRates, errors, 'ecXa2', i)
		cnot(errorRates, errors, 'ecZa2', i, 'Da3', i)
		measX(errorRates, errors, 'ecZa2', i)
		
		cnot(errorRates, errors, 'Db3', i, 'ecXb2', i)
		measZ(errorRates, errors, 'ecXb2', i)
		cnot(errorRates, errors, 'ecZb2', i, 'Db3', i)
		measX(errorRates, errors, 'ecZb2', i)

	cXa2 = corrector.correctXError(errors['X']['ecXa2'])
	cZa2 = corrector.correctZError(errors['Z']['ecZa2'])
	cXb2 = corrector.correctXError(errors['X']['ecXb2'])
	cZb2 = corrector.correctZError(errors['Z']['ecZb2'])
	errors['X']['Da3'] ^= cXa2
	errors['Z']['Da3'] ^= cZa2
	errors['X']['Db3'] ^= cXb2
	errors['Z']['Db3'] ^= cZb2

	
	return errors

def sample1RecErrors(errorRates):
	return decode1RecErrors(sampleExRec(errorRates), corrector)

def simExRecErrors(sampleSize, errorRates, pool=DummyPool()):
	xErrors = [0] * 4
	zErrors = [0] * 4
	numFailures = 0
	
	results = [pool.apply_async(sample1RecErrors, [errorRates]) for _ in xrange(sampleSize)]
	
	for result in results: 
		xError, zError = result.get()
		xErrors[xError] += 1
		zErrors[zError] += 1
		numFailures += (xError or zError)
		
	return xErrors, zErrors, numFailures


# Initialize the corrections that are applied for each possible error.  
# Also initialize the ideal decoder.  
corrector = Corrector()

errorRates = sim.simulateUtils.ErrorRates()
errorRates.cnot = 3e-3
errorRates.prepX = errorRates.prepZ = errorRates.measX = errorRates.measZ = 4/15. * errorRates.cnot
errorRates.rest = 12/15. * errorRates.cnot
#errorRates.rest = 0


random.seed()
#random.setstate(	# use this to restore a previous session
#(2, (-2147483648, -2048355095, -867227090, 534793798, -1758285229, -826820770, 38506121, 1114308135, -1135163709, -4639976, 1365929529, 1510617060, 2140362810, -750698076, 764378756, 35777602, 1044397518, 225373878, -1255788122, -1329210713, 1684674155, 1395517120, -640969024, 492246621, -1531045018, 1337165509, -367682072, 794818743, 1542060637, -404286835, -782360204, 1487404984, 386729505, -1936284797, -569072332, -2062465719, 2063509930, 231321487, 1965312421, -1030345738, -1186513723, 1811150750, 952141736, -1263218452, -94779527, 1413302211, -220822439, -1954839982, 229945563, -2044693658, 1584471840, 1304268517, -205157626, 2014607069, 751889655, -121922700, 809416665, -255220373, 660762657, 703568864, 1717572488, 640353435, -1206231461, -1277471129, -2007262243, -1383979378, 561877721, -177062326, 204150767, -386654639, -2004069462, -657641553, -1531370317, 615656992, -1909422697, 47539038, 429949898, 602657737, -891911509, 1639400427, 702248563, 1031931460, 1408926269, -1457296267, 264691184, -1713237605, 1273784928, -102442559, -1813447287, 1902061802, 1241317271, 2003232190, -682875032, -863570962, 486087581, -451386957, 21508435, -949250790, -1558829750, 982935171, 952736395, 880631756, 1087152853, -347372375, -551735000, -977230977, 544386597, 1113078946, 1751058162, -1267875743, -1596276710, -1701111232, -1684414751, -646510360, 1847221775, 715531049, -347354698, -1755683349, -1266527347, 533818084, 995441798, -1832853504, -1897448694, 1138335975, 1115622880, -1716657729, 1422497302, -1537245462, -665826465, 871177855, -1899996366, -77421505, 307469838, 1898538100, 1196404346, 753803651, 426410740, 1966136209, -1695073541, 182805065, -2116804760, -1042704670, -928030599, 1827181483, -2142903973, 934908921, -1178201733, -910899416, -1355666298, 1784494142, 1135805516, -788918681, -1232482673, 1587406881, -1365449815, -544826609, 1420178510, -1571203109, -1606315709, -1953200945, 1679822721, -1478242974, 195515121, -136856713, 1881203536, 2095726709, 1313772323, -837951772, -1200214520, -556002046, -1174998930, 1998003214, -1128679144, -230027981, -853296331, 476616207, 1524402770, -1196052144, -1496677533, -1919882779, 1804811759, 105025611, -1934727099, 197150635, -55465518, 687291450, 1810068472, -1417214414, -677184647, 1285096329, 356065487, 1298956984, -1472449028, -1135641826, -42866992, 1392933328, -204175220, -483086844, -1260974485, 514243374, -632519612, -2005923285, 1543983192, -358305859, -1028950752, 965012846, -2052359718, 2072049190, 1756755681, 910385311, 2020600751, -1637503149, 1656954405, 943076836, 920941282, 52060806, 497498833, -2133751175, -85405342, -1094858280, 1384875855, 163157318, 141138689, -771143777, -1040079328, 915779, -315375214, -682424210, -1043482215, -1125364350, 2113287504, 1837488683, 896867912, -1978125602, 199560499, 318801395, 1943932456, 1655073376, -2023719922, 1301396100, -779382479, 2138503732, 285536018, 1044683380, 2060213012, -208187709, 285947827, 1863634479, -1192115825, 231039599, -1756054390, 1745587654, 1150127316, 761652307, -1442408308, -161087406, -779719946, 1102212901, 1320745906, 833329051, 1264563339, 194422441, -230913517, 643747227, -175753376, 133443159, -763618397, -1372021793, 270250490, -1827092917, 2082901201, 389216482, -405217683, -1403163930, -650450756, 12121034, -1644447506, -201566826, 1348837162, 1911851992, 1288593099, -1468688942, -30913810, -1457169864, 1777881245, -1023857853, 673480136, 296707683, -1357242429, 1128891982, -443173541, -1583655606, 2128771887, 1832324610, -1110612938, -1668289757, -1186148623, -1160154791, -1000534122, 1004420194, 60375650, 846216349, -1593532901, -622209753, 787931386, 750656998, -292464561, 555582556, 198744244, 958077222, -1058742838, -37476263, 1387482063, 951169597, -402915923, -1726665318, 1883781355, 1658181791, 1297942582, 661374004, -1472313698, 115168309, -620570396, -1226238559, 1706867812, -594045534, -275383505, -986723484, -1007768804, -1958620892, -1666405312, -251862087, -2034219957, -998409933, -1099867779, -670878539, 335784102, -110914108, 862517114, 295161105, -2088759981, -37958886, -2026997081, -1498865109, 1136658251, 982592272, -533095050, 364776187, 715465087, -1445586588, 1898168612, -1984079354, 708836796, 468914343, 1218665768, -499913052, -1304093744, -1181233769, -900695587, -1276808037, 1746393301, 355914644, -1702649159, 776239991, 253273886, -261832850, 1055530440, 1022277998, 185673866, 564976759, -1944377144, 1737275767, -977480280, 96937403, 1574485237, 415735413, -1780304209, -2133120928, -1762075980, -84419887, -70942773, 2130034380, 19742978, -1636584661, -2063105268, -1351394370, -889045068, 1154992222, 1509485043, 301334141, -522321119, -1970191130, -360840691, -1021119420, 1189803598, 1865340337, -1825929687, -43362241, 772658117, 532348110, 185543576, 1282081071, -1017380486, -560009962, -566896695, 871150800, 295328809, -1014231440, 1493347022, -1456583190, -2036382469, 1908515162, 506223686, 112715925, 2035041853, 509242950, 219026191, 319170796, 1447250329, -1881229121, 266220803, -1765492456, 1314981838, 923593487, 643664239, 1557864641, -366698088, 904438975, 603328719, -673447328, -1052075494, 1904817875, 2140244797, -885780090, -1312527852, 669681321, -1250409203, -423657852, 1767782970, 1730154204, -914675972, -588477677, 1762389904, 870845872, -1658495672, -331849169, -1485741622, -367541735, -1487059468, 1868802391, 1911697761, -834347377, 516923548, 372284591, 403920721, 1922682661, -567959312, -1696635395, 18042930, -1910496001, -577406164, -404308106, -916169780, 288331114, 447470550, -169010527, -1064160463, 1237236858, -1582502460, 1714391359, 1856798124, 1754796901, -2042369966, -1529728823, 1635343365, -1211908487, -808676591, 973461322, -1794491533, -1203196772, 1309795545, 1587853607, -1138902729, 674382745, 1717694925, 1743925170, 2008481623, 294748460, 1255345654, -1297342768, 420881773, 1763741096, 769812504, -1231039874, 436960348, -746267297, 1318602185, 960768666, -1640932996, -1438463940, 1982873135, -249707681, 769717125, 1142505824, -1416646680, -1009207389, 2080502854, -1559713002, -1731624852, 1408898148, -483407560, -324531184, 390637825, 282221498, -31972397, 1776481464, -768986039, 1412626062, -1757144231, 558965029, 884810266, -160061135, 2011240051, 1744021533, -1747186872, 1313709018, 606946181, -1840787370, 1555652170, -374381388, 889692031, -1111091917, 1849372615, -1867334632, 1134241754, 827268277, 1121041937, 1224807073, -1705201131, -1430915116, -446461225, 1821701128, -955395501, -209839851, 1310766608, -1338326064, 1754291595, -83764866, 1980988547, 1100975706, -614741509, -1307074881, 401778142, -841188137, 951965790, 892307475, 71169236, -1895595952, 1181092356, 1500976014, 1107489927, -544879266, 1051248725, 1056636179, -18384763, 1383389442, 493489518, -879758552, 201545641, 262886423, 1468908908, -1956020530, 183146614, -225487610, 440511337, 1471768414, 1373963546, 1494774388, 2039588081, -229926534, 35687512, 1715809712, -1981169122, 1606636423, -910779543, 1350823866, 1552453262, 1504934899, -890148507, -1468280016, -1441642451, -151165921, 1254587017, -1235954158, 1372835484, -1756877201, -1809777854, -1430458385, -83139650, -359662832, 1314625394, -1161264689, 1388888038, -1332830229, -876302256, 1832654842, -1511835418, -1585926318, 1038108665, -605597439, 1385405439, 968045019, 1855615592, -264347373, -1111939733, 1007410871, 1022920250, 1276722540, 84146206, 1909183786, 797334952, -1216785692, -466032198, 1264605614, 624), None)
#)
print random.getstate()

numTrials = 100000
numFailures = 0
numXFailures = numZFailures = 0

#percentComplete = 0
#for trial in range(numTrials): 
#	pc = int(100 * trial / float(numTrials))
#	if pc > percentComplete:
#		percentComplete = pc
#		print "%3d%% complete, so far error rate is %f" % (percentComplete, numFailures / float(trial))
#	
#	errors = sampleExRec(errorRates)
#	
#	isXMalignant, isZMalignant = checkMalignantUsingDecoder(errors, corrector)
#	if isXMalignant:
#		numXFailures += 1
#	if isZMalignant:
#		numZFailures += 1
#	if isXMalignant or isZMalignant:
#		numFailures += 1
#		#break
#
#print errorRates
#print numFailures, "failures out of", numTrials, "trials is", 100*numFailures/float(numTrials), "%"
#print "numXFailures=", numXFailures, ", numZFailures=", numZFailures


def findPseudoThresh(pStart=0.0009, pStep=0.0001):
	p1 = 0
	p = pStart
	while p1 <= p:
		errorRates = sim.simulateUtils.ErrorRates()
		errorRates.cnot = p
		errorRates.prepX = errorRates.prepZ = errorRates.measX = errorRates.measZ = 4/15. * errorRates.cnot
		errorRates.rest = 12/15. * errorRates.cnot
		#errorRates.rest = 0
	
		xErrors, zErrors, numFailures = simExRecErrors(numTrials, errorRates, pool)
		p1 = float(numFailures)/numTrials
		print p, p1, xErrors, zErrors, numFailures
		
		p += pStep

if __name__ == "__main__":
	import sys
	
	nSlots = 1	
	if 1 < len(sys.argv):
		nSlots = int(sys.argv[1])
	
	pool = configureMultiProcess(nSlots)
	
	findPseudoThresh(0.0015)
	
#	for i in range(len(golayCode.goodSteaneRandomPreps)):
#		
#		def prepareSteaneRandom_i(*args, **kwargs):
#			return prepareAncillaSteaneRandom(i, *args, **kwargs)
#		
#		prepareAncilla = prepareSteaneRandom_i
#		
#		print 'Steane Random ', i
#	
#		findPseudoThresh(0.0015)
		
	print 'Done'