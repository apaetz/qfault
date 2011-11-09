#
# Some basic routines for preparing encoded stabilizer states.
#
# Adam Paetznick (9/12/2011)
from util.counterUtils import locXprep, locZprep, loccnot, locrest, \
	propagateAllErrors


def permuteCircuit(circuit, rounds):
	return [circuit[r] for r in rounds]

def ancillaZPrep(schedule, roundPermutation=None, name='0', reverseQubitOrder=True):
	"""
	Returns the circuit for preparing encoded |0> based on the given CNOT schedule (given as a list of rounds).
	By default, the qubit ordering of the circuit is reversed.  For example, if there are n qubits, then 
	n <-> 0, (n-1) <-> 1, etc.  This is usually desired because the location counting code uses a
	descending (little endian) bit ordering.  That is, an error on qubit k is represented by a 1 on bit
	k.  However, it is more natural to express operators, stabilizer generators in particular, in ascending
	(big endian) bit ordering.  For example XII is an X on qubit 0, not on qubit 2.
	To disable qubit reversal, set reverseQubitOrder=False.
	"""
	if None != roundPermutation:
		cnots = permuteCircuit(schedule, roundPermutation)
	else:
		cnots = schedule
	
	locations = []
	def addLoc(loc):
		locations.append(loc)

	#flatcnots = reduce(lambda x,y: x+y, cnots)
	flatcnots = [cnot for round in cnots for cnot in round]
	allQubits = set(qubit for cnot in flatcnots for qubit in cnot)
	
	numQubits = max(allQubits)
	bitLookup = [i for i in range(numQubits+1)]
	if reverseQubitOrder:
		bitLookup = list(reversed(bitLookup))
	
	# A qubit that is a target of any CNOT in the prep circuit should
	# be prepared as |0>.  All other qubits should be prepared as |+>
	prepAsZ  = set([c[1] for c in flatcnots])
	prepAsX = allQubits - prepAsZ

	# the initial qubit preparations
	for i in prepAsX:
		addLoc(locXprep(name, bitLookup[i]))
	for i in prepAsZ:
		addLoc(locZprep(name, bitLookup[i]))

	touched = set()
	
	for rnd in xrange(len(cnots)):			# simulate each round sequentially
		# Rests will only be added for qubits which have already been touched.
		resting = set(touched)
		
		for src, tgt in cnots[rnd]:
			cnotQubits = set([src, tgt])
			addLoc(loccnot(name, bitLookup[src], name, bitLookup[tgt]))
			resting -= cnotQubits
			touched.update(cnotQubits)
			
		# now add the rest locations
		for r in resting:
			addLoc(locrest(name, bitLookup[r]))
			
	#propagateAllErrors(locations)
	return locations